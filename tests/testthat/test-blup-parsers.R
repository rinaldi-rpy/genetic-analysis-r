make_native_fixture <- function(k = 3L) {
  set.seed(20)
  p <- data.frame(ID = c("a", "b", "c", "d", "e"), Sire = c(NA, NA, "a", "a", "c"), Dam = c(NA, NA, "b", "b", "d"))
  d <- data.frame(ID = p$ID)
  for (i in seq_len(k)) d[[paste0("t", i)]] <- rnorm(5)
  job <- prepare_blupf90(p, d, paste0("t", seq_len(k)))
  G <- diag(k) * 3 + matrix(.5, k, k)
  E <- diag(k) * 2 + matrix(.2, k, k)
  lines <- c(
    "In round 20 convergence = 1.0D-14", "delta convergence = 1.0D-09",
    "Final Estimates", paste("Genetic variance(s) for effect", job$animal_effect),
    matrix_lines(G), "Residual variance(s)", matrix_lines(E)
  )
  for (i in seq_len(k)) {
    v <- G[i, i] / (G[i, i] + E[i, i])
    lines <- c(
      lines, paste0("h2_", i, " - Function : test"), paste(" Mean :", v),
      paste(" Sample Mean :", v), " Sample SD : 3.0D-02"
    )
    if (i < k) {
      for (j in (i + 1):k) {
        v <- G[i, j] / sqrt(G[i, i] * G[j, j])
        lines <- c(
          lines, paste0("rg_", i, "_", j, " - Function : test"), paste(" Mean :", v),
          paste(" Sample Mean :", v), " Sample SD : 4.0D-02"
        )
      }
    }
  }
  log <- tempfile()
  sol <- tempfile()
  writeLines(lines, log)
  s <- expand.grid(trait = seq_len(k), effect = job$animal_effect, level = seq_len(nrow(p)))
  s$solution <- s$trait / 10
  s$se <- .2
  write.table(s, sol, row.names = FALSE, col.names = TRUE, quote = FALSE)
  list(job = job, log = log, sol = sol, lines = lines, G = G, E = E)
}

test_that("native outputs preserve any trait count and full pedigree ID mapping", {
  for (k in c(1L, 2L, 3L, 6L)) {
    x <- make_native_fixture(k)
    on.exit(unlink(c(x$log, x$sol)), add = TRUE)
    r <- read_blupf90_results(x$job, x$log, x$sol)
    expect_equal(nrow(r$ebv), 5 * k)
    expect_equal(r$ebv$ID, rep(x$job$id_map$ID, k))
    expect_equal(nrow(r$genetic_correlations), k * (k - 1) / 2)
    expect_equal(r$heritability$SE, rep(.03, k))
    expect_equal(r$genetic_correlations$SE, rep(.04, k * (k - 1) / 2))
    expect_equal(unname(r$genetic_covariance), x$G)
  }
})

test_that("native parser rejects nonconvergence and incomplete outputs", {
  x <- make_native_fixture()
  on.exit(unlink(c(x$log, x$sol)), add = TRUE)
  writeLines(sub("1.0D-14", "1.0D-02", sub("1.0D-09", "1.0D-02", x$lines)), x$log)
  expect_error(read_blupf90_results(x$job, x$log, x$sol), "converge")
  writeLines(x$lines[!grepl("Sample SD", x$lines)], x$log)
  expect_error(read_blupf90_results(x$job, x$log, x$sol), "SE output")
  writeLines(x$lines, x$log)
  s <- readLines(x$sol)
  writeLines(s[-2], x$sol)
  expect_error(read_blupf90_results(x$job, x$log, x$sol), "every pedigree")
  expect_error(read_blupf90_results(x$job, "/does/not/exist", x$sol), "missing")
})

test_that("sparse precision equals dense inverse for randomized inbred pedigrees", {
  for (seed in 1:4) {
    set.seed(seed)
    n <- 60
    p <- data.frame(ID = as.character(seq_len(n)), Sire = NA_character_, Dam = NA_character_)
    for (i in 5:n) {
      parent <- sample(seq_len(i - 1), 2)
      p$Sire[i] <- parent[1]
      p$Dam[i] <- parent[2]
    }
    a <- pedigree_precision(p)
    b <- pedigree_relationship(p)
    expect_equal(a$inbreeding, b$inbreeding, ignore_attr = TRUE)
    expect_equal(as.matrix(a$A_inverse), b$A_inverse, tolerance = 1e-10)
  }
})

test_that("job writer uses isolated paths and refuses existing directories", {
  x <- make_native_fixture(2)
  on.exit(unlink(c(x$log, x$sol)), add = TRUE)
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  old <- getwd()
  write_blupf90_job(x$job, path)
  expect_identical(getwd(), old)
  expect_true(all(file.exists(file.path(path, c("analysis.par", "ainverse.txt", "phenotypes.txt", "id_map.csv")))))
  expect_error(write_blupf90_job(x$job, path), "existing")
  expect_error(prepare_blupf90(x$job$pedigree, data.frame(ID = letters[1:5], w = 1:5), "w", max_traits = 0), "controls")
  expect_error(install_blupf90(tempfile(), url = "http://evil.invalid/test",allow_emulation=TRUE), "official")
})
