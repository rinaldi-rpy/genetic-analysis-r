test_that("real BLUPF90+ fits scalable trait models on Linux", {
  exe <- Sys.getenv("BLUPF90_PLUS", unset = "")
  skip_if(!nzchar(exe), "Set BLUPF90_PLUS for real executable integration tests")
  for (k in c(1L, 2L, 3L, 5L)) {
    x <- make_example(k)
    r <- run_genetic_analysis(x$p, x$d, x$traits,
      fixed = ~Sex, impute = FALSE,
      backend = "blupf90", blupf90_control = list(executable = exe, timeout = 300)
    )
    expect_equal(nrow(r$ebv), nrow(x$p) * k)
    expect_equal(nrow(r$heritability), k)
    expect_equal(nrow(r$genetic_correlations), k * (k - 1) / 2)
    expect_true(all(is.finite(r$ebv$SE)))
    expect_true(all(is.finite(r$heritability$SE)))
    expect_equal(r$diagnostics$backend, "blupf90")
    expect_true(r$diagnostics$converged)
  }
})

test_that("known-variance native predictions match exact marginal equations", {
  exe <- Sys.getenv("BLUPF90_PLUS", unset = "")
  skip_if(!nzchar(exe), "Set BLUPF90_PLUS for real executable integration tests")
  p <- data.frame(
    ID = letters[1:8], Sire = c(NA, NA, "a", "a", "c", "c", "e", "e"),
    Dam = c(NA, NA, "b", "b", "d", "d", "f", "f")
  )
  d <- data.frame(ID = letters[4:8], w = c(10, 12, 9, 15, 13))
  job <- prepare_blupf90(p, d, "w")
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  write_blupf90_job(job, path)
  G <- matrix(2)
  E <- matrix(3)
  writeLines(blup_parameter_lines(job, G, E, FALSE), file.path(path, "prediction.par"))
  run_native_blup(resolve_blup_executable(exe), "prediction.par", path, "prediction", 120, 1L)
  native <- parse_blup_solutions(job, file.path(path, "solutions"))
  rel <- pedigree_relationship(p)
  long <- reshape_traits(d, "w")
  exact <- full_pedigree_blup(long, matrix(1, 5, 1), rel$A, G, E, "w")
  expect_equal(native$EBV, exact$EBV, tolerance = 1e-6)
  expect_equal(native$SE, exact$SE, tolerance = 1e-6)
})
