test_that("three-column offspring records gain unique founders before children", {
  p <- data.frame(ID = c("C2", "C1"), Sire = c("S", "S"), Dam = c("D2", "D1"))
  s <- sort_pedigree(p)
  expect_setequal(s$ID, c("C1", "C2", "S", "D1", "D2"))
  expect_equal(nrow(s), 5)
  expect_true(all(is.na(s$Sire[s$ID %in% c("S", "D1", "D2")])))
  expect_true(all(is.na(s$Dam[s$ID %in% c("S", "D1", "D2")])))
  expect_true(all(is.na(s$Birth_Date)))
  expect_true(all(match(c("S", "D1", "D2"), s$ID) < match("C1", s$ID)))
  expect_equal(sort_pedigree(s), s)
  expect_equal(sort_pedigree(p[2:1, ])$ID, s$ID)
})
test_that("added founder metadata stay missing and dated offspring define ordering", {
  p <- data.frame(
    ID = c("C1", "C2"), Sire = c("SA", "SB"), Dam = c("DA", "DB"),
    Birth_Date = as.Date(c("2020-10-10", "2020-10-08")), Group = factor(c("a", "b")), Score = c(1, 2)
  )
  s <- sort_pedigree(p)
  expect_equal(s$ID, c("DB", "SB", "DA", "SA", "C2", "C1"))
  expect_s3_class(s$Birth_Date, "Date")
  expect_s3_class(s$Group, "factor")
  expect_true(all(is.na(s$Group[1:4])))
  expect_equal(s$Score[match(p$ID, s$ID)], p$Score)
})
test_that("inferred founders preserve dense and sparse relationships and resource guards", {
  p <- data.frame(ID = c("c", "d", "e"), Sire = c("a", "a", "c"), Dam = c("b", "b", "d"))
  dense <- pedigree_relationship(p)
  sparse <- pedigree_precision(p)
  expect_equal(nrow(dense$A), 5)
  expect_equal(dense$A["c", "d"], .5)
  expect_equal(dense$inbreeding$F_percent[dense$inbreeding$ID == "e"], 25)
  expect_equal(as.matrix(sparse$A_inverse), dense$A_inverse, tolerance = 1e-10)
  expect_error(pedigree_relationship(p, max_animals = 4), "guard")
  expect_error(sort_pedigree(data.frame(ID = c("a", "b"), Sire = c("b", "a"), Dam = c("missing", NA))), "cycle")
})
test_that("native BW6 analysis accepts an offspring-only three-column pedigree", {
  exe <- Sys.getenv("BLUPF90_PLUS", unset = "")
  skip_if(!nzchar(exe), "Set BLUPF90_PLUS for native integration test")
  x <- make_example(1L)
  p <- x$p[x$p$ID %in% x$d$ID, c("ID", "Sire", "Dam")]
  d <- x$d
  names(d)[names(d) == "w1"] <- "BW6"
  d$Age_6Month <- seq(170, 190, length.out = nrow(d))
  expect_error(run_genetic_analysis(p, d, "BW6",
    backend = "blupf90", max_animals = nrow(p),
    blupf90_control = list(executable = exe)
  ), "guard")
  r <- run_genetic_analysis(p, d,
    traits = "BW6", fixed = ~ Age_6Month + Sex, backend = "blupf90",
    impute = FALSE, blupf90_control = list(executable = exe)
  )
  expect_equal(nrow(r$pedigree), 142)
  expect_equal(nrow(r$ebv), 142)
  expect_equal(r$diagnostics$added_founder_count, 22)
  expect_setequal(r$diagnostics$added_founder_ids, setdiff(x$p$ID, p$ID))
  expect_true(r$diagnostics$converged)
})
