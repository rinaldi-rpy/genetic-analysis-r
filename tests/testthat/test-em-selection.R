test_that("all-selected-traits missing rows are excluded without imputing them", {
  d <- data.frame(ID = letters[1:6], BW6 = c(1, 2, 4, 3, NA, NA), BL6 = c(3, 7, 6, 8, NA, NA))
  r <- sequential_em(method="sequential", d, c("BW6", "BL6"))
  expect_equal(r$data, d[1:4, ])
  expect_equal(r$excluded_row_indices, 5:6)
  expect_equal(r$excluded_ids, c("e", "f"))
  expect_equal(r$n_excluded, 2L)
  expect_true(all(r$diagnostics$n_missing == 0))
  expect_true(all(r$diagnostics$iterations == 0))
  expect_equal(nrow(r$missing), 4L)
})
test_that("more observed values determine chain order while model trait order stays supplied", {
  set.seed(456)
  x <- rnorm(40)
  y <- 2 * x + rnorm(40)
  d <- data.frame(ID = paste0("a", 1:40), BW6 = y, BL6 = x)
  d$BW6[1:12] <- NA
  d$BL6[40] <- NA
  r <- sequential_em(method="sequential", d, c("BW6", "BL6"))
  expect_equal(r$imputation_order, c("BL6", "BW6"))
  expect_equal(r$diagnostics$n_observed, c(39L, 28L))
  expect_equal(r$diagnostics$predictor_traits, c("intercept_only", "BL6"))
  expect_equal(colnames(r$missing), c("BW6", "BL6"))
  expect_equal(names(r$data), names(d))
  expect_false(anyNA(r$data[c("BW6", "BL6")]))
  expect_equal(r$data$BW6[13:40], d$BW6[13:40])
  expect_equal(sequential_em(method="sequential", d, c("BW6", "BL6"), trait_order = "supplied")$imputation_order, c("BW6", "BL6"))
})
test_that("multi-trait all-missing exclusion, ties, and insufficient records are explicit", {
  d <- data.frame(a = c(1, 2, 4, 3, NA), b = c(5, 4, 8, 7, NA), c = c(7, 3, 6, 8, NA))
  r <- sequential_em(method="sequential", d, c("b", "a", "c"))
  expect_equal(r$imputation_order, c("b", "a", "c"))
  expect_equal(r$n_excluded, 1L)
  expect_error(sequential_em(method="sequential", data.frame(a = rep(NA_real_, 5), b = rep(NA_real_, 5)), c("a", "b")), "No observations")
  expect_error(sequential_em(method="sequential", data.frame(a = c(1, 2, NA, NA), b = c(1, 3, NA, NA)), c("a", "b")), "Fewer than 3")
})
test_that("pipeline excludes absent measurements but keeps full-pedigree predictions", {
  exe <- Sys.getenv("BLUPF90_PLUS", unset = "")
  skip_if(!nzchar(exe), "Set BLUPF90_PLUS for native test")
  x <- make_example(2L)
  d <- x$d
  d[1:4, x$traits] <- NA_real_
  expect_warning_free <- function(impute) {
    run_genetic_analysis(x$p, d, x$traits,
      fixed = ~Sex, backend = "blupf90", impute = impute, blupf90_control = list(executable = exe)
    )
  }
  r <- expect_warning_free(TRUE)
  expect_equal(nrow(r$completed_data), 116L)
  expect_equal(r$diagnostics$n_excluded_all_missing_traits, 4L)
  expect_setequal(r$diagnostics$excluded_phenotype_ids, d$ID[1:4])
  expect_false(r$diagnostics$single_imputation)
  expect_true(all(d$ID[1:4] %in% r$ebv$ID))
  expect_equal(nrow(r$ebv), nrow(x$p) * 2L)
  r2 <- expect_warning_free(FALSE)
  expect_equal(r2$diagnostics$n_excluded_all_missing_traits, 4L)
})
