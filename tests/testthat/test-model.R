test_that("long-format covariance mapping and SEs match analytic expressions", {
  cmp <- list(
    G = matrix(c(4, 1, 1, 9), 2), E = diag(c(6, 3)),
    gi = matrix(c(1L, 2L, 2L, 3L), 2), ei = matrix(c(4L, 5L, 5L, 6L), 2), count = 6L
  )
  S <- diag(6) * 0.04
  S[1, 4] <- S[4, 1] <- 0.01
  par <- genetic_parameters(cmp, S, c("a", "b"))
  expect_equal(par$h2$h2, c(.4, .75))
  expect_equal(par$h2$SE[1], sqrt(.06^2 * .04 + .04^2 * .04 - 2 * .06 * .04 * .01))
  expect_equal(par$rg$rg, 1 / 6)
  expect_equal(par$rg$SE, sqrt(.04 * ((1 / 6)^2 + (1 / 48)^2 + (1 / 108)^2)))
  cmp$G[1, 1] <- 0
  expect_true(is.na(genetic_parameters(cmp, S, c("a", "b"))$rg$SE))
})

test_that("full-pedigree prediction has correct unrecorded animal behavior", {
  A <- diag(3)
  dimnames(A) <- list(c("a", "b", "c"), c("a", "b", "c"))
  d <- data.frame(ID = c("a", "b"), .ga_y = c(-2, 2))
  r <- full_pedigree_blup(d, matrix(1, 2, 1), A, matrix(2), matrix(1), "w")
  expect_equal(r$EBV, c(-4 / 3, 4 / 3, 0))
  expect_equal(r$SE[3], sqrt(2))
})

test_that("real sommer single-trait pipeline covers all pedigree animals", {
  skip_if_not_installed("sommer")
  x <- make_example(1L)
  r <- run_genetic_analysis(backend = "sommer", x$p, x$d, x$traits, fixed = ~Sex, impute = FALSE, keep_model = TRUE)
  expect_s3_class(r, "genetic_analysis")
  expect_equal(nrow(r$ebv), nrow(x$p))
  expect_true(all(is.finite(r$ebv$SE)))
  expect_equal(nrow(r$genetic_correlations), 0L)
  expect_true(r$diagnostics$converged)
  expect_true(r$heritability$h2 >= 0 && r$heritability$h2 <= 1)
  expect_output(print(r), "Animal model")
  expect_error(run_genetic_analysis(backend = "sommer", x$p, x$d, x$traits, model_control = list(bad = 1)), "control")
})

test_that("real sommer multi-trait EM pipeline is finite and auditable", {
  skip_if_not_installed("sommer")
  x <- make_example()
  x$d$w2[c(3, 9, 21)] <- NA
  expect_warning(r <- run_genetic_analysis(backend = "sommer", x$p, x$d, x$traits, fixed = ~Sex), "imputation")
  expect_equal(nrow(r$ebv), nrow(x$p) * 2)
  expect_equal(nrow(r$heritability), 2L)
  expect_equal(nrow(r$genetic_correlations), 1L)
  expect_true(all(is.finite(r$heritability$SE)))
  expect_true(all(is.finite(r$genetic_correlations$SE)))
  expect_true(all(abs(r$genetic_correlations$rg) <= 1 + 1e-6))
  expect_true(r$diagnostics$single_imputation)
  expect_null(r$model)
  expect_equal(sum(r$imputation$missing), 3)
  expect_error(run_genetic_analysis(backend = "sommer", x$p, x$d, x$traits, impute = FALSE), "Missing")
})

test_that("four trait complete-data model produces all six correlations", {
  skip_if_not_installed("sommer")
  x <- make_example(4L)
  r <- run_genetic_analysis(
    backend = "sommer", x$p, x$d, x$traits,
    impute = FALSE,
    model_control = list(n_iters = 200L, tol = 1e-5)
  )
  expect_equal(nrow(r$genetic_correlations), 6L)
  expect_equal(nrow(r$ebv), 4 * nrow(x$p))
})

test_that("unsupported covariance layouts and invalid estimates fail explicitly", {
  expect_error(covariance_components(list(sigma = list(matrix(1))), 2), "Unsupported")
  bad <- list(sigma = setNames(list(matrix(1)), "unknown"))
  expect_error(covariance_components(bad, 1), "Unexpected")
  bad$sigma <- setNames(list(matrix(1)), "Q:.ga_animal")
  expect_error(covariance_components(bad, 2), "Cannot map")
  bad$sigma <- setNames(list(matrix(1)), "T1:.ga_animal")
  expect_error(covariance_components(bad, 2), "Incomplete")
  bad$sigma <- setNames(list(matrix(1), matrix(-1)), c("u:.ga_animal", "units"))
  expect_error(covariance_components(bad, 1), "admissible")
  cmp <- list(G = matrix(1), E = matrix(1), gi = matrix(1L), ei = matrix(2L), count = 2L)
  expect_error(genetic_parameters(cmp, NULL, "a"), "Unsupported")
  expect_error(genetic_parameters(cmp, -diag(2), "a"), "Invalid")
})
