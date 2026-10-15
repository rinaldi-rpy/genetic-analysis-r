test_that("accuracy uses prediction error variance and inbreeding-adjusted prior variance", {
  ebv <- data.frame(ID = c("b", "a"), trait = "BW6", EBV = c(3, -1), SE = c(3, 3))
  G <- matrix(25, 1, 1, dimnames = list("BW6", "BW6"))
  F <- data.frame(ID = c("a", "b"), F = c(0, .25))
  out <- add_ebv_accuracy(ebv, G, F)
  expect_equal(out$PEV, c(9, 9))
  expect_equal(out$Reliability, c(1 - 9 / 31.25, .64))
  expect_equal(out$Accuracy, c(sqrt(1 - 9 / 31.25), .8))
  expect_equal(out$Accuracy_percent, 100 * out$Accuracy)
  expect_equal(out[c("ID", "trait", "EBV", "SE")], ebv)
})
test_that("accuracy maps multiple traits by name rather than row order", {
  G <- diag(c(4, 25))
  dimnames(G) <- list(c("YW", "BW"), c("YW", "BW"))
  ebv <- data.frame(ID = c("a", "a"), trait = c("BW", "YW"), EBV = 0, SE = c(3, 1))
  out <- add_ebv_accuracy(ebv, G, data.frame(ID = "a", F = 0))
  expect_equal(out$Reliability, c(.64, .75))
})
test_that("accuracy handles no information, perfect prediction, boundaries and rounding", {
  G <- matrix(4, 1, 1, dimnames = list("w", "w"))
  e <- data.frame(ID = c("a", "b", "c"), trait = "w", EBV = 0, SE = c(2, 0, sqrt(4 * (1 + 1e-5))))
  F <- data.frame(ID = c("a", "b", "c"), F = 0)
  r <- add_ebv_accuracy(e, G, F)
  expect_equal(r$Accuracy, c(0, 1, 0))
  G[1, 1] <- 0
  expect_true(all(is.na(add_ebv_accuracy(e, G, F)$Accuracy)))
  G[1, 1] <- 4
  e$SE[1] <- 3
  expect_warning(r <- add_ebv_accuracy(e, G, F), "exceeds")
  expect_true(is.na(r$Accuracy[1]))
  expect_error(add_ebv_accuracy(e, G, data.frame(ID = "z", F = 0)), "ID")
  expect_error(add_ebv_accuracy(e, unname(G), F), "trait")
})
test_that("real model outputs include accuracy for BLUPF90+", {
  x <- make_example(1L)
  exe <- Sys.getenv("BLUPF90_PLUS", unset = "")
  if (nzchar(exe)) {
    r <- run_genetic_analysis(x$p, x$d, x$traits,
      backend = "blupf90", impute = FALSE,
      blupf90_control = list(executable = exe)
    )
    expect_true(all(c("PEV", "Reliability", "Accuracy", "Accuracy_percent") %in% names(r$ebv)))
    expect_equal(r$ebv$PEV, r$ebv$SE^2)
    expect_true(all(r$ebv$Accuracy >= 0 & r$ebv$Accuracy <= 1))
  }
})
