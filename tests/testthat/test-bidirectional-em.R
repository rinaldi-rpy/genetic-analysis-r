test_that("bidirectional EM fills both directions and waits at least ten cycles", {
  set.seed(18)
  x <- rnorm(150)
  y <- 5 + 2 * x + rnorm(150, .0, .5)
  d <- data.frame(ID = seq_len(150), x = x, y = y)
  d$x[1:20] <- NA
  d$y[21:50] <- NA
  d[149:150, c("x", "y")] <- NA_real_
  r <- sequential_em(d, c("x", "y"), method = "bidirectional")
  expect_equal(r$n_excluded, 2L)
  expect_false(anyNA(r$data[c("x", "y")]))
  expect_equal(r$data$x[21:148], x[21:148])
  expect_equal(r$data$y[c(1:20, 51:148)], y[c(1:20, 51:148)])
  expect_true(all(r$diagnostics$iterations >= 10))
  expect_equal(r$diagnostics$predictor_traits, c("y", "x"))
  expect_true(tail(r$correlation_history$corr_diff, 1) <= 1e-6)
  expect_true(all(diff(r$correlation_history$log_likelihood) >= -1e-7))
  mu <- r$distribution_mean
  S <- r$distribution_covariance
  expect_equal(r$data$x[1:20], mu["x"] + S["x", "y"] / S["y", "y"] * (d$y[1:20] - mu["y"]), ignore_attr = TRUE)
  expect_equal(r$data$y[21:50], mu["y"] + S["y", "x"] / S["x", "x"] * (d$x[21:50] - mu["x"]), ignore_attr = TRUE)
})
test_that("no imputation skips loops and unidentified covariance is rejected", {
  d <- data.frame(x = 1:5, y = c(2, 5, 1, 4, 3))
  r <- sequential_em(d, c("x", "y"))
  expect_true(all(r$diagnostics$iterations == 0L))
  expect_equal(r$data, d)
  expect_error(sequential_em(d, c("x", "y"), max_iter = 9), "controls")
  d <- data.frame(x = c(1, 2, 3, NA, NA, NA), y = c(NA, NA, NA, 4, 5, 6))
  expect_error(sequential_em(d, c("x", "y")), "overlap")
})
test_that("joint EM scales to three traits with non-monotone missingness", {
  set.seed(65)
  x <- rnorm(100)
  y <- x + rnorm(100)
  z <- y + rnorm(100)
  d <- data.frame(x = x, y = y, z = z)
  d$x[1:8] <- NA
  d$y[9:18] <- NA
  d$z[19:32] <- NA
  r <- sequential_em(d, c("x", "y", "z"))
  expect_false(anyNA(r$data))
  expect_equal(nrow(r$correlation_pairs), 3L)
  expect_true(all(r$diagnostics$iterations >= 10))
})

test_that('joint EM agrees with analytic normal ML for one-sided missing outcomes', {
  set.seed(44)
  x <- rnorm(150); y <- 4+1.5*x+rnorm(150,sd=.4)
  d <- data.frame(x=x,y=y)
  d$y[order(x,decreasing=TRUE)[1:30]] <- NA
  observed <- !is.na(d$y)
  fit <- stats::lm.fit(cbind(1,x[observed]),y[observed])
  b <- fit$coefficients[2]; a <- fit$coefficients[1]
  sx <- mean((x-mean(x))^2); residual <- mean(fit$residuals^2)
  expected <- matrix(c(sx,b*sx,b*sx,b*b*sx+residual),2)
  r <- sequential_em(d,c('x','y'))
  expect_equal(unname(r$distribution_covariance),expected,tolerance=1e-6)
  expect_equal(unname(r$distribution_mean),unname(c(mean(x),a+b*mean(x))),tolerance=1e-6)
  expect_equal(r$data$y[!observed],a+b*x[!observed],ignore_attr=TRUE,tolerance=1e-6)
})
