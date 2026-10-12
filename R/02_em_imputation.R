#' Complete ordered traits using sequential Gaussian EM
#' @param data Wide phenotype data frame.
#' @param traits Unique numeric trait column names in temporal order.
#' @param max_iter Maximum EM iterations per trait.
#' @param tol Relative parameter convergence tolerance.
#' @param predictors Use all earlier traits or just the previous trait.
#' @return List with completed data, original missing mask, and diagnostics.
#' @details The first trait uses an intercept-only normal model. Each later
#'   trait uses completed earlier traits as fixed predictors. Observed values
#'   remain unchanged. Conditional variances enter the EM second moments.
#'   This assumes ignorable missingness and linear Gaussian relationships.
#' @export
sequential_em <- function(data, traits, max_iter = 500L, tol = 1e-8,
                          predictors = c("previous", "all")) {
  predictors <- match.arg(predictors)
  validate_traits(data, traits)
  if (length(max_iter) != 1L || !is.numeric(max_iter) || !is.finite(max_iter) ||
    max_iter < 1 || max_iter != floor(max_iter) || length(tol) != 1L ||
    !is.numeric(tol) || !is.finite(tol) || tol <= 0) {
    stop("Invalid EM controls")
  }
  out <- data
  mask <- is.na(as.matrix(data[traits]))
  logs <- vector("list", length(traits))
  for (j in seq_along(traits)) {
    y <- data[[traits[j]]]
    missing <- is.na(y)
    prev <- if (j == 1L) integer() else if (predictors == "all") seq_len(j - 1L) else j - 1L
    X <- cbind(Intercept = 1, as.matrix(out[traits[prev]]))
    observed <- !missing
    q <- ncol(X)
    if (sum(observed) <= q) stop("Insufficient observations for trait ", traits[j])
    if (qr(X[observed, , drop = FALSE])$rank < q) stop("Predictor rank deficiency for trait ", traits[j])
    b <- qr.solve(X[observed, , drop = FALSE], y[observed])
    v <- mean((y[observed] - drop(X[observed, , drop = FALSE] %*% b))^2)
    if (any(missing) && (!is.finite(v) || v <= .Machine$double.eps)) stop("Zero residual variance for trait ", traits[j])
    converged <- !any(missing)
    delta <- 0
    iter <- 0L
    ll <- numeric()
    if (any(missing)) {
      for (iter in seq_len(max_iter)) {
        mu <- drop(X %*% b)
        ey <- y
        ey[missing] <- mu[missing]
        bn <- qr.solve(X, ey)
        vn <- (sum((ey - drop(X %*% bn))^2) + sum(missing) * v) / length(y)
        delta <- max(abs(c(bn - b, vn - v)) / (1 + abs(c(b, v))))
        b <- bn
        v <- vn
        ll <- c(ll, sum(dnorm(y[observed], drop(X[observed, , drop = FALSE] %*% b), sqrt(v), log = TRUE)))
        if (delta < tol) {
          converged <- TRUE
          break
        }
      }
    }
    if (!converged) stop("EM failed to converge for trait ", traits[j], "; increase max_iter")
    out[[traits[j]]][missing] <- drop(X[missing, , drop = FALSE] %*% b)
    logs[[j]] <- data.frame(
      trait = traits[j], n_missing = sum(missing), iterations = iter,
      converged = converged, delta = delta, residual_variance = v
    )
    attr(logs[[j]], "log_likelihood") <- ll
  }
  list(
    data = out, missing = mask, diagnostics = do.call(rbind, logs),
    log_likelihood = lapply(logs, attr, "log_likelihood")
  )
}

validate_traits <- function(data, traits) {
  if (!is.data.frame(data) || nrow(data) < 3L || anyDuplicated(names(data)) || !is.character(traits) ||
    !length(traits) || anyNA(traits) || anyDuplicated(traits) ||
    !all(traits %in% names(data))) {
    stop("Provide data and unique existing trait names")
  }
  for (t in traits) {
    if (!is.numeric(data[[t]]) || any(!is.finite(data[[t]]) & !is.na(data[[t]]))) {
      stop("Traits must be numeric and finite or NA")
    }
  }
}
