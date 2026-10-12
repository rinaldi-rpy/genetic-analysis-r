#' Complete partially observed traits with bidirectional Gaussian EM
#' @param data Wide phenotype data frame.
#' @param traits Selected numeric trait names in model/output order.
#' @param max_iter Maximum iteration cycles.
#' @param tol Relative parameter-change tolerance; also legacy sequential tolerance.
#' @param predictors Predictor chaining for method='sequential' only.
#' @param trait_order Most-observed-first processing or supplied order.
#' @param method Bidirectional joint Gaussian EM (default) or legacy sequential EM.
#' @param min_iter Minimum bidirectional cycles before convergence, at least ten.
#' @param cor_tol Maximum absolute change in pairwise correlations for convergence.
#' @return Completed retained data, missing mask, per-trait diagnostics, excluded
#'   records and iteration/correlation history. Observed values remain unchanged.
#' @details Entirely unmeasured rows are excluded. Bidirectional cycles use two
#'   Gaussian EM updates, processing missingness patterns forward then in reverse.
#'   Conditional means and covariances enter the sufficient statistics; imputed
#'   values are never treated as observed predictors. Conditional predictions use
#'   the other observed traits in the same row. Stopping requires both correlation
#'   and parameter stability after min_iter cycles. Complete data skip EM.
#' @export
sequential_em <- function(
  data, traits, max_iter = 500L, tol = 1e-8,
  predictors = c("previous", "all"), trait_order = c("most_observed", "supplied"),
  method = c("bidirectional", "sequential"), min_iter = 10L, cor_tol = 1e-6
) {
  method <- match.arg(method)
  predictors <- match.arg(predictors)
  trait_order <- match.arg(trait_order)
  if (method == "sequential") {
    result <- sequential_em_legacy(data, traits, max_iter, tol, predictors, trait_order)
    result$method <- "sequential"
    return(result)
  }
  validate_traits(data, traits)
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
    max_iter != floor(max_iter) || !is.numeric(min_iter) || length(min_iter) != 1L ||
    !is.finite(min_iter) || min_iter != floor(min_iter) || min_iter < 10L || max_iter < min_iter ||
    !is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
    !is.numeric(cor_tol) || length(cor_tol) != 1L || !is.finite(cor_tol) || cor_tol <= 0) {
    stop("Invalid bidirectional EM controls; max_iter >= min_iter >= 10 required")
  }
  result <- bidirectional_gaussian_em(data, traits, max_iter, tol, trait_order, min_iter, cor_tol)
  result$method <- "bidirectional"
  result
}

sequential_em_legacy <- function(data, traits, max_iter = 500L, tol = 1e-8,
                                 predictors = c("previous", "all"),
                                 trait_order = c("most_observed", "supplied")) {
  predictors <- match.arg(predictors)
  validate_traits(data, traits)
  trait_order <- match.arg(trait_order)
  if (length(max_iter) != 1L || !is.numeric(max_iter) || !is.finite(max_iter) ||
    max_iter < 1 || max_iter != floor(max_iter) || length(tol) != 1L ||
    !is.numeric(tol) || !is.finite(tol) || tol <= 0) {
    stop("Invalid EM controls")
  }
  selected <- select_observed_phenotypes(data, traits)
  data <- selected$data
  out <- data
  mask <- is.na(as.matrix(data[traits]))
  observed_counts <- colSums(!mask)
  chain <- if (trait_order == "most_observed") traits[order(-observed_counts, seq_along(traits))] else traits
  logs <- vector("list", length(chain))
  for (j in seq_along(chain)) {
    y <- data[[chain[j]]]
    missing <- is.na(y)
    prev <- if (j == 1L) integer() else if (predictors == "all") seq_len(j - 1L) else j - 1L
    X <- cbind(Intercept = 1, as.matrix(out[chain[prev]]))
    observed <- !missing
    q <- ncol(X)
    if (sum(observed) <= q) stop("Insufficient observations for trait ", chain[j])
    if (qr(X[observed, , drop = FALSE])$rank < q) stop("Predictor rank deficiency for trait ", chain[j])
    b <- qr.solve(X[observed, , drop = FALSE], y[observed])
    v <- mean((y[observed] - drop(X[observed, , drop = FALSE] %*% b))^2)
    if (any(missing) && (!is.finite(v) || v <= .Machine$double.eps)) stop("Zero residual variance for trait ", chain[j])
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
    if (!converged) stop("EM failed to converge for trait ", chain[j], "; increase max_iter")
    out[[chain[j]]][missing] <- drop(X[missing, , drop = FALSE] %*% b)
    logs[[j]] <- data.frame(
      sequence_order = j, trait = chain[j], n_observed = sum(observed),
      predictor_traits = if (length(prev)) paste(chain[prev], collapse = ",") else "intercept_only",
      n_missing = sum(missing), iterations = iter,
      converged = converged, delta = delta, residual_variance = v
    )
    attr(logs[[j]], "log_likelihood") <- ll
  }
  c(list(
    data = out, missing = mask, diagnostics = do.call(rbind, logs),
    log_likelihood = lapply(logs, attr, "log_likelihood"),
    imputation_order = chain, observed_counts = observed_counts
  ), selected$audit)
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

# Drop measurement-free rows only; pedigree inclusion is handled independently.
select_observed_phenotypes <- function(data, traits) {
  excluded <- which(rowSums(!is.na(as.matrix(data[traits]))) == 0L)
  retained <- setdiff(seq_len(nrow(data)), excluded)
  if (!length(retained)) stop("No observations remain after excluding all-missing trait rows; EM was not run")
  if (length(retained) < 3L) stop("Fewer than 3 phenotype rows remain after excluding all-missing trait rows")
  list(data = data[retained, , drop = FALSE], audit = list(
    n_input = nrow(data), n_retained = length(retained), n_excluded = length(excluded),
    retained_row_indices = retained, excluded_row_indices = excluded,
    excluded_ids = if ("ID" %in% names(data)) as.character(data$ID[excluded]) else character(),
    excluded_rows = data[excluded, , drop = FALSE]
  ))
}
