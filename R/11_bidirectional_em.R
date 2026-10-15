bidirectional_gaussian_em <- function(data, traits, max_iter, tol, trait_order, min_iter, cor_tol) {
  selected <- select_observed_phenotypes(data, traits)
  out <- selected$data
  mask <- is.na(as.matrix(out[traits]))
  counts <- colSums(!mask)
  chain <- if (trait_order == "most_observed") traits[order(-counts, seq_along(traits))] else traits
  Y <- as.matrix(out[chain])
  missing <- is.na(Y)
  k <- length(chain)
  pairs <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)
  pair_map <- data.frame(
    column = if (nrow(pairs)) paste0("r_T", pairs[, 1], "_T", pairs[, 2]) else character(),
    trait1 = chain[pairs[, 1]], trait2 = chain[pairs[, 2]]
  )
  iterations <- 0L
  delta <- corr_diff <- 0
  history <- list()
  likelihood <- numeric()
  if (!any(missing)) {
    mu <- colMeans(Y)
    S <- crossprod(sweep(Y, 2L, mu, "-")) / nrow(Y)
    EY <- Y
    estimated_r <- gaussian_correlation(S)
    completed_r <- gaussian_correlation(S)
  } else {
    centers <- colMeans(Y, na.rm = TRUE)
    scales <- sqrt(vapply(seq_len(k), function(j) mean((Y[!missing[, j], j] - centers[j])^2), numeric(1)))
    if (any(counts < 3L) || any(!is.finite(scales)) || any(scales <= 0)) {
      stop("Insufficient observations or zero variance for bidirectional EM")
    }
    overlaps <- crossprod(!missing)
    if (any(overlaps[upper.tri(overlaps)] < 3L)) {
      stop("Bidirectional EM requires at least 3 overlapping observed records per trait pair")
    }
    Yz <- sweep(sweep(Y, 2L, centers, "-"), 2L, scales, "/")
    groups <- split(seq_len(nrow(Yz)), apply(missing, 1L, paste, collapse = ""))
    mu <- rep(0, k)
    S <- stats::cor(Yz, use = "pairwise.complete.obs")
    if (any(!is.finite(S))) stop("Unidentifiable trait covariance on overlapping records")
    # Shrink only the initializer toward identity when pairwise covariances are not PD.
    while (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < 1e-4) S <- (S + diag(k)) / 2
    old_r <- gaussian_correlation(S)
    old_completed_r <- stats::cor(gaussian_expectation(Yz, mu, S, groups)$mean)
    converged <- FALSE
    for (iterations in seq_len(max_iter)) {
      old <- c(mu, as.vector(S))
      for (direction in list(groups, rev(groups))) {
        moments <- gaussian_expectation(Yz, mu, S, direction)
        mu <- colMeans(moments$mean)
        centered <- sweep(moments$mean, 2L, mu, "-")
        S <- (crossprod(centered) + moments$conditional_sum) / nrow(Yz)
        S <- (S + t(S)) / 2
        if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) <= 1e-12) {
          stop("Singular trait covariance during bidirectional EM; review redundant traits")
        }
      }
      EY <- gaussian_expectation(Yz, mu, S, groups)$mean
      estimated_r <- gaussian_correlation(S)
      completed_r <- stats::cor(EY)
      corr_diff <- max(
        abs(estimated_r[upper.tri(S)] - old_r[upper.tri(S)]),
        abs(completed_r[upper.tri(S)] - old_completed_r[upper.tri(S)])
      )
      delta <- max(abs(c(mu, as.vector(S)) - old) / (1 + abs(old)))
      ll <- gaussian_observed_loglik(Yz, mu, S, groups) - sum(colSums(!missing) * log(scales))
      h <- data.frame(iteration = iterations, corr_diff = corr_diff, parameter_delta = delta, log_likelihood = ll)
      if (nrow(pair_map)) h[pair_map$column] <- as.list(estimated_r[pairs])
      history[[iterations]] <- h
      likelihood <- c(likelihood, ll)
      if (iterations >= min_iter && corr_diff <= cor_tol && delta <= tol) {
        converged <- TRUE
        break
      }
      old_r <- estimated_r
      old_completed_r <- completed_r
    }
    if (!converged) stop("Bidirectional EM failed to converge; increase max_iter or review data")
    mu <- centers + scales * mu
    S <- S * outer(scales, scales)
    # Final conditional values correspond to the converged joint parameters.
    EY <- sweep(sweep(EY, 2L, scales, "*"), 2L, centers, "+")
  }
  names(mu) <- chain
  dimnames(S) <- list(chain, chain)
  dimnames(estimated_r) <- dimnames(completed_r) <- list(chain, chain)
  for (j in seq_len(k)) out[[chain[j]]][missing[, j]] <- EY[missing[, j], j]
  precision <- if (any(missing)) solve(S) else NULL
  logs <- lapply(seq_len(k), function(j) {
    data.frame(
      sequence_order = j, trait = chain[j], n_observed = counts[chain[j]],
      predictor_traits = if (any(missing[, j])) paste(setdiff(chain, chain[j]), collapse = ",") else "not_needed",
      n_missing = sum(missing[, j]), iterations = iterations, converged = TRUE, delta = delta,
      corr_diff = corr_diff, residual_variance = if (is.null(precision)) S[j, j] else 1 / precision[j, j]
    )
  })
  c(list(
    data = out, missing = mask, diagnostics = do.call(rbind, logs),
    imputation_order = chain, observed_counts = counts,
    log_likelihood = stats::setNames(rep(list(likelihood), k), chain),
    correlation_history = if (length(history)) {
      do.call(rbind, history)
    } else {
      data.frame(iteration = integer(), corr_diff = numeric(), parameter_delta = numeric(), log_likelihood = numeric())
    },
    correlation_pairs = pair_map, distribution_mean = mu, distribution_covariance = S,
    estimated_correlation = estimated_r, completed_correlation = completed_r
  ), selected$audit)
}

gaussian_expectation <- function(Y, mu, S, groups) {
  EY <- Y
  extra <- matrix(0, ncol(Y), ncol(Y))
  for (rows in groups) {
    observed <- which(!is.na(Y[rows[1L], ]))
    missing <- which(is.na(Y[rows[1L], ]))
    if (!length(missing)) next
    coefficient <- solve(S[observed, observed, drop = FALSE], S[observed, missing, drop = FALSE])
    EY[rows, missing] <- sweep(
      sweep(Y[rows, observed, drop = FALSE], 2L, mu[observed], "-") %*% coefficient,
      2L, mu[missing], "+"
    )
    conditional <- S[missing, missing, drop = FALSE] - t(S[observed, missing, drop = FALSE]) %*% coefficient
    extra[missing, missing] <- extra[missing, missing, drop = FALSE] + length(rows) * conditional
  }
  list(mean = EY, conditional_sum = extra)
}

gaussian_observed_loglik <- function(Y, mu, S, groups) {
  total <- 0
  for (rows in groups) {
    observed <- which(!is.na(Y[rows[1L], ]))
    C <- chol(S[observed, observed, drop = FALSE])
    centered <- sweep(Y[rows, observed, drop = FALSE], 2L, mu[observed], "-")
    z <- forwardsolve(t(C), t(centered))
    total <- total - .5 * (length(rows) * (length(observed) * log(2 * pi) + 2 * sum(log(diag(C)))) + sum(z^2))
  }
  total
}

gaussian_correlation <- function(S) {
  R <- S / outer(sqrt(pmax(0, diag(S))), sqrt(pmax(0, diag(S))))
  R[!is.finite(R)] <- NA_real_
  diag(R) <- 1
  R
}
