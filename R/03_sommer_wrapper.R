#' Reshape wide phenotypes to trait-major long form
#' @param data Wide phenotype data, with one row per ID.
#' @param traits Ordered numeric trait names.
#' @return Long data with .ga_trait, .ga_y, .ga_record, and .ga_animal.
#' @export
reshape_traits <- function(data, traits) {
  validate_traits(data, traits)
  if (!"ID" %in% names(data) || anyNA(data$ID) || anyDuplicated(as.character(data$ID))) {
    stop("Phenotype ID must be present, unique and nonmissing")
  }
  if (any(startsWith(names(data), ".ga_"))) stop("Reserved .ga_ column names in data")
  n <- nrow(data)
  k <- length(traits)
  long <- data[rep(seq_len(n), k), setdiff(names(data), traits), drop = FALSE]
  long$.ga_trait <- factor(rep(paste0("T", seq_len(k)), each = n), levels = paste0("T", seq_len(k)))
  long$.ga_y <- unlist(data[traits], use.names = FALSE)
  long$.ga_record <- factor(rep(paste0("R", seq_len(n)), k), levels = paste0("R", seq_len(n)))
  long$.ga_animal <- factor(as.character(long$ID))
  rownames(long) <- NULL
  long
}

#' Fit an unstructured multi-trait pedigree animal model
#' @param data Completed wide phenotypes (one row per animal).
#' @param traits Ordered trait names.
#' @param A Named additive relationship matrix covering every phenotype ID.
#' @param fixed One-sided formula; effects are fitted separately for each trait.
#' @param n_iters Maximum sommer iterations.
#' @param tol Model likelihood convergence tolerance.
#' @param max_records Limit on number of long records and pedigree coefficients.
#' @return List with raw model, long data, covariance matrices, EBVs, prediction
#'   SEs, accuracy/reliability and genetic parameters.
#' @details Uses mmer's long-format usr structures for genetic and within-animal
#'   residual covariance. Full-pedigree BLUPs and prediction-error SEs are
#'   calculated from fitted covariances by marginal mixed-model equations.
#' @export
fit_animal_model <- function(data, traits, A, fixed = ~1, n_iters = 100L, tol = 1e-6,
                             max_records = 2000L) {
  if (!requireNamespace("sommer", quietly = TRUE)) stop("Install sommer to use backend = sommer")
  long <- reshape_traits(data, traits)
  k <- length(traits)
  n <- nrow(data)
  if (!inherits(fixed, "formula") || length(fixed) != 2L) stop("fixed must be a one-sided formula")
  if (length(n_iters) != 1L || !is.numeric(n_iters) || !is.finite(n_iters) ||
    n_iters < 2 || n_iters != floor(n_iters) || length(tol) != 1L ||
    !is.numeric(tol) || !is.finite(tol) || tol <= 0) {
    stop("Invalid model controls")
  }
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) != ncol(A) ||
    is.null(rownames(A)) || anyDuplicated(rownames(A)) ||
    !identical(rownames(A), colnames(A)) || any(!is.finite(A)) ||
    !isSymmetric(A) || !all(as.character(data$ID) %in% rownames(A))) {
    stop("A must be symmetric, finite, named and cover phenotype IDs")
  }
  if (length(max_records) != 1L || !is.numeric(max_records) || !is.finite(max_records) ||
    max_records < 1 || n * k > max_records || nrow(A) * k > max_records) {
    stop("Model exceeds max_records dense-memory guard")
  }
  tryCatch(chol(A), error = function(e) stop("A must be positive definite"))
  if (anyNA(long$.ga_y)) stop("Complete traits with sequential_em before fitting")
  if (any(all.vars(fixed) %in% c(traits, "ID", names(long)[startsWith(names(long), ".ga_")]))) {
    stop("Fixed predictors cannot be traits, ID, or reserved columns")
  }
  mf <- model.frame(fixed, data, na.action = na.pass)
  if (anyNA(mf)) stop("Missing fixed predictors; supply complete covariates")
  baseX <- model.matrix(fixed, mf)
  if (!ncol(baseX) || any(!is.finite(baseX)) || qr(baseX)$rank < ncol(baseX) || n <= ncol(baseX)) {
    stop("Fixed design is rank deficient or saturated")
  }
  # A model-matrix column is used for each fixed coefficient and trait.
  # This supports transformations/interactions without parsing user text.
  X <- kronecker(diag(k), baseX)
  xnames <- paste0(".ga_X", seq_len(ncol(X)))
  long[xnames] <- as.data.frame(X)
  ff <- reformulate(xnames, response = ".ga_y", intercept = FALSE)
  A_obs <- A[as.character(data$ID), as.character(data$ID), drop = FALSE]
  long$.ga_animal <- factor(as.character(long$ID), levels = as.character(data$ID))
  # mmer evaluates variance expressions in its caller environment.
  vsr <- long_variance_structure
  usr <- sommer::usr
  random <- if (k == 1L) ~ vsr(.ga_animal, Gu = A_obs) else ~ vsr(usr(.ga_trait), .ga_animal, Gu = A_obs)
  residual <- if (k == 1L) ~units else ~ vsr(usr(.ga_trait), .ga_record)
  # Start from positive-definite trait covariances on the response scale.
  initial <- stats::cov(as.matrix(data[traits]))
  pair_values <- unlist(lapply(seq_len(k), function(i) initial[i, seq_len(i)]))
  init <- lapply(c(0.3 * pair_values, 0.7 * pair_values), function(v) matrix(v, 1L, 1L))
  fit <- sommer::mmer(
    fixed = ff, random = random, rcov = residual, data = long, init = init,
    nIters = n_iters, tolParConvLL = tol, getPEV = FALSE, verbose = FALSE, dateWarning = FALSE,
    stepWeight = rep(0.5, n_iters)
  )
  if (!isTRUE(fit$convergence)) stop("sommer model did not converge; increase n_iters or review model")
  components <- covariance_components(fit, k)
  dimnames(components$G) <- dimnames(components$E) <- list(traits, traits)
  ebv <- full_pedigree_blup(long, X, A, components$G, components$E, traits)
  ebv <- add_ebv_accuracy(
    ebv, components$G,
    data.frame(ID = rownames(A), F = diag(A) - 1)
  )
  pars <- genetic_parameters(components, fit$sigmaSE, traits)
  list(
    model = fit, long_data = long, genetic_covariance = components$G,
    residual_covariance = components$E, ebv = ebv, heritability = pars$h2,
    genetic_correlations = pars$rg, variance_components = summary(fit)$varcomp
  )
}

covariance_components <- function(fit, k) {
  G <- E <- matrix(0, k, k)
  gi <- ei <- matrix(0L, k, k)
  nn <- names(fit$sigma)
  if (is.null(nn) || any(vapply(fit$sigma, length, integer(1)) != 1L)) {
    stop("Unsupported sommer covariance output layout")
  }
  for (s in seq_along(nn)) {
    genetic <- grepl(".ga_animal", nn[s], fixed = TRUE)
    residual <- grepl(".ga_record", nn[s], fixed = TRUE) || (k == 1L && nn[s] == "units")
    if (!genetic && !residual) stop("Unexpected sommer component: ", nn[s])
    if (k == 1L) {
      pair <- c(1L, 1L)
    } else {
      labels <- regmatches(nn[s], gregexpr("T[0-9]+", nn[s]))[[1]]
      pair <- as.integer(sub("T", "", labels, fixed = TRUE))
      if (length(pair) == 1L) pair <- rep(pair, 2L)
      if (length(pair) != 2L || any(!pair %in% seq_len(k))) {
        stop("Cannot map sommer trait covariance: ", nn[s])
      }
    }
    i <- pair[1]
    j <- pair[2]
    value <- as.numeric(fit$sigma[[s]])
    if (genetic) {
      G[i, j] <- G[j, i] <- value
      gi[i, j] <- gi[j, i] <- s
    } else {
      E[i, j] <- E[j, i] <- value
      ei[i, j] <- ei[j, i] <- s
    }
  }
  if (any(gi == 0L) || any(ei == 0L)) stop("Incomplete sommer covariance mapping")
  if (min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) < -1e-7 ||
    min(eigen(E, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Estimated covariance matrices are not admissible")
  }
  list(G = G, E = E, gi = gi, ei = ei, count = length(nn))
}

full_pedigree_blup <- function(long, X, A, G, E, traits) {
  k <- length(traits)
  m <- nrow(A)
  n <- nrow(long) / k
  C <- kronecker(G, A)
  selected <- unlist(lapply(seq_len(k), function(j) {
    match(as.character(long$ID[seq_len(n)]), rownames(A)) + (j - 1L) * m
  }))
  B <- C[, selected, drop = FALSE]
  V <- C[selected, selected, drop = FALSE] + kronecker(E, diag(n))
  Vi <- chol2inv(chol(V))
  VX <- Vi %*% X
  P <- Vi - VX %*% solve(crossprod(X, VX), t(VX))
  u <- drop(B %*% P %*% long$.ga_y)
  pev <- diag(C) - rowSums((B %*% P) * B)
  if (any(pev < -1e-6)) stop("Negative prediction error variance")
  data.frame(ID = rep(rownames(A), k), trait = rep(traits, each = m), EBV = u, SE = sqrt(pmax(0, pev)))
}

genetic_parameters <- function(cmp, S, traits) {
  if (!is.matrix(S) || !identical(dim(S), rep(cmp$count, 2L)) || any(!is.finite(S))) {
    stop("Unsupported sommer variance-parameter covariance layout")
  }
  delta_se <- function(gradient) {
    v <- drop(t(gradient) %*% S %*% gradient)
    if (v < -1e-8) stop("Invalid delta-method variance")
    sqrt(max(0, v))
  }
  h <- vector("list", length(traits))
  r <- list()
  for (i in seq_along(traits)) {
    a <- cmp$G[i, i]
    e <- cmp$E[i, i]
    g <- numeric(cmp$count)
    g[cmp$gi[i, i]] <- e / (a + e)^2
    g[cmp$ei[i, i]] <- -a / (a + e)^2
    h[[i]] <- data.frame(trait = traits[i], h2 = a / (a + e), SE = delta_se(g))
    if (i < length(traits)) {
      for (j in (i + 1L):length(traits)) {
        b <- cmp$G[j, j]
        c <- cmp$G[i, j]
        g <- numeric(cmp$count)
        if (a <= 0 || b <= 0) {
          value <- se <- NA_real_
        } else {
          value <- c / sqrt(a * b)
          g[cmp$gi[i, j]] <- 1 / sqrt(a * b)
          g[cmp$gi[i, i]] <- -value / (2 * a)
          g[cmp$gi[j, j]] <- -value / (2 * b)
          se <- delta_se(g)
        }
        r[[length(r) + 1L]] <- data.frame(trait1 = traits[i], trait2 = traits[j], rg = value, SE = se)
      }
    }
  }
  list(h2 = do.call(rbind, h), rg = if (length(r)) {
    do.call(rbind, r)
  } else {
    data.frame(trait1 = character(), trait2 = character(), rg = numeric(), SE = numeric())
  })
}

# mmer treats residual Z entries as covariance kernels. A record incidence
# matrix is rectangular, so expand it to a square within-record kernel.
long_variance_structure <- function(...) {
  call <- match.call()
  call[[1L]] <- quote(sommer::vsr)
  structure <- eval(call, parent.frame())
  if (".ga_record" %in% structure$terms) {
    structure$Z <- Map(function(Z, K) Z %*% K %*% t(Z), structure$Z, structure$K)
  }
  structure
}
