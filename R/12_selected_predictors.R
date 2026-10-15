# Separate observed-data Gaussian models avoid treating filled outcomes as predictors.
selected_predictor_em <- function(data, traits, predictor_traits, max_iter, tol, trait_order, min_iter, cor_tol) {
  validate_traits(data, traits)
  selected <- select_observed_phenotypes(data, traits)
  original <- selected$data
  out <- original
  if (is.character(predictor_traits)) {
    mapping <- stats::setNames(rep(list(predictor_traits), length(traits)), traits)
  } else {
    mapping <- predictor_traits
    if (!is.list(mapping) || is.null(names(mapping)) || anyDuplicated(names(mapping)) ||
        !setequal(names(mapping), traits)) stop("Predictor list must be named for every target trait")
  }
  models <- logs <- histories <- list()
  for (target in traits) {
    predictors <- mapping[[target]]
    if (!is.character(predictors) || anyNA(predictors) || anyDuplicated(predictors)) {
      stop("Each predictor selection must be a unique character vector")
    }
    predictors <- setdiff(predictors, target)
    columns <- c(target, predictors)
    validate_traits(original, columns)
    # Fit each target with only the explicitly selected observed predictors.
    fit <- sequential_em(original, columns, max_iter=max_iter, tol=tol,
      trait_order=trait_order, min_iter=min_iter, cor_tol=cor_tol)
    mu <- fit$distribution_mean[columns]
    S <- fit$distribution_covariance[columns, columns, drop=FALSE]
    rows <- which(is.na(original[[target]]))
    for (i in rows) {
      observed <- predictors[!is.na(unlist(original[i, predictors, drop=FALSE], use.names=FALSE))]
      value <- mu[target]
      if (length(observed)) value <- value + drop(
        (as.numeric(original[i, observed, drop=TRUE]) - mu[observed]) %*%
          solve(S[observed, observed, drop=FALSE], S[observed, target, drop=FALSE]))
      out[[target]][i] <- value
    }
    log <- fit$diagnostics[fit$diagnostics$trait==target, , drop=FALSE]
    log$predictor_traits <- if (length(predictors)) paste(predictors, collapse=",") else "intercept_only"
    log$sequence_order <- match(target, traits)
    log$n_observed <- sum(!is.na(original[[target]]))
    log$n_missing <- length(rows)
    logs[[target]] <- log
    models[[target]] <- fit
    h <- fit$correlation_history
    if (nrow(h)) { h$target <- target; histories[[target]] <- h }
  }
  # Correlation column meanings differ per target; keep histories in target_models.
  c(list(data=out, missing=is.na(as.matrix(original[traits])),
    diagnostics=do.call(rbind,logs), imputation_order=traits,
    observed_counts=colSums(!is.na(as.matrix(original[traits]))),
    method="selected_predictor_gaussian_em", predictor_traits=mapping,
    target_models=models, correlation_history=NULL, correlation_pairs=NULL), selected$audit)
}
