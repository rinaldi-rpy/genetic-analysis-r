#' Run genetic analysis with BLUPF90+
#' @param pedigree Pedigree data frame.
#' @param data Wide phenotype data frame with one row per ID.
#' @param traits Ordered numeric trait names for model/output. EM defaults to
#'   joint bidirectional EM with most-observed-first processing priority.
#' @param fixed One-sided fixed-effect formula, with trait-specific coefficients.
#' @param impute Complete partially missing traits by sequential EM. Rows missing
#'   every selected trait are excluded from fitting with either value of impute.
#' @param em_control Named sequential_em controls.
#' @param model_control Deprecated; use blupf90_control.
#' @param max_animals Pedigree resource guard; default 100000.
#' @param keep_model Retain backend model/audit metadata.
#' @param backend Only blupf90 is supported.
#' @param blupf90_control Named fit_blupf90_model controls: executable, job_root,
#'   keep_files, timeout, threads, n_iters, tol, em_rounds, max_pairs, max_traits.
#' @param predictor_traits Optional numeric auxiliary EM predictor columns.
#' @param random Named list of additional pedigree_effect or iid_effect terms.
#' @return A genetic_analysis list with consistent EBV, heritability, correlation,
#'   inbreeding, covariance, phenotype and convergence tables.
#' @export
run_genetic_analysis <- function(
  pedigree, data, traits, fixed = ~1, impute = TRUE,
  em_control = list(), model_control = list(), max_animals = NULL, keep_model = FALSE,
  backend = "blupf90", blupf90_control = list(),
  predictor_traits = NULL, random = list()
) {
  if (!identical(backend, "blupf90")) stop("This package wraps BLUPF90+ only; backend must be blupf90")
  validate_traits(data, traits)
  validate_phenotype_ids(data)
  for (flag in list(impute, keep_model)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("Invalid logical pipeline flag")
  }
  check_control(em_control, c("max_iter", "tol", "predictors", "trait_order", "method", "min_iter", "cor_tol"))
  check_control(model_control, c("n_iters", "tol", "max_records"))
  check_control(blupf90_control, c(
    "executable", "job_root", "keep_files", "timeout", "threads",
    "n_iters", "tol", "em_rounds", "max_pairs", "max_traits"
  ))
  if (length(model_control)) stop("Use blupf90_control for backend blupf90")
  if (is.null(max_animals)) max_animals <- 100000L
  if (!is.data.frame(pedigree) || !is.numeric(max_animals) || length(max_animals) != 1L ||
    !is.finite(max_animals) || max_animals < 1 || nrow(pedigree) > max_animals) {
    stop("Pedigree exceeds max_animals resource guard or is invalid")
  }
  if (backend == "blupf90") {
    executable <- if (is.null(blupf90_control$executable)) Sys.getenv("BLUPF90_PLUS", unset = "blupf90+") else blupf90_control$executable
    blupf90_control$executable <- resolve_blup_executable(executable)
  }
  p <- sort_pedigree(pedigree)
  if (nrow(p) > max_animals) stop("Expanded pedigree exceeds max_animals resource guard")
  added_founders <- setdiff(p$ID, trimws(as.character(pedigree$ID)))
  if (!all(as.character(data$ID) %in% p$ID)) stop("Every phenotype ID must occur in the pedigree")
  if (impute) {
    completed <- do.call(sequential_em, c(list(data = data, traits = traits, predictor_traits = predictor_traits), em_control))
    if (any(completed$missing)) warning("Single imputation: SEs exclude imputation uncertainty; genetic estimates may be biased", call. = FALSE)
  } else {
    selected <- select_observed_phenotypes(data, traits)
    retained <- selected$data
    if (anyNA(retained[traits])) stop("Missing traits with impute=FALSE")
    completed <- c(list(
      data = retained, missing = is.na(as.matrix(retained[traits])),
      diagnostics = NULL, imputation_order = character(),
      observed_counts = colSums(!is.na(as.matrix(retained[traits])))
    ), selected$audit)
  }
  fitted <- do.call(fit_blupf90_model, c(list(pedigree = p, data = completed$data,
    traits = traits, fixed = fixed, random = random), blupf90_control))
  inbreeding <- fitted$inbreeding
  diagnostics <- fitted$diagnostics
  diagnostics$EBV_accuracy_method <- "sqrt(1-PEV/((1+F)*genetic_variance))"
  diagnostics$accuracy_includes_imputation_uncertainty <- FALSE
  diagnostics$single_imputation <- any(completed$missing)
  diagnostics$SE_includes_imputation_uncertainty <- FALSE
  diagnostics$n_traits <- length(traits)
  diagnostics$n_input_phenotypes <- completed$n_input
  diagnostics$n_used_phenotypes <- completed$n_retained
  diagnostics$n_excluded_all_missing_traits <- completed$n_excluded
  diagnostics$excluded_phenotype_ids <- completed$excluded_ids
  diagnostics$imputation_order <- completed$imputation_order
  diagnostics$imputation_method <- if (impute) completed$method else "none"
  diagnostics$added_founder_count <- length(added_founders)
  diagnostics$added_founder_ids <- added_founders
  structure(list(
    ebv = fitted$ebv, random_solutions = fitted$random_solutions,
    random_covariances = fitted$random_covariances,
    direct_maternal_covariance = fitted$direct_maternal_covariance,
    joint_genetic_covariance = fitted$joint_genetic_covariance,
    direct_maternal_parameters = fitted$direct_maternal_parameters,
    variance_ratios = fitted$variance_ratios,
    heritability = fitted$heritability,
    genetic_correlations = fitted$genetic_correlations, inbreeding = inbreeding,
    genetic_covariance = fitted$genetic_covariance, residual_covariance = fitted$residual_covariance,
    pedigree = p, completed_data = completed$data, imputation = completed, diagnostics = diagnostics,
    model = if (keep_model) fitted$model else NULL
  ), class = c("genetic_analysis", "list"))
}

check_control <- function(control, allowed) {
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
    anyDuplicated(names(control)) || any(!names(control) %in% allowed)))) {
    stop("Invalid or unknown control arguments")
  }
}

#' @export
print.genetic_analysis <- function(x, ...) {
  cat("Animal model [", x$diagnostics$backend, "]: ", nrow(x$pedigree), " pedigree animals; ",
    nrow(x$heritability), " traits\n",
    sep = ""
  )
  print(x$heritability, row.names = FALSE)
  invisible(x)
}
