#' Prepare a sparse BLUPF90+ animal-model job for any trait count
#' @param pedigree Pedigree with ID, Sire, Dam and optional Birth_Date.
#' @param data Completed wide phenotype data, one row per ID.
#' @param traits Ordered numeric trait columns; no fixed 1-4 trait assumption.
#' @param fixed One-sided fixed formula, fitted separately within each trait.
#' @param random Named list of additional pedigree_effect/iid_effect terms.
#' @param n_iters Maximum AI-REML rounds.
#' @param tol Native REML convergence criterion.
#' @param em_rounds EM-REML warmup rounds before AI-REML.
#' @param max_pairs Maximum pedigree relationship cache entries.
#' @param max_traits Configurable resource limit on trait count.
#' @return A blupf90_job with numeric data, ID map, sparse inverse triplets,
#'   generated parameter lines, and audit metadata. Does not run or download.
#' @export
prepare_blupf90 <- function(pedigree, data, traits, fixed = ~1, n_iters = 1000L, tol = 1e-10,
                            em_rounds = 10L, max_pairs = 2000000L, max_traits = 32L, random = list()) {
  validate_traits(data, traits)
  validate_blup_controls(n_iters, tol, em_rounds, max_traits)
  k <- length(traits)
  if (k > max_traits) stop("Trait count exceeds max_traits; increase the resource limit deliberately")
  validate_phenotype_ids(data)
  rel <- pedigree_precision(pedigree, max_pairs = max_pairs)
  if (!all(as.character(data$ID) %in% rel$pedigree$ID)) stop("Every phenotype ID must occur in pedigree")
  X <- animal_fixed_design(data, traits, fixed)
  Y <- as.matrix(data[traits])
  if (anyNA(Y)) stop("Complete traits with sequential_em before BLUPF90 preparation")
  scales <- apply(Y, 2L, stats::var)
  if (any(scales <= 0)) stop("Traits must have positive observed variance")
  # Nonzero off-diagonals are essential: BLUPF90 fixes initial zeros to zero.
  S <- outer(sqrt(scales), sqrt(scales)) * (diag(k) * 0.95 + matrix(0.05, k, k))
  G <- S * 0.3
  E <- S * 0.7
  id_map <- data.frame(code = seq_len(nrow(rel$pedigree)), ID = rel$pedigree$ID)
  random_effects <- build_random_effects(random, data, rel, X, S)
  numeric_data <- cbind(Y, X, do.call(cbind,lapply(random_effects, `[[`, "codes")))
  colnames(numeric_data)[k+ncol(X)+seq_along(random_effects)] <- names(random_effects)
  missing <- -999L
  while (any(Y == missing)) missing <- missing - 1L
  triplet <- methods::as(rel$A_inverse, "TsparseMatrix")
  inverse <- data.frame(i = triplet@i + 1L, j = triplet@j + 1L, x = triplet@x)
  inverse <- as.data.frame(inverse[inverse$i <= inverse$j, c("i", "j", "x")])
  fixed_positions <- k + seq_len(ncol(X))
  effect <- ncol(X) + 1L
  effect_lines <- c(
    vapply(fixed_positions, function(pos) {
      paste(c(rep(pos, k), 1L, "cov"), collapse = " ")
    }, character(1)),
    vapply(seq_along(random_effects), function(j) {
      paste(c(rep(k+ncol(X)+j,k), nrow(random_effects[[j]]$map), "cross"),collapse=" ")
    }, character(1))
  )
  job <- structure(list(
    data = numeric_data, pedigree = rel$pedigree,
    inbreeding = rel$inbreeding, id_map = id_map, inverse = inverse, traits = traits,
    n_traits = k, animal_effect = effect, random_effects = random_effects,
    n_effects = ncol(X)+length(random_effects), fixed_columns = colnames(X), G_start = G, E_start = E,
    missing_code = missing, n_iters = n_iters, tol = tol, em_rounds = em_rounds,
    effect_lines = effect_lines, parameter_map = do.call(rbind, lapply(random_effects, function(r) {
      map <- blup_parameter_map(k,r$effect); map <- map[map$component=="G",]; map$component <- r$name; map
    }))
  ), class = "blupf90_job")
  residual_map <- blup_parameter_map(k,effect)
  job$parameter_map <- rbind(job$parameter_map,residual_map[residual_map$component=="R",])
  job <- configure_direct_maternal(job)
  job$parameter_lines <- blup_parameter_lines(job, G, E, estimate = TRUE)
  job
}

validate_phenotype_ids <- function(data) {
  if (!"ID" %in% names(data) || anyNA(data$ID) || anyDuplicated(as.character(data$ID)) ||
    any(trimws(as.character(data$ID)) == "")) {
    stop("Phenotype ID must be present, unique and nonmissing")
  }
  if (any(startsWith(names(data), ".ga_"))) stop("Reserved .ga_ column names in data")
}

animal_fixed_design <- function(data, traits, fixed) {
  if (!inherits(fixed, "formula") || length(fixed) != 2L) stop("fixed must be a one-sided formula")
  if (any(all.vars(fixed) %in% c(traits, "ID"))) stop("Fixed predictors cannot be traits or ID")
  mf <- model.frame(fixed, data, na.action = na.pass)
  if (anyNA(mf)) stop("Missing fixed predictors")
  X <- model.matrix(fixed, mf)
  if (!ncol(X) || any(!is.finite(X)) || qr(X)$rank < ncol(X) || nrow(X) <= ncol(X)) {
    stop("Fixed design is rank deficient or saturated")
  }
  X
}

validate_blup_controls <- function(n_iters, tol, em_rounds, max_traits) {
  for (x in list(n_iters = n_iters, em_rounds = em_rounds, max_traits = max_traits)) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x != floor(x)) {
      stop("Invalid BLUPF90 integer controls")
    }
  }
  if (n_iters < 2L || max_traits < 1L || em_rounds >= n_iters ||
    !is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("Invalid BLUPF90 controls")
  }
}

blup_parameter_map <- function(k, effect) {
  pairs <- which(upper.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  # Labels, not guessed AI positions, identify all estimated parameters.
  rbind(
    data.frame(
      component = "G", trait1 = pairs[, 1], trait2 = pairs[, 2],
      label = paste0("G_", effect, "_", effect, "_", pairs[, 1], "_", pairs[, 2])
    ),
    data.frame(
      component = "R", trait1 = pairs[, 1], trait2 = pairs[, 2],
      label = paste0("R_", pairs[, 1], "_", pairs[, 2])
    )
  )
}

matrix_lines <- function(M) apply(M, 1L, function(row) paste(sprintf("%.17g", row), collapse = " "))

blup_parameter_lines <- function(job, G, E, estimate, random_covariances = NULL, direct_maternal_covariance = NULL) {
  k <- job$n_traits
  lines <- c(
    "DATAFILE", "phenotypes.txt", "NUMBER_OF_TRAITS", as.character(k),
    "NUMBER_OF_EFFECTS", as.character(job$n_effects), "OBSERVATION(S)",
    paste(seq_len(k), collapse = " "), "WEIGHT(S)", "", "EFFECTS:", job$effect_lines,
    "RANDOM_RESIDUAL VALUES", matrix_lines(E)
  )
  for (name in names(job$random_effects)) {
    if(!is.null(job$direct_maternal_effect) && name==job$direct_maternal_effect) next
    r <- job$random_effects[[name]]
    covariance <- if (name=="animal") G else if (!is.null(random_covariances)) random_covariances[[name]] else r$start
    if (is.null(covariance)) stop("Missing estimated random covariance: ",name)
    group <- as.character(r$effect)
    if(name=="animal" && !is.null(job$direct_maternal_effect)) {
      maternal <- job$random_effects[[job$direct_maternal_effect]]
      M <- if(is.null(random_covariances)) maternal$start else random_covariances[[maternal$name]]
      C <- if(is.null(direct_maternal_covariance)) {
        if(!estimate) stop("Supply converged direct_maternal_covariance for correlated prediction")
        job$direct_maternal_start
      } else direct_maternal_covariance
      covariance <- direct_maternal_block(G,M,C)
      group <- paste(r$effect,maternal$effect)
    }
    lines <- c(lines,"RANDOM_GROUP",group,"RANDOM_TYPE",
      if(r$type=="pedigree") "user_file" else "diagonal", "FILE",
      if(r$type=="pedigree") "ainverse.txt" else "", "(CO)VARIANCES",matrix_lines(covariance))
  }
  lines <- c(lines,paste("OPTION missing",job$missing_code),"OPTION sol se","OPTION solv_method FSPAK")
  if (!estimate) {
    return(c(lines, "OPTION method BLUP"))
  }
  lines <- c(
    lines, "OPTION method VCE", paste("OPTION maxrounds", job$n_iters),
    paste("OPTION conv_crit", sprintf("%.17g", job$tol)), paste("OPTION EM-REML", job$em_rounds)
  )
  g <- function(i, j) paste0("G_", job$animal_effect, "_", job$animal_effect, "_", i, "_", j)
  for (i in seq_len(k)) {
    denominator <- paste(c(vapply(job$random_effects,function(r) paste0("G_",r$effect,"_",r$effect,"_",i,"_",i),character(1)),paste0("R_",i,"_",i)),collapse="+")
    if(!is.null(job$direct_maternal_effect)) denominator <- paste0(denominator,"+",dm_label(job,i,i))
    expr <- paste0(g(i, i), "/(",denominator,")")
    lines <- c(lines, paste("OPTION se_covar_function", paste0("h2_", i), expr))
    for (name in setdiff(names(job$random_effects),"animal")) {
      r <- job$random_effects[[name]]
      expr <- paste0("G_",r$effect,"_",r$effect,"_",i,"_",i,"/(",denominator,")")
      lines <- c(lines,paste("OPTION se_covar_function",paste0("vprop_",r$effect,"_",i),expr))
    }
    if (i < k) {
      for (j in (i + 1L):k) {
        expr <- paste0(g(i, j), "/sqrt(", g(i, i), "*", g(j, j), ")")
        lines <- c(lines, paste("OPTION se_covar_function", paste0("rg_", i, "_", j), expr))
      }
    }
  }
  c(lines,direct_maternal_se_lines(job))
}

#' Write a prepared BLUPF90+ job into a fresh directory
#' @param job A blupf90_job from prepare_blupf90.
#' @param directory New directory to create; existing directories are rejected.
#' @return Normalized directory path, invisibly.
#' @export
write_blupf90_job <- function(job, directory) {
  if (!inherits(job, "blupf90_job")) stop("job must come from prepare_blupf90")
  if (!is.character(directory) || length(directory) != 1L || is.na(directory) ||
    file.exists(directory)) {
    stop("Choose a new job directory; existing paths are never overwritten")
  }
  if (!dir.create(directory, recursive = TRUE, mode = "0700")) stop("Cannot create job directory")
  directory <- normalizePath(directory, mustWork = TRUE)
  write.table(job$data, file.path(directory, "phenotypes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(job$inverse, file.path(directory, "ainverse.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  audit_pedigree <- data.frame(
    ID = job$id_map$code,
    Sire = match(job$pedigree$Sire, job$id_map$ID, nomatch = 0L),
    Dam = match(job$pedigree$Dam, job$id_map$ID, nomatch = 0L)
  )
  write.table(audit_pedigree, file.path(directory, "pedigree.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  for (name in names(job$random_effects)) write.csv(job$random_effects[[name]]$map, file.path(directory,paste0("effect_",name,"_map.csv")),row.names=FALSE)
  write.csv(job$id_map, file.path(directory, "id_map.csv"), row.names = FALSE)
  write.csv(job$parameter_map, file.path(directory, "parameter_map.csv"), row.names = FALSE)
  write.csv(job$inbreeding, file.path(directory, "inbreeding.csv"), row.names = FALSE)
  writeLines(job$parameter_lines, file.path(directory, "analysis.par"))
  invisible(directory)
}
