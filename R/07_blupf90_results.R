#' Read and validate native BLUPF90+ results
#' @param job Prepared blupf90_job carrying trait and ID mappings.
#' @param log_file BLUPF90+ variance-estimation stdout log.
#' @param solutions_file Final BLUP solutions with OPTION sol se (five columns).
#' @return Backend-neutral EBV, covariance, h2 and rg tables and diagnostics.
#' @details Genetic parameter SEs are native Monte-Carlo approximations from
#'   se_covar_function, not guessed inverse-AI parameter ordering. Missing or
#'   incompatible output, nonconvergence, and inadmissible matrices are errors.
#' @export
read_blupf90_results <- function(job, log_file, solutions_file) {
  if (!inherits(job, "blupf90_job")) stop("job must come from prepare_blupf90")
  estimates <- parse_blup_estimates(job, log_file)
  ebv <- parse_blup_solutions(job, solutions_file)
  ebv <- add_ebv_accuracy(ebv, estimates$G, job$inbreeding)
  list(
    ebv = ebv,
    random_solutions = lapply(setdiff(names(job$random_effects),"animal"), function(name) {
      r <- job$random_effects[[name]]
      values <- parse_effect_solutions(job,solutions_file,r)
      if(r$type=="pedigree") values <- add_ebv_accuracy(values,estimates$random_covariances[[name]],job$inbreeding)
      values$effect <- name
      values
    }) |> stats::setNames(setdiff(names(job$random_effects),"animal")),
    random_covariances = estimates$random_covariances,
    direct_maternal_covariance = estimates$direct_maternal_covariance,
    joint_genetic_covariance = estimates$joint_genetic_covariance,
    direct_maternal_parameters = estimates$direct_maternal_parameters,
    variance_ratios = estimates$variance_ratios,
    heritability = estimates$heritability,
    genetic_correlations = estimates$genetic_correlations,
    genetic_covariance = estimates$G, residual_covariance = estimates$E,
    variance_components = NULL, model = NULL, pedigree = job$pedigree, inbreeding = job$inbreeding,
    diagnostics = list(
      backend = "blupf90", converged = TRUE, iterations = estimates$iterations,
      direct_maternal_correlated = !is.null(job$direct_maternal_effect),
      phenotypic_variance_convention = "reference_noninbred_offspring_dam_A=0.5",
      parameter_SE_method = "native_monte_carlo", EBV_SE_method = "prediction_error_sd"
    )
  )
}

f90_numeric <- function(text) suppressWarnings(as.numeric(gsub("[dD]", "E", text)))

parse_blup_estimates <- function(job, log_file) {
  if (!file.exists(log_file)) stop("BLUPF90 log file missing")
  lines <- readLines(log_file, warn = FALSE)
  rounds <- grep("In round.*convergence[[:space:]]*=", lines, ignore.case = TRUE)
  if (!length(rounds)) stop("No BLUPF90 REML convergence diagnostics")
  last <- tail(rounds, 1L)
  criterion <- f90_numeric(sub(".*convergence[[:space:]]*=[[:space:]]*([^[:space:]]+).*", "\\1", lines[last], ignore.case = TRUE))
  delta_lines <- grep("delta convergence[[:space:]]*=", lines, ignore.case = TRUE)
  delta <- if (length(delta_lines)) f90_numeric(sub(".*=[[:space:]]*([^[:space:]]+).*", "\\1", lines[tail(delta_lines, 1L)])) else Inf
  if (!is.finite(criterion) || !(criterion < job$tol || (is.finite(delta) && delta < sqrt(job$tol)))) {
    stop("BLUPF90 REML did not converge under the requested criterion")
  }
  iteration <- as.integer(sub(".*In round[[:space:]]+([0-9]+).*", "\\1", lines[last], ignore.case = TRUE))
  final <- grep("Final Estimates", lines, ignore.case = TRUE)
  if (!length(final)) stop("BLUPF90 final estimates missing")
  lines <- lines[seq.int(tail(final, 1L), length(lines))]
  all_headers <- grep("Genetic.*variance.*effect", lines, ignore.case = TRUE)
  random <- read_random_covariances(job,lines,all_headers)
  random_covariances <- random$covariances
  G <- random_covariances$animal
  r_header <- grep("Residual.*variance", lines, ignore.case = TRUE)
  if (length(r_header) != 1L) stop("Unsupported BLUPF90 covariance headers")
  E <- read_native_matrix(lines, r_header, job$n_traits)
  if (!isSymmetric(G, tol = 1e-6) || !isSymmetric(E, tol = 1e-6) ||
    min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) < -1e-7 ||
    min(eigen(E, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("BLUPF90 estimated covariance matrices are not admissible")
  }
  dimnames(G) <- dimnames(E) <- list(job$traits, job$traits)
  h <- vector("list", job$n_traits)
  r <- list()
  ratios <- list()
  total_variance <- Reduce(`+`,lapply(random_covariances,diag))+diag(E)
  if(!is.null(random$cross)) total_variance <- total_variance+diag(random$cross)
  if(any(!is.finite(total_variance)) || any(total_variance<=0)) stop("Nonpositive reference phenotypic variance")
  maternal_parameters <- read_direct_maternal_parameters(job,lines,G,
    if(!is.null(job$direct_maternal_effect)) random_covariances[[job$direct_maternal_effect]] else NULL,random$cross)
  for (i in seq_len(job$n_traits)) {
    value <- G[i, i] / total_variance[i]
    se <- read_native_function_se(lines, paste0("h2_", i), value)
    h[[i]] <- data.frame(trait = job$traits[i], h2 = value, SE = se)
    for(name in setdiff(names(random_covariances),"animal")) {
      value_r <- random_covariances[[name]][i,i]/total_variance[i]
      effect <- job$random_effects[[name]]$effect
      ratios[[length(ratios)+1L]] <- data.frame(effect=name,trait=job$traits[i],
        proportion=value_r, SE=read_native_function_se(lines,paste0("vprop_",effect,"_",i),value_r))
    }
    if (i < job$n_traits) {
      for (j in (i + 1L):job$n_traits) {
        if (G[i, i] <= 0 || G[j, j] <= 0) {
          value <- se <- NA_real_
        } else {
          value <- G[i, j] / sqrt(G[i, i] * G[j, j])
          se <- read_native_function_se(lines, paste0("rg_", i, "_", j), value)
        }
        r[[length(r) + 1L]] <- data.frame(trait1 = job$traits[i], trait2 = job$traits[j], rg = value, SE = se)
      }
    }
  }
  list(
    G = G, E = E, random_covariances=random_covariances,
    direct_maternal_covariance=random$cross, joint_genetic_covariance=random$joint,
    direct_maternal_parameters=maternal_parameters,
    variance_ratios=if(length(ratios)) do.call(rbind,ratios) else data.frame(effect=character(),trait=character(),proportion=numeric(),SE=numeric()), iterations = iteration, heritability = do.call(rbind, h),
    genetic_correlations = if (length(r)) {
      do.call(rbind, r)
    } else {
      data.frame(trait1 = character(), trait2 = character(), rg = numeric(), SE = numeric())
    }
  )
}

read_native_matrix <- function(lines, header, k) {
  values <- numeric()
  if (header == length(lines)) stop("Truncated BLUPF90 covariance matrix")
  for (i in seq.int(header + 1L, length(lines))) {
    text <- trimws(lines[i])
    if (!nzchar(text)) next
    nums <- f90_numeric(strsplit(text, "[[:space:]]+")[[1]])
    if (any(!is.finite(nums))) break
    values <- c(values, nums)
    if (length(values) >= k * k) break
  }
  if (length(values) != k * k) stop("Truncated or unsupported BLUPF90 covariance matrix")
  matrix(values, k, k, byrow = TRUE)
}

read_native_function_se <- function(lines, label, estimate) {
  start <- grep(paste0("^[[:space:]]*", label, "[[:space:]]+-[[:space:]]*Function"), lines, ignore.case = TRUE)
  if (length(start) != 1L) stop("Missing BLUPF90 genetic parameter SE for ", label)
  headers <- grep("^[[:space:]]*[^[:space:]]+[[:space:]]+-[[:space:]]*Function", lines, ignore.case = TRUE)
  next_header <- headers[headers > start]
  end <- if (length(next_header)) next_header[1L] - 1L else length(lines)
  section <- lines[seq.int(start, min(start + 15L, end))]
  sd <- grep("Sample SD[[:space:]]*:", section, ignore.case = TRUE)
  mean <- grep("^[[:space:]]*Mean[[:space:]]*:", section, ignore.case = TRUE)
  if (length(sd) != 1L || length(mean) != 1L) stop("Unsupported BLUPF90 function SE output for ", label)
  se <- f90_numeric(trimws(sub(".*:", "", section[sd])))
  native <- f90_numeric(trimws(sub(".*:", "", section[mean])))
  if (!is.finite(se) || se < 0 || !is.finite(native) || abs(native - estimate) > 1e-3 * (1 + abs(estimate))) {
    stop("Invalid BLUPF90 function SE or trait mapping for ", label)
  }
  se
}

parse_blup_solutions <- function(job, solutions_file) {
  parse_effect_solutions(job,solutions_file,job$random_effects$animal)
}

parse_effect_solutions <- function(job, solutions_file, effect) {
  if (!file.exists(solutions_file)) stop("BLUPF90 solutions file missing")
  s <- utils::read.table(solutions_file, header = FALSE, skip = 1L, colClasses = "character", quote = "", comment.char = "")
  if (ncol(s) != 5L) stop("BLUPF90 solutions require five columns; enable OPTION sol se")
  s[] <- lapply(s, f90_numeric)
  if (any(!is.finite(as.matrix(s)))) stop("Invalid BLUPF90 solution values")
  names(s) <- c("trait", "effect", "level", "EBV", "SE")
  if (any(s$trait != floor(s$trait) | s$effect != floor(s$effect) | s$level != floor(s$level))) {
    stop("Invalid BLUPF90 solution indices")
  }
  s <- s[s$effect == effect$effect, , drop = FALSE]
  m <- nrow(effect$map)
  k <- job$n_traits
  if (any(!s$trait %in% seq_len(k)) || any(!s$level %in% seq_len(m)) || any(s$SE < 0)) {
    stop("BLUPF90 solution indices or SE outside valid range")
  }
  keys <- paste(s$trait, s$level, sep = ":")
  wanted <- paste(rep(seq_len(k), each = m), rep(seq_len(m), k), sep = ":")
  if (anyDuplicated(keys) || length(keys) != length(wanted) || !setequal(keys, wanted)) {
    stop("BLUPF90 solutions do not cover every pedigree ID and trait")
  }
  s <- s[match(wanted, keys), , drop = FALSE]
  data.frame(ID = effect$map$ID[s$level], trait = job$traits[s$trait], EBV = s$EBV, SE = s$SE)
}
