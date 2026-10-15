#' Fit a sparse pedigree animal model using BLUPF90+
#' @param pedigree Pedigree data frame.
#' @param data Completed wide phenotype data.
#' @param traits Ordered numeric trait names.
#' @param fixed One-sided fixed formula.
#' @param random Named additional random effects.
#' @param executable BLUPF90+ executable path/name; BLUPF90_PLUS environment
#'   variable is used when set. Must be a trusted server-controlled executable.
#' @param job_root Parent directory for an isolated temporary job.
#' @param keep_files Retain logs and native files for audit.
#' @param timeout Maximum seconds per native execution.
#' @param threads Positive OpenMP/MKL thread count; default one avoids worker oversubscription.
#' @param n_iters Maximum variance-estimation rounds.
#' @param tol REML convergence criterion.
#' @param em_rounds EM warmup rounds.
#' @param max_pairs Inbreeding relationship-cache limit.
#' @param max_traits Configurable trait-count resource limit.
#' @return Backend-neutral result list, including full-pedigree EBVs and SEs.
#' @details Performs VCE then a final BLUP pass using the converged covariances.
#'   External processes use isolated working directories without changing R's
#'   working directory. Downloads never occur implicitly during analysis.
#' @export
fit_blupf90_model <- function(
  pedigree, data, traits, fixed = ~1,
  executable = Sys.getenv("BLUPF90_PLUS", unset = "blupf90+"), job_root = tempdir(),
  keep_files = FALSE, timeout = 600, threads = 1L, n_iters = 1000L, tol = 1e-10,
  em_rounds = 10L, max_pairs = 2000000L, max_traits = 32L, random = list()
) {
  executable <- resolve_blup_executable(executable)
  if (!is.logical(keep_files) || length(keep_files) != 1L || is.na(keep_files)) stop("Invalid keep_files")
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout <= 0 ||
    !is.numeric(threads) || length(threads) != 1L || !is.finite(threads) || threads < 1 || threads != floor(threads)) {
    stop("Invalid timeout or threads")
  }
  if (!is.character(job_root) || length(job_root) != 1L || is.na(job_root) || !dir.exists(job_root)) {
    stop("job_root must be an existing directory")
  }
  job <- prepare_blupf90(pedigree, data, traits, fixed, n_iters, tol, em_rounds, max_pairs, max_traits, random)
  directory <- tempfile("blupf90-job-", tmpdir = normalizePath(job_root, mustWork = TRUE))
  write_blupf90_job(job, directory)
  if (!keep_files) on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  run_native_blup(executable, "analysis.par", directory, "vce", timeout, threads)
  log <- file.path(directory, "vce.stdout.log")
  estimates <- parse_blup_estimates(job, log)
  writeLines(
    blup_parameter_lines(job, estimates$G, estimates$E, estimate = FALSE, random_covariances = estimates$random_covariances,
      direct_maternal_covariance = estimates$direct_maternal_covariance),
    file.path(directory, "prediction.par")
  )
  if (file.exists(file.path(directory, "solutions"))) unlink(file.path(directory, "solutions"))
  run_native_blup(executable, "prediction.par", directory, "prediction", timeout, threads)
  result <- read_blupf90_results(job, log, file.path(directory, "solutions"))
  result$diagnostics$executable <- executable
  result$diagnostics$threads <- threads
  result$diagnostics$job_directory <- if (keep_files) directory else NULL
  result$model <- if (keep_files) {
    list(
      job_directory = directory, id_map = job$id_map,
      parameter_map = job$parameter_map
    )
  } else {
    NULL
  }
  result
}

resolve_blup_executable <- function(executable) {
  if (!is.character(executable) || length(executable) != 1L || is.na(executable) || !nzchar(executable)) {
    stop("Provide a BLUPF90+ executable")
  }
  if (identical(executable, "blupf90+")) {
    installed <- file.path(tools::R_user_dir("geneticAnalysisR", "data"), "bin", if (.Platform$OS.type=="windows") "blupf90+.exe" else "blupf90+")
    if (!nzchar(Sys.which(executable)) && file.exists(installed)) executable <- installed
  }
  if (!file.exists(executable)) executable <- unname(Sys.which(executable))
  if (!nzchar(executable) || !file.exists(executable) || dir.exists(executable) || (.Platform$OS.type!="windows" && file.access(executable, 1L) != 0L)) {
    stop("BLUPF90+ executable not found; call install_blupf90() or supply its trusted path", call. = FALSE)
  }
  normalizePath(executable, mustWork = TRUE)
}

run_native_blup <- function(executable, parameter, directory, prefix, timeout, threads) {
  env <- c("current",
    OMP_NUM_THREADS = as.character(threads), MKL_NUM_THREADS = as.character(threads),
    OPENBLAS_NUM_THREADS = as.character(threads)
  )
  res <- tryCatch(processx::run(executable,
    args = parameter, wd = directory,
    stdout = file.path(directory, paste0(prefix, ".stdout.log")),
    stderr = file.path(directory, paste0(prefix, ".stderr.log")), timeout = timeout,
    env = env, error_on_status = FALSE
  ), error = function(e) {
    stop("BLUPF90+ process could not complete (timeout, executable or runtime-library error)", call. = FALSE)
  })
  if (res$status != 0L) stop("BLUPF90+ failed; inspect retained job logs and runtime libraries", call. = FALSE)
  invisible(res$status)
}
