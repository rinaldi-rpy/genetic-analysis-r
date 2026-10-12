#' Download and install BLUPF90+ for Linux x86-64
#' @param directory Installation directory; defaults to the package user-data bin.
#' @param url Official University of Georgia Linux binary HTTPS URL.
#' @param overwrite Replace an existing executable only when TRUE.
#' @param ca_bundle Optional trusted CA bundle path for HTTPS verification.
#' @param install_runtimes Install pinned Intel MKL/OpenMP libraries into a private
#'   runtime directory using Python pip when native libraries are missing.
#' @param check_runtime Run --help to verify executable/runtime compatibility.
#' @return List with executable path, source URL, checksum and runtime check.
#' @details Installation is explicit and never happens inside model fitting.
#'   Downloads require HTTPS on the official host. A verified bundled certificate
#'   chain can fill missing server intermediates; TLS verification stays enabled.
#'   Intel oneAPI runtime libraries may need installation separately. External
#'   BLUPF90+ licensing terms apply; no executable is bundled in this R package.
#' @export
install_blupf90 <- function(
  directory = file.path(tools::R_user_dir("geneticAnalysisR", "data"), "bin"),
  url = "https://nce.ads.uga.edu/html/projects/programs/Linux/64bit/blupf90+",
  overwrite = FALSE, ca_bundle = NULL, check_runtime = TRUE, install_runtimes = TRUE
) {
  if (Sys.info()[["sysname"]] != "Linux" || !Sys.info()[["machine"]] %in% c("x86_64", "amd64")) {
    stop("Automatic installation currently supports Linux x86-64 only")
  }
  if (!is.character(url) || length(url) != 1L || is.na(url) ||
    !grepl("^https://nce[.]ads[.]uga[.]edu/html/projects/programs/Linux/", url)) {
    stop("Use an official BLUPF90 Linux HTTPS binary URL")
  }
  for (flag in list(overwrite, check_runtime, install_runtimes)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("Invalid installation flag")
  }
  if (!is.character(directory) || length(directory) != 1L || is.na(directory) || !nzchar(directory)) {
    stop("Provide an installation directory")
  }
  directory <- path.expand(directory)
  target <- file.path(directory, "blupf90+")
  if (file.exists(target) && !overwrite) stop("BLUPF90+ already exists; use overwrite=TRUE deliberately")
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, mode = "0700")) {
    stop("Cannot create installation directory")
  }
  directory <- normalizePath(directory, mustWork = TRUE)
  target <- file.path(directory, "blupf90+")
  part <- tempfile(".blupf90-download-", tmpdir = directory)
  on.exit(unlink(part), add = TRUE)
  handle <- curl::new_handle(timeout = 120, followlocation = TRUE, ssl_verifypeer = TRUE, ssl_verifyhost = 2L)
  if (!is.null(ca_bundle)) {
    if (!is.character(ca_bundle) || length(ca_bundle) != 1L || !file.exists(ca_bundle)) stop("Invalid CA bundle")
    curl::handle_setopt(handle, cainfo = normalizePath(ca_bundle, mustWork = TRUE))
  }
  tryCatch(curl::curl_download(url, part, quiet = TRUE, handle = handle), error = function(e) {
    if (!is.null(ca_bundle) || !grepl("certificate|SSL", conditionMessage(e), ignore.case = TRUE)) {
      stop("BLUPF90 download failed: ", conditionMessage(e), call. = FALSE)
    }
    bundle <- blup_verified_ca_bundle()
    on.exit(unlink(bundle), add = TRUE)
    curl::handle_setopt(handle, cainfo = bundle)
    curl::curl_download(url, part, quiet = TRUE, handle = handle)
  })
  con <- file(part, "rb")
  magic <- readBin(con, "raw", n = 20L)
  close(con)
  if (length(magic) < 20L || !identical(magic[1:4], as.raw(c(127, 69, 76, 70))) ||
    magic[5] != as.raw(2) || magic[19] != as.raw(62) || magic[20] != as.raw(0)) {
    stop("Download is not a Linux x86-64 ELF executable")
  }
  Sys.chmod(part, mode = "0755")
  runtime <- NULL
  private_runtime <- FALSE
  runtime_dir <- file.path(directory, "runtime")
  if (check_runtime) {
    runtime <- tryCatch(processx::run(part, "--help", timeout = 30, error_on_status = FALSE),
      error = function(e) stop("Downloaded executable cannot start; inspect Linux runtime libraries", call. = FALSE)
    )
    if (runtime$status != 0L && install_runtimes && grepl("shared libraries", runtime$stderr, fixed = TRUE)) {
      install_blup_runtimes(runtime_dir)
      private_runtime <- TRUE
      ld <- paste(c(file.path(runtime_dir, "lib"), Sys.getenv("LD_LIBRARY_PATH")), collapse = ":")
      runtime <- processx::run(part, "--help",
        timeout = 30, error_on_status = FALSE,
        env = c("current", LD_LIBRARY_PATH = ld)
      )
    }
    if (runtime$status != 0L) stop("Downloaded BLUPF90+ needs compatible runtime libraries; install Intel oneAPI runtimes", call. = FALSE)
  }
  checksum_command <- unname(Sys.which("sha256sum"))
  checksum <- if (nzchar(checksum_command)) strsplit(processx::run(checksum_command, part)$stdout, "[[:space:]]+")[[1]][1] else NA_character_
  if (file.exists(target) && overwrite) unlink(target)
  native_target <- if (private_runtime) file.path(directory, "blupf90+.native") else target
  if (!file.rename(part, native_target)) stop("Cannot finalize BLUPF90 installation")
  if (private_runtime) write_blup_runtime_launcher(target)
  result <- list(executable = target, url = url, SHA256 = checksum, runtime_checked = check_runtime, private_runtime = private_runtime)
  writeLines(
    c(paste("URL", url), paste("SHA256", checksum), paste("Runtime checked", check_runtime)),
    file.path(directory, "blupf90-installation.txt")
  )
  result
}

blup_verified_ca_bundle <- function() {
  cert_dir <- system.file("certs", package = "geneticAnalysisR")
  chain <- file.path(cert_dir, "incommon-chain.pem")
  intermediate <- file.path(cert_dir, "incommon-OVG2C.pem")
  roots <- c(
    Sys.getenv("CURL_CA_BUNDLE"), "/etc/ssl/certs/ca-certificates.crt",
    "/etc/pki/tls/certs/ca-bundle.crt", "/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem"
  )
  roots <- roots[nzchar(roots) & file.exists(roots)]
  openssl <- unname(Sys.which("openssl"))
  if (!length(roots) || !nzchar(openssl) || !file.exists(chain)) {
    stop("Missing trusted CA roots/openssl; supply ca_bundle with a valid certificate chain")
  }
  verified <- processx::run(openssl, c("verify", "-CAfile", roots[1], "-untrusted", chain, intermediate),
    timeout = 15, error_on_status = FALSE
  )
  if (verified$status != 0L) stop("Bundled BLUPF90 certificate chain is not trusted by this device")
  bundle <- tempfile("blupf90-ca-", fileext = ".pem")
  writeLines(c(readLines(roots[1], warn = FALSE), readLines(chain, warn = FALSE)), bundle)
  bundle
}

install_blup_runtimes <- function(directory) {
  python <- unname(Sys.which("python3"))
  if (!nzchar(python)) stop("Python3/pip is needed to install private Intel runtime libraries")
  res <- processx::run(python, c(
    "-m", "pip", "install", "--ignore-installed", "--prefix", directory,
    "--no-warn-script-location", "mkl==2025.3.0", "intel-openmp==2025.3.0"
  ),
  timeout = 1800, error_on_status = FALSE
  )
  if (res$status != 0L) stop("Private Intel runtime installation failed; check pip/network availability")
  invisible(directory)
}

write_blup_runtime_launcher <- function(target) {
  writeLines(c(
    "#!/bin/sh",
    'BLUPF90_BINDIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd) || exit 1',
    'LD_LIBRARY_PATH="$BLUPF90_BINDIR/runtime/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"',
    "export LD_LIBRARY_PATH",
    'exec "$BLUPF90_BINDIR/blupf90+.native" "$@"'
  ), target)
  Sys.chmod(target, mode = "0755")
  invisible(target)
}
