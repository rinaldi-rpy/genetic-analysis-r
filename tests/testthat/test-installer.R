test_that("explicit installer checks downloads and never overwrites without consent", {
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  system <- Sys.info()[["sysname"]]
  magic <- switch(system,
    Linux={ m<-as.raw(rep(0,20));m[1:5]<-as.raw(c(127,69,76,70,2));m[19]<-as.raw(62);m },
    Windows={ m<-as.raw(rep(0,512));m[1:2]<-charToRaw("MZ");m[61]<-as.raw(128);m[129:134]<-as.raw(c(80,69,0,0,100,134));m },
    Darwin=as.raw(c(207,250,237,254,7,0,0,1)))
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    writeBin(magic, destfile)
    invisible(destfile)
  }, .package = "curl")
  res <- install_blupf90(path, check_runtime = FALSE,allow_emulation=TRUE)
  expect_true(file.exists(res$executable))
  expect_false(res$runtime_checked)
  expect_equal(nchar(res$SHA256), 64)
  expect_error(install_blupf90(path, check_runtime = FALSE,allow_emulation=TRUE), "already exists")
  expect_error(install_blupf90(tempfile(),allow_emulation=TRUE, overwrite = NA), "flag")
  expect_error(install_blupf90(NULL,allow_emulation=TRUE), "directory")
  expect_error(install_blupf90(tempfile(),allow_emulation=TRUE, ca_bundle = "/missing/ca"), "CA bundle")
})

test_that("installer rejects non-executable responses before any execution", {
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  local_mocked_bindings(curl_download = function(url, destfile, ...) {
    writeLines("<html>not a binary</html>", destfile)
    invisible(destfile)
  }, .package = "curl")
  expect_error(install_blupf90(path, check_runtime = FALSE,allow_emulation=TRUE), "compatible")
  expect_false(file.exists(file.path(path, "blupf90+")))
})

test_that("runner rejects unsafe controls and detects native process failures", {
  x <- make_example(1)
  expect_error(fit_blupf90_model(x$p, x$d, x$traits, executable = NA), "executable")
  expect_error(fit_blupf90_model(x$p, x$d, x$traits, executable = file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript"), timeout = 0), "timeout")
  expect_error(fit_blupf90_model(x$p, x$d, x$traits, executable = file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript"), keep_files = NA), "keep_files")
  expect_error(fit_blupf90_model(x$p, x$d, x$traits, executable = file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript"), job_root = "/missing"), "job_root")
  path <- tempfile()
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  local_mocked_bindings(run=function(...) list(status=1L),.package="processx")
  expect_error(run_native_blup("trusted", "bad.par", path, "failure", 5, 1), "failed")
  local_mocked_bindings(run=function(...) stop("cannot execute"),.package="processx")
  expect_error(run_native_blup("/missing/executable", "bad.par", path, "failure", 5, 1), "could not complete")
  expect_error(run_genetic_analysis(x$p, x$d, x$traits, backend = "blupf90", model_control = list(tol = 1)), "blupf90_control")
  expect_error(run_genetic_analysis(x$p, x$d, x$traits, backend = "sommer"), "BLUPF90")
  expect_error(run_genetic_analysis(x$p, x$d, x$traits, backend = "blupf90", impute = NA), "logical")
  expect_error(run_genetic_analysis(x$p, x$d, x$traits, backend = "blupf90", max_animals = 1), "max_animals")
})
