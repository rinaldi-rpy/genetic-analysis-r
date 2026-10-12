# Contributing

Keep pedigree processing, imputation, model fitting and orchestration separate.
Add regression tests before changing statistical behavior. Preserve observed
phenotypes and deterministic ordering. Do not coerce invalid covariance matrices
or report nonconverged fits as successful analyses.

Use R >= 4.1. Sommer is optional; its tested version is 4.4.87.
Set BLUPF90_PLUS to a verified BLUPF90+ executable to run native integration tests.
Native tests must run before modifying parameter-file or log-parser behavior.
Keep downloaded executables and runtime libraries outside this repository. Run:

```r
roxygen2::roxygenise()
testthat::test_local()
coverage <- covr::package_coverage()
stopifnot(covr::percent_coverage(coverage) >= 80)
```

Build and check the source package using R CMD build and R CMD check --no-manual.
Changes to solver compatibility require real end-to-end single-trait and
multi-trait tests, including residual covariance and unphenotyped pedigree IDs.
Document any change to assumptions, uncertainty definitions or output schemas.
Use synthetic data for issues and tests; do not upload private animal records.
