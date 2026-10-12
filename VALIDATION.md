# Release 1.1.0 validation

Validated on Linux x86-64 on 12 October 2026, with R 4.6.1, sommer 4.4.87,
and the official BLUPF90+ Linux binary downloaded during this task.

- R CMD build: passed.
- R CMD check --no-manual: Status OK, no errors, warnings or notes.
- Testthat: 256 assertions passed, zero failures, warnings or skips.
- covr package line coverage: 87.87% (real native tests enabled).
- Real BLUPF90+ REML/BLUP pipeline: 1, 2, 3 and 5 traits converged and returned
  all pedigree EBVs/SEs, trait heritabilities/SEs, and every pairwise rg/SE.
- Preparation and output-mapping tests also cover 6 and 10 traits.
- Native known-variance EBVs and SEs match exact R marginal-model equations
  on an inbred pedigree to tolerance 1e-6.
- Sparse A-inverse agrees with the dense inverse on randomized inbred pedigrees
  to tolerance 1e-10, with exact F agreement.
- Existing real sommer single-, two-, and four-trait tests pass.
- Native parser rejects missing SEs, nonconvergence and incomplete ID coverage.
- Explicit installer end-to-end: verified HTTPS binary download, ELF validation,
  private Intel runtime installation, executable launcher and --help check passed.
- This-device setup and a three-trait example with sequential EM ran successfully.

## Validation environment

BLUPF90+ and private Intel MKL/OpenMP 2025.3.0 runtime libraries were installed
separately from this repository for real native tests. The official executable
passed its --help check and the Linux installer completed successfully. Users
install the R package and external native backend separately; no local device
paths or installed libraries are included in the repository or source archive.

Native binary SHA256 used during validation:
aa4039297a3d6b578b9ac5f711d9b6f1ed8c0fedc04fa470f01421c4779dd499

## Limits

BLUPF90+ uses native Monte Carlo approximations for genetic parameter SEs; sommer
uses first-order delta SEs. Both report prediction-error SDs as EBV SEs. Neither
propagates single-imputation uncertainty. No performance speedup is guaranteed;
no large national-dataset benchmark or remote GitHub Actions run was performed.

Trait count is configurable, with a default native resource guard of 32 traits.
Unstructured covariance estimation and exact inbreeding can still become costly
or unidentifiable. Failed/inadmissible fits raise errors. Linux installer and
current native output layout were tested; other operating systems are not
validated by this delivery. No GitHub repository was pushed or web service deployed.

## Automatic founder completion regression

An offspring-only ID/Sire/Dam table now expands to include absent named parents
as unknown-ancestry founders and sorts them before offspring. Tests verify
unique parent insertion, typed missing metadata, optional dates, idempotence,
full-sib mating F, dense/sparse inverse agreement, and expanded-pedigree guards.
A real BW6 native analysis with fixed=~Age_6Month+Sex passed on simulated data
with 120 offspring-only rows, automatically adding 22 parents and returning
142 pedigree predictions. This regression is included in release 1.1.0.

## EBV accuracy regression

EBV tables for both backends now include PEV, Reliability, Accuracy and
Accuracy_percent. The helper uses SE squared and the inbreeding-adjusted prior
variance (1+F)*G_tt. Analytic tests verify 80% accuracy when G_tt=25, SE=3 and
F=0, and the correct adjustment for F=0.25. Tests cover trait/ID order mapping,
zero information, perfect prediction, zero genetic variance, rounded reliability
and inconsistent PEV. Both real backends return these fields. Existing model
results can be enriched with add_ebv_accuracy() without a new fit.

## Release metadata

Version 1.1.0 is the requested public release numbering for this consolidated
implementation. Maintainer name/email and GitHub account remain editable
placeholders at the user's request. No GitHub upload or remote CI run was
performed by this task.
