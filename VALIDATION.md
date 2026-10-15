# Release 1.2.0 validation

Validated on Linux x86-64 on 15 October 2026 with R 4.6.1 and the separately
installed official BLUPF90+ executable. This release has no sommer dependency
or fitting implementation.

- R CMD build: passed.
- R CMD check --no-manual: Status OK.
- Testthat with the native executable enabled: 328 assertions passed; no failures,
  warnings or skips.
- covr line coverage: 88.39% (native integration tests enabled).
- Real native direct animal REML/BLUP runs cover 1, 2, 3 and 5 traits; preparation
  and parser tests also cover 6 and 10 traits.
- Known-variance direct, maternal genetic and independent-group BLUPs and SEs
  match exact mixed-model marginal equations to tolerance 1e-6, including a
  two-trait maternal model with nonzero cross-trait covariances.
- Real maternal genetic plus maternal permanent-environment REML completes and
  returns all effect covariances, solutions, direct h2 and extra variance ratios
  with native approximate SEs. The direct h2 denominator includes all components.
- Selected EM predictor tests verify immutable observations, omission of unselected
  auxiliary traits, per-target mappings, missing-predictor fallback, and diagnostics.
- Existing bidirectional Gaussian EM analytic checks, pedigree founder completion,
  sparse/dense A-inverse agreement, exclusion audits and EBV accuracy tests pass.
- Installer tests mock downloads and validate Linux ELF, Windows PE x86-64 and
  Intel macOS Mach-O headers, platform selection and architecture rejection.
- The previously installed Linux binary and private Intel runtimes are used for
  native tests. No new Windows/macOS binary download or native execution was
  validated on this device; those platforms need suitable external executables
  and runtimes. Apple Silicon Intel execution requires Rosetta.

SHA256 of the Linux native binary used:
`aa4039297a3d6b578b9ac5f711d9b6f1ed8c0fedc04fa470f01421c4779dd499`

CI is configured for R package checks on Linux, Windows and macOS, plus a Linux
native-execution/coverage job. Those remote jobs have not been run in this chat.
No GitHub upload or large-dataset performance benchmark was performed.

Direct-maternal genetic covariance can be estimated explicitly; other added
random terms remain independent. A single-record maternal simulation demonstrates execution and
output consistency, not accurate component recovery for every family structure.
Single imputation and variance-component estimation uncertainty are not propagated
into prediction SEs or accuracy. Genetic parameter SEs use BLUPF90+'s native
Monte Carlo approximation. Failed convergence or inadmissible output is an error.

The source archive and repository exclude installed R libraries, native binaries,
Intel runtimes, private datasets and retained analysis jobs. Maintainer identity
remains a placeholder for the repository owner.


## Correlated direct-maternal extension (same release)

- Known-variance native BLUPs/SEs agree with exact joint mixed-model equations
  within 1e-6 for one and two traits, positive and negative cross-effect covariance,
  and a nonsymmetric two-trait cross-effect block.
- Real single- and two-trait REML runs on a simulated six-generation pedigree
  converge and return nonzero cross-effect covariances, direct-maternal correlations
  and native approximate SEs. All full-pedigree direct/maternal predictions are retained.
- Tests verify covariance-group ordering, parameter labels, effect reordering,
  invalid/non-PSD blocks, and covariance-adjusted heritability denominators.
- The native joint matrix is ordered effect first, trait within effect. The parser
  supports BLUPF90+'s header naming the first effect of a correlated group.
- Reported heritabilities use a reference non-inbred offspring/dam relationship
  of 0.5, adding C_dm once to phenotypic variance; the BLUP model uses full A.
