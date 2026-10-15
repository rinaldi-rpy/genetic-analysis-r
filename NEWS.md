# geneticAnalysisR 1.2.0

- BLUPF90+ is the sole fitting engine; remove sommer code and dependency.
- Add shared/per-target numerical EM predictor selection; NULL preserves joint EM.
- Add maternal genetic, maternal environment and other group terms,
  full effect solutions/covariances and component variance ratios with native SEs.
- Estimate full direct-maternal covariance with pedigree_effect("Dam",
  correlate_with="animal"), including cross-trait terms and native covariance/
  correlation SEs; propagate the joint covariance into final BLUP predictions.
- Adjust direct/maternal heritability denominators for direct-maternal covariance
  and all additional random components under the reference non-inbred convention.
- Add Windows executable discovery, platform-specific download selection, binary
  architecture validation and portable SHA256; explicit Rosetta opt-in on Apple Silicon.
- Configure Linux, Windows and macOS R checks and separate Linux native tests.
- Keep historical version 1.1.2 artifacts separately for reproducibility.
