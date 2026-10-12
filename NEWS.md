# geneticAnalysisR 1.1.2

* Default to bidirectional joint Gaussian EM for partially observed traits.
* Perform two E/M updates per forward/reverse cycle, including conditional covariance.
* Require at least ten cycles before correlation-difference convergence, with a parameter guard.
* Return correlation histories, distribution moments and phenotype correlation matrices.
* Preserve all-missing-row exclusion and optional legacy sequential chaining.

# geneticAnalysisR 1.1.1

* Exclude phenotype rows missing every selected trait before EM/fitting, while
  retaining their animals in the full pedigree.
* Order EM by initial observed trait counts, with deterministic supplied-order ties.
* Add trait_order="supplied" for explicit biological/specified chaining.
* Return excluded IDs/rows/indices, retained indices, observation counts and chain diagnostics.
* Preserve model trait order and all observed measurements.

# geneticAnalysisR 1.1.0

* Select BLUPF90+ or sommer through one high-level analysis API.
* Support configurable trait counts and ordered sequential Gaussian EM.
* Accept ID/Sire/Dam tables with automatically completed founder rows.
* Construct exact sparse pedigree precision and inbreeding for native analysis.
* Return EBVs, SEs, PEVs, reliabilities, accuracies and accuracy percentages.
* Export add_ebv_accuracy() for existing results without refitting.
* Return trait heritabilities and pairwise genetic correlations with SEs.
* Add an explicit Linux installer with private Intel runtimes and verified HTTPS.
* Isolate native jobs with timeouts, thread controls and auditable ID mappings.
* Document CSV inputs, data dictionaries, database schemas and GitHub publishing.
* Include generated R documentation, tests, validation and CI configuration.
