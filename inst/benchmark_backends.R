# Run example_script.R first to define pedigree, phenotypes and traits.
library(geneticAnalysisR)
completed <- sequential_em(phenotypes,traits)$data
native_time <- system.time(native <- run_genetic_analysis(
  pedigree,completed,traits,fixed=~Sex,impute=FALSE,backend="blupf90"))
sommer_time <- system.time(alternative <- run_genetic_analysis(
  pedigree,completed,traits,fixed=~Sex,impute=FALSE,backend="sommer"))
print(rbind(BLUPF90=native_time,sommer=sommer_time))
print(native$heritability)
print(alternative$heritability)
# Speed depends on trait count, pedigree, solver settings and thread counts.
# SE methods differ, as recorded in diagnostics.
