# Use your own prepared data frames. Keep the full pedigree including ancestors.
# Maternal models require known dams for every phenotype row used for fitting.
library(geneticAnalysisR)

# phenotypes must include the outcomes and all fixed predictors.
# Example with explicit EM predictor selection and correlated maternal genetics and independent maternal environment:
# result <- run_genetic_analysis(
#   pedigree, phenotypes,
#   traits=c("BW6", "BL6", "CG6"), fixed=~Sex+Age_6Month,
#   predictor_traits=list(BW6="BL6", BL6="BW6", CG6=c("BW6", "BL6")),
#   random=list(maternal=pedigree_effect("Dam", correlate_with="animal"),
#               maternal_environment=iid_effect("Dam")),
#   blupf90_control=list(keep_files=TRUE)
# )
# result$random_solutions$maternal
# result$variance_ratios
# result$imputation$target_models$CG6$diagnostics

# result$direct_maternal_parameters
# result$direct_maternal_covariance
