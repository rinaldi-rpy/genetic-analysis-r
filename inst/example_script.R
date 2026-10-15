library(geneticAnalysisR)
set.seed(103)
k <- 3L # Change to 1, 2, 5, or another trait count.
pedigree <- data.frame(
  ID=c(paste0("S",1:10),paste0("D",1:12),paste0("O",1:120)),
  Sire=c(rep(NA,22),rep(paste0("S",1:10),12)),
  Dam=c(rep(NA,22),rep(paste0("D",1:12),each=10)))
# A is used here only to SIMULATE correlated breeding values, not by native fitting.
rel <- pedigree_relationship(pedigree)
G <- diag(k)*2+matrix(.5,k,k)
u <- t(chol(G)) %*% matrix(rnorm(nrow(pedigree)*k),k) %*% chol(rel$A)
idx <- match(paste0("O",1:120),rel$pedigree$ID)
phenotypes <- data.frame(ID=paste0("O",1:120),Sex=rep(c("F","M"),60))
traits <- paste0("Trait",seq_len(k))
for (j in seq_len(k)) phenotypes[[traits[j]]] <- 20*j+u[j,idx]+rnorm(120)
if (k>1L) phenotypes[[traits[2L]]][c(3,9,21)] <- NA
# Install once if necessary: install_blupf90()
result <- run_genetic_analysis(pedigree,phenotypes,traits,fixed=~Sex,backend="blupf90")
print(result)
head(result$ebv)
result$genetic_correlations
result$imputation$diagnostics
