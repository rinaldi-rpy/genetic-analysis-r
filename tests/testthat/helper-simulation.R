make_example <- function(k = 2L) {
  set.seed(if (k == 2L) 103 else 102)
  p <- data.frame(
    ID = c(paste0("S", 1:10), paste0("D", 1:12), paste0("O", 1:120)),
    Sire = c(rep(NA, 22), rep(paste0("S", 1:10), 12)),
    Dam = c(rep(NA, 22), rep(paste0("D", 1:12), each = 10))
  )
  rel <- pedigree_relationship(p)
  G <- diag(k) * 2 + matrix(0.5, k, k)
  u <- t(chol(G)) %*% matrix(rnorm(nrow(p) * k), k) %*% chol(rel$A)
  ix <- match(paste0("O", 1:120), rel$pedigree$ID)
  d <- data.frame(ID = paste0("O", 1:120), Sex = rep(c("F", "M"), 60))
  for (j in seq_len(k)) d[[paste0("w", j)]] <- 20 * j + u[j, ix] + rnorm(120)
  list(p = p, d = d, rel = rel, traits = paste0("w", seq_len(k)))
}
