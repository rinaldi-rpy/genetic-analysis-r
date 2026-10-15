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

full_pedigree_blup <- function(long, X, A, G, E, traits) {
  k <- length(traits)
  m <- nrow(A)
  n <- nrow(long) / k
  C <- kronecker(G, A)
  selected <- unlist(lapply(seq_len(k), function(j) {
    match(as.character(long$ID[seq_len(n)]), rownames(A)) + (j - 1L) * m
  }))
  B <- C[, selected, drop = FALSE]
  V <- C[selected, selected, drop = FALSE] + kronecker(E, diag(n))
  Vi <- chol2inv(chol(V))
  VX <- Vi %*% X
  P <- Vi - VX %*% solve(crossprod(X, VX), t(VX))
  u <- drop(B %*% P %*% long$.ga_y)
  pev <- diag(C) - rowSums((B %*% P) * B)
  if (any(pev < -1e-6)) stop("Negative prediction error variance")
  data.frame(ID = rep(rownames(A), k), trait = rep(traits, each = m), EBV = u, SE = sqrt(pmax(0, pev)))
}

