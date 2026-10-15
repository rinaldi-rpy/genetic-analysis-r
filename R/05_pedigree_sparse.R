# O(n log n) priority-queue topological ordering; no repeated pedigree scans.
pedigree_topological_order <- function(sire, dam, ranks) {
  n <- length(sire)
  degree <- as.integer(!is.na(sire)) + as.integer(!is.na(dam))
  children <- split(rep(seq_len(n), 2L), c(sire, dam))
  heap <- integer(n)
  size <- 0L
  push <- function(v) {
    size <<- size + 1L
    pos <- size
    while (pos > 1L) {
      parent <- pos %/% 2L
      if (ranks[heap[parent]] < ranks[v]) break
      heap[pos] <<- heap[parent]
      pos <- parent
    }
    heap[pos] <<- v
  }
  pop <- function() {
    top <- heap[1L]
    last <- heap[size]
    size <<- size - 1L
    pos <- 1L
    while (pos * 2L <= size) {
      child <- pos * 2L
      if (child < size && ranks[heap[child + 1L]] < ranks[heap[child]]) child <- child + 1L
      if (ranks[last] <= ranks[heap[child]]) break
      heap[pos] <<- heap[child]
      pos <- child
    }
    heap[pos] <<- last
    top
  }
  for (v in which(degree == 0L)) push(v)
  out <- integer(n)
  count <- 0L
  while (size > 0L) {
    v <- pop()
    count <- count + 1L
    out[count] <- v
    for (child in children[[as.character(v)]]) {
      degree[child] <- degree[child] - 1L
      if (degree[child] == 0L) push(child)
    }
  }
  if (count != n) stop("Pedigree contains a cycle")
  out
}

#' Compute exact pedigree inbreeding without a dense relationship matrix
#' @param pedigree Pedigree data frame.
#' @param max_pairs Maximum memoized ancestral relationships (memory guard).
#' @return Data frame with sorted ID, F, and F_percent.
#' @details Uses iterative tabular recursion with a bounded pair cache. Memory
#'   depends on pedigree connectivity; worst-case computation remains quadratic.
#' @export
pedigree_inbreeding <- function(pedigree, max_pairs = 2000000L) {
  p <- sort_pedigree(pedigree)
  sparse_inbreeding_sorted(p, max_pairs)
}

sparse_inbreeding_sorted <- function(p, max_pairs) {
  if (length(max_pairs) != 1L || !is.numeric(max_pairs) || !is.finite(max_pairs) || max_pairs < 1) {
    stop("Invalid max_pairs")
  }
  sire <- match(p$Sire, p$ID, nomatch = 0L)
  dam <- match(p$Dam, p$ID, nomatch = 0L)
  F <- numeric(nrow(p))
  cache <- new.env(hash = TRUE, parent = emptyenv())
  count <- 0L
  pair <- function(i, j) c(max(i, j), min(i, j))
  key <- function(v) paste(v, collapse = ":")
  known <- function(v) v[2L] == 0L || v[1L] == v[2L] || exists(key(v), cache, inherits = FALSE)
  value <- function(v) {
    if (v[2L] == 0L) {
      return(0)
    }
    if (v[1L] == v[2L]) {
      return(1 + F[v[1L]])
    }
    get(key(v), cache, inherits = FALSE)
  }
  relationship <- function(i, j) {
    target <- pair(i, j)
    stack <- list(target)
    while (length(stack)) {
      v <- stack[[length(stack)]]
      if (known(v)) {
        stack[[length(stack)]] <- NULL
        next
      }
      left <- pair(sire[v[1L]], v[2L])
      right <- pair(dam[v[1L]], v[2L])
      if (!known(left)) {
        stack[[length(stack) + 1L]] <- left
        next
      }
      if (!known(right)) {
        stack[[length(stack) + 1L]] <- right
        next
      }
      count <<- count + 1L
      if (count > max_pairs) stop("Inbreeding exceeds max_pairs cache guard")
      assign(key(v), (value(left) + value(right)) / 2, cache)
      stack[[length(stack)]] <- NULL
    }
    value(target)
  }
  for (i in seq_len(nrow(p))) if (sire[i] > 0L && dam[i] > 0L) F[i] <- relationship(sire[i], dam[i]) / 2
  data.frame(ID = p$ID, F = F, F_percent = 100 * F)
}

#' Construct a sparse exact Henderson pedigree precision matrix
#' @param pedigree Pedigree data frame.
#' @param max_pairs Maximum cached pairs for exact inbreeding.
#' @return List with sorted pedigree, sparse A_inverse and inbreeding.
#' @export
pedigree_precision <- function(pedigree, max_pairs = 2000000L) {
  p <- sort_pedigree(pedigree)
  n <- nrow(p)
  inb <- sparse_inbreeding_sorted(p, max_pairs)
  sire <- match(p$Sire, p$ID, nomatch = 0L)
  dam <- match(p$Dam, p$ID, nomatch = 0L)
  F <- c(0, inb$F)
  D <- 1 - (!is.na(p$Sire)) * (1 + F[sire + 1L]) / 4 - (!is.na(p$Dam)) * (1 + F[dam + 1L]) / 4
  if (any(D <= 0)) stop("Non-positive Mendelian sampling variance")
  ii <- jj <- integer(9L * n)
  xx <- numeric(9L * n)
  end <- 0L
  for (a in seq_len(n)) {
    parents <- c(sire[a], dam[a])
    parents <- parents[parents > 0L]
    ids <- c(a, parents)
    w <- c(1, rep(-0.5, length(parents)))
    len <- length(ids)^2L
    at <- end + seq_len(len)
    ii[at] <- rep(ids, times = length(ids))
    jj[at] <- rep(ids, each = length(ids))
    xx[at] <- as.vector(outer(w, w) / D[a])
    end <- end + len
  }
  at <- seq_len(end)
  Ai <- Matrix::sparseMatrix(i = ii[at], j = jj[at], x = xx[at], dims = c(n, n), dimnames = list(p$ID, p$ID))
  list(pedigree = p, A_inverse = Ai, inbreeding = inb)
}
