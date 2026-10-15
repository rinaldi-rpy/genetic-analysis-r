#' Sort a pedigree with parents preceding progeny
#' @param pedigree Data frame with ID, Sire, Dam, and optional Birth_Date.
#'   Unknown parents are NA, empty strings, or "0". IDs are character strings.
#' @return A sorted data frame, including automatically added missing parents.
#'   Missing dates sort last; ID breaks ties.
#' @details Named parents absent from ID are added as founders with unknown
#'   ancestry and missing metadata. Known relationships are preserved.
#' @export
sort_pedigree <- function(pedigree) {
  if (!is.data.frame(pedigree) || !nrow(pedigree) || anyDuplicated(names(pedigree)) ||
    !all(c("ID", "Sire", "Dam") %in% names(pedigree))) {
    stop("pedigree requires ID, Sire and Dam and at least one row")
  }
  p <- pedigree
  for (v in c("ID", "Sire", "Dam")) p[[v]] <- trimws(as.character(p[[v]]))
  if (anyNA(p$ID) || any(p$ID %in% c("", "0")) || anyDuplicated(p$ID)) {
    stop("ID must be nonmissing, nonzero and unique")
  }
  for (v in c("Sire", "Dam")) p[[v]][p[[v]] %in% c("", "0")] <- NA_character_
  if (any(p$ID == p$Sire, na.rm = TRUE) || any(p$ID == p$Dam, na.rm = TRUE)) {
    stop("An individual cannot be its own parent")
  }
  if (any(p$Sire == p$Dam, na.rm = TRUE)) stop("Sire and Dam must differ")
  parents <- unique(na.omit(c(p$Sire, p$Dam)))
  missing_parents <- setdiff(parents, p$ID)
  if (!"Birth_Date" %in% names(p)) p$Birth_Date <- as.Date(NA)
  if (!inherits(p$Birth_Date, "Date")) {
    raw <- p$Birth_Date
    if (!is.character(raw) && !is.factor(raw)) stop("Birth_Date must be Date or ISO date text")
    dates <- suppressWarnings(as.Date(as.character(raw), format = "%Y-%m-%d"))
    if (any(!is.na(raw) & (is.na(dates) | as.character(dates) != as.character(raw)))) stop("Invalid Birth_Date; use YYYY-MM-DD")
    p$Birth_Date <- dates
  }
  if (length(missing_parents)) {
    # Typed NA rows preserve Date, factor and other metadata column classes.
    added <- p[rep(NA_integer_, length(missing_parents)), , drop = FALSE]
    added$ID <- missing_parents
    p <- rbind(p, added)
  }
  for (v in c("Sire", "Dam")) {
    pd <- p$Birth_Date[match(p[[v]], p$ID)]
    if (any(pd >= p$Birth_Date, na.rm = TRUE)) stop("Parent birth date must precede offspring")
  }
  sire <- match(p$Sire, p$ID)
  dam <- match(p$Dam, p$ID)
  founders <- is.na(sire) & is.na(dam)
  key <- as.numeric(p$Birth_Date)
  parent_index <- c(sire, dam)
  child_dates <- rep(key, 2L)
  valid <- !is.na(parent_index) & !is.na(child_dates)
  if (any(valid)) {
    oldest <- tapply(child_dates[valid], parent_index[valid], min)
    missing_founders <- which(founders & is.na(key))
    key[missing_founders] <- oldest[as.character(missing_founders)]
  }
  ranks <- integer(nrow(p))
  ranks[order(!founders, key, p$ID, na.last = TRUE, method = "radix")] <- seq_len(nrow(p))
  idx <- pedigree_topological_order(sire, dam, ranks)
  p <- p[idx, , drop = FALSE]
  rownames(p) <- NULL
  p
}

#' Construct additive relationships and inbreeding
#' @param pedigree Pedigree accepted by sort_pedigree.
#' @param inverse Whether to also compute the dense inverse.
#' @param max_animals Dense-memory guard; increase deliberately.
#' @return List with pedigree, A, optional A_inverse, and inbreeding table.
#' @export
pedigree_relationship <- function(pedigree, inverse = TRUE, max_animals = 2000L) {
  if (!is.data.frame(pedigree)) stop("pedigree requires a data frame")
  n <- nrow(pedigree)
  if (!is.numeric(max_animals) || length(max_animals) != 1L ||
    !is.finite(max_animals) || max_animals < 1 || n > max_animals) {
    stop("Pedigree exceeds max_animals dense-memory guard")
  }
  p <- sort_pedigree(pedigree)
  n <- nrow(p)
  if (n > max_animals) stop("Expanded pedigree exceeds max_animals dense-memory guard")
  A <- matrix(0, n, n, dimnames = list(p$ID, p$ID))
  sire <- match(p$Sire, p$ID)
  dam <- match(p$Dam, p$ID)
  for (i in seq_len(n)) {
    known <- c(sire[i], dam[i])
    known <- known[!is.na(known)]
    if (i > 1L && length(known)) {
      j <- seq_len(i - 1L)
      A[i, j] <- colSums(A[known, j, drop = FALSE]) / 2
      A[j, i] <- A[i, j]
    }
    A[i, i] <- 1 + if (length(known) == 2L) A[sire[i], dam[i]] / 2 else 0
  }
  Ai <- if (inverse) chol2inv(chol(A)) else NULL
  if (inverse) dimnames(Ai) <- dimnames(A)
  list(
    pedigree = p, A = A, A_inverse = Ai,
    inbreeding = data.frame(ID = p$ID, F = diag(A) - 1, F_percent = 100 * (diag(A) - 1))
  )
}
