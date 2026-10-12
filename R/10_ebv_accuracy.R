#' Add prediction-error based EBV reliability and accuracy
#' @param ebv Data frame with ID, trait, EBV and SE. SE must be the
#'   prediction-error standard deviation, as returned by this package.
#' @param genetic_covariance Named genetic covariance matrix; its diagonal
#'   contains additive genetic variances on each trait's input scale.
#' @param inbreeding Data frame with unique ID and F (fraction, not percentage).
#' @param tolerance Dimensionless tolerance for small negative reliability
#'   caused by native printed covariance/SE rounding.
#' @return Input table with PEV, Reliability, Accuracy and Accuracy_percent.
#'   Existing columns and row order are preserved. Zero genetic variance gives NA.
#' @details PEV = SE^2; reliability = 1 - PEV / ((1 + F) * genetic variance);
#'   accuracy is sqrt(reliability). Only minor rounding excursions are clipped.
#'   Materially negative reliability is undefined (NA), with a warning. These
#'   model-based values condition on fitted covariance components and do not
#'   account for single-imputation or variance-component estimation uncertainty.
#' @export
add_ebv_accuracy <- function(ebv, genetic_covariance, inbreeding, tolerance = 1e-4) {
  if (!is.data.frame(ebv) || !all(c("ID", "trait", "EBV", "SE") %in% names(ebv)) ||
    !is.numeric(ebv$SE) || any(!is.finite(ebv$SE) & !is.na(ebv$SE)) ||
    any(ebv$SE < 0, na.rm = TRUE)) {
    stop("Provide an EBV table with valid prediction-error SEs")
  }
  G <- genetic_covariance
  if (!is.matrix(G) || !is.numeric(G) || nrow(G) != ncol(G) ||
    is.null(colnames(G)) || anyDuplicated(colnames(G)) ||
    !identical(rownames(G), colnames(G)) || !isSymmetric(G) ||
    any(!is.finite(G)) || any(diag(G) < 0) || anyNA(ebv$trait) ||
    !all(as.character(ebv$trait) %in% colnames(G))) {
    stop("Provide a named genetic covariance matrix covering every EBV trait")
  }
  if (!is.data.frame(inbreeding) || !all(c("ID", "F") %in% names(inbreeding)) ||
    anyNA(inbreeding$ID) || anyDuplicated(inbreeding$ID) ||
    !is.numeric(inbreeding$F) || any(!is.finite(inbreeding$F)) ||
    any(inbreeding$F < 0 | inbreeding$F > 1) || anyNA(ebv$ID) ||
    !all(as.character(ebv$ID) %in% as.character(inbreeding$ID))) {
    stop("Provide inbreeding fractions for every EBV ID")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) ||
    tolerance < 0 || tolerance >= 1) {
    stop("Invalid accuracy rounding tolerance")
  }
  F <- inbreeding$F[match(as.character(ebv$ID), as.character(inbreeding$ID))]
  variance <- diag(G)[match(as.character(ebv$trait), colnames(G))]
  prior <- (1 + F) * variance
  pev <- ebv$SE^2
  reliability <- rep(NA_real_, nrow(ebv))
  valid <- prior > 0 & is.finite(prior) & is.finite(pev)
  reliability[valid] <- 1 - pev[valid] / prior[valid]
  invalid <- !is.na(reliability) & reliability < -tolerance
  if (any(invalid)) warning("EBV PEV exceeds prior genetic variance; accuracy is NA for affected rows", call. = FALSE)
  reliability[invalid] <- NA_real_
  reliability <- pmin(1, pmax(0, reliability))
  out <- ebv
  out$PEV <- pev
  out$Reliability <- reliability
  out$Accuracy <- sqrt(reliability)
  out$Accuracy_percent <- 100 * out$Accuracy
  out
}
