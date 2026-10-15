#' Reshape wide phenotypes to trait-major long form
#' @param data Wide phenotype data, with one row per ID.
#' @param traits Ordered numeric trait names.
#' @return Long data with .ga_trait, .ga_y, .ga_record, and .ga_animal.
#' @export
reshape_traits <- function(data, traits) {
  validate_traits(data, traits)
  if (!"ID" %in% names(data) || anyNA(data$ID) || anyDuplicated(as.character(data$ID))) {
    stop("Phenotype ID must be present, unique and nonmissing")
  }
  if (any(startsWith(names(data), ".ga_"))) stop("Reserved .ga_ column names in data")
  n <- nrow(data)
  k <- length(traits)
  long <- data[rep(seq_len(n), k), setdiff(names(data), traits), drop = FALSE]
  long$.ga_trait <- factor(rep(paste0("T", seq_len(k)), each = n), levels = paste0("T", seq_len(k)))
  long$.ga_y <- unlist(data[traits], use.names = FALSE)
  long$.ga_record <- factor(rep(paste0("R", seq_len(n)), k), levels = paste0("R", seq_len(n)))
  long$.ga_animal <- factor(as.character(long$ID))
  rownames(long) <- NULL
  long
}

