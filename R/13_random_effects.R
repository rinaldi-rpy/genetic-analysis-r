#' Define an additional pedigree genetic effect
#' @param column Phenotype column holding pedigree IDs, or Sire/Dam looked up from pedigree.
#' @param correlate_with NULL for independence, or 'animal' to estimate direct-maternal covariance.
#' @return A validated random-effect specification for the random argument.
#' @details Use correlate_with='animal' on a Dam effect to estimate its full
#'   cross-trait covariance with direct genetics. Other effects stay independent.
#' @export
pedigree_effect <- function(column = "Dam", correlate_with = NULL) {
  if (!is.null(correlate_with) && !identical(correlate_with,"animal")) stop("correlate_with must be NULL or animal")
  if (!is.null(correlate_with) && !identical(column,"Dam")) stop("Direct-maternal correlation requires the pedigree Dam effect")
  spec <- random_effect_spec(column,"pedigree")
  spec$correlate_with <- correlate_with
  spec
}

#' Define an independent group or permanent-environment effect
#' @param column Group ID column, or Sire/Dam looked up from pedigree.
#' @return A random-effect specification for the random argument.
#' @export
iid_effect <- function(column) random_effect_spec(column, "iid")

random_effect_spec <- function(column, type) {
  if (!is.character(column) || length(column)!=1L || is.na(column) || !nzchar(column)) {
    stop("Provide one effect column name")
  }
  structure(list(column=column,type=type), class="blupf90_effect")
}

build_random_effects <- function(random, data, rel, X, S) {
  if (!is.list(random) || (length(random) && (is.null(names(random)) ||
      anyDuplicated(names(random)) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$",names(random))) ||
      "animal" %in% names(random)))) stop("random must be a uniquely named effect list; animal is reserved")
  correlated <- names(random)[vapply(random,function(r) inherits(r,"blupf90_effect") && identical(r$correlate_with,"animal"),logical(1))]
  if(length(correlated)>1L) stop("Only one maternal term can correlate with animal")
  random <- random[c(correlated,setdiff(names(random),correlated))]
  specs <- c(list(animal=pedigree_effect("ID")), random)
  effects <- list()
  k <- ncol(S)
  for (name in names(specs)) {
    spec <- specs[[name]]
    if (!inherits(spec,"blupf90_effect") || any(!names(spec) %in% c("column","type","correlate_with")) ||
        length(spec$type)!=1L || !spec$type %in% c("pedigree","iid")) stop("Invalid random effect specification")
    random_effect_spec(spec$column,spec$type)
    if(!is.null(spec$correlate_with) && (spec$type!="pedigree" || spec$column!="Dam" || !identical(spec$correlate_with,"animal"))) stop("Invalid direct-maternal effect specification")
    column <- spec$column
    if (column %in% c("Sire","Dam")) {
      values <- rel$pedigree[[column]][match(as.character(data$ID),rel$pedigree$ID)]
    } else {
      if (!column %in% names(data)) stop("Missing random effect column: ",column)
      values <- as.character(data[[column]])
    }
    values <- trimws(values)
    # Unknown dams are never pooled into a fictitious common mother.
    if (anyNA(values) || any(values %in% c("","0"))) stop("Missing IDs for random effect ",name,"; exclude these phenotype rows or provide known groups")
    levels <- if (spec$type=="pedigree") rel$pedigree$ID else sort(unique(values))
    if (!all(values %in% levels)) stop("Random pedigree IDs must occur in the pedigree: ",name)
    if (spec$type=="iid" && (length(levels)<2L || length(levels)>=nrow(data))) {
      stop("IID effects need at least two groups and replicated levels: ",name)
    }
    effect <- ncol(X)+length(effects)+1L
    effects[[name]] <- list(name=name,type=spec$type,column=column,effect=effect,correlate_with=spec$correlate_with,
      codes=match(values,levels),map=data.frame(code=seq_along(levels),ID=levels),
      start=S*(if(name=="animal") .3 else .1/length(random)))
  }
  effects
}
