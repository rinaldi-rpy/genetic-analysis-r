# Direct and maternal groups are ordered effect-major, with traits within effect.
configure_direct_maternal <- function(job) {
  names <- names(job$random_effects)[vapply(job$random_effects,
    function(r) identical(r$correlate_with,"animal"),logical(1))]
  job$direct_maternal_effect <- if(length(names)) names[[1L]] else NULL
  if(!length(names)) return(job)
  m <- job$random_effects[[names[[1L]]]]
  # A small nonzero PD initializer allows every cross-effect covariance to be estimated.
  job$direct_maternal_start <- job$G_start * (.05*sqrt(m$start[1,1]/job$G_start[1,1]))
  pairs <- expand.grid(trait1=seq_len(job$n_traits),trait2=seq_len(job$n_traits))
  cross <- data.frame(component="direct_maternal",pairs,
    label=paste0("G_",job$animal_effect,"_",m$effect,"_",pairs$trait1,"_",pairs$trait2))
  job$parameter_map <- rbind(job$parameter_map,cross)
  job
}

dm_label <- function(job,i,j) {
  paste0("G_",job$animal_effect,"_",job$random_effects[[job$direct_maternal_effect]]$effect,"_",i,"_",j)
}

direct_maternal_block <- function(G,M,C) {
  k <- nrow(G)
  if(!is.matrix(M) || !is.matrix(C) || !identical(dim(M),dim(G)) ||
      !identical(dim(C),dim(G)) || any(!is.finite(c(G,M,C)))) stop("Invalid direct-maternal covariance dimensions/values")
  block <- rbind(cbind(G,C),cbind(t(C),M))
  if(!isSymmetric(block,tol=1e-6) || min(eigen(block,symmetric=TRUE,only.values=TRUE)$values)< -1e-7) stop("Joint direct-maternal covariance is not admissible")
  block
}

direct_maternal_se_lines <- function(job) {
  if(is.null(job$direct_maternal_effect)) return(character())
  m <- job$random_effects[[job$direct_maternal_effect]]$effect
  a <- job$animal_effect
  lines <- character()
  for(i in seq_len(job$n_traits)) for(j in seq_len(job$n_traits)) {
    covariance <- dm_label(job,i,j)
    correlation <- paste0(covariance,"/sqrt(G_",a,"_",a,"_",i,"_",i,"*G_",m,"_",m,"_",j,"_",j,")")
    lines <- c(lines,paste("OPTION se_covar_function",paste0("covam_",i,"_",j),covariance),
      paste("OPTION se_covar_function",paste0("ram_",i,"_",j),correlation))
  }
  lines
}

read_random_covariances <- function(job,lines,headers) {
  covariance <- list()
  C <- joint <- NULL
  for(name in names(job$random_effects)) {
    r <- job$random_effects[[name]]
    if(!is.null(job$direct_maternal_effect) && name==job$direct_maternal_effect) next
    group <- r$effect
    correlated <- name=="animal" && !is.null(job$direct_maternal_effect)
    if(correlated) group <- c(group,job$random_effects[[job$direct_maternal_effect]]$effect)
    header <- headers[vapply(headers,function(h) {
      numbers<-regmatches(lines[h],gregexpr("[0-9]+",lines[h]))[[1]]
      identical(as.integer(numbers),as.integer(group)) ||
        (correlated && identical(as.integer(numbers),as.integer(group[1L])))
    },logical(1))]
    if(length(header)!=1L) stop("Unsupported BLUPF90 covariance headers for effect ",name)
    M <- read_native_matrix(lines,header,job$n_traits*length(group))
    if(!isSymmetric(M,tol=1e-6) || min(eigen(M,symmetric=TRUE,only.values=TRUE)$values)< -1e-7) stop("Random covariance matrix is not admissible")
    if(correlated) {
      k <- job$n_traits; a <- seq_len(k); b <- k+a
      covariance$animal <- M[a,a,drop=FALSE]
      covariance[[job$direct_maternal_effect]] <- M[b,b,drop=FALSE]
      C <- M[a,b,drop=FALSE]
      joint <- M
    } else covariance[[name]] <- M
  }
  covariance <- covariance[names(job$random_effects)]
  covariance <- lapply(covariance,function(M) {dimnames(M)<-list(job$traits,job$traits);M})
  if(!is.null(C)) dimnames(C)<-list(job$traits,job$traits)
  if(!is.null(joint)) {
    labels<-c(paste0("direct:",job$traits),paste0("maternal:",job$traits))
    dimnames(joint)<-list(labels,labels)
  }
  list(covariances=covariance,cross=C,joint=joint)
}

read_direct_maternal_parameters <- function(job,lines,G,M,C) {
  empty <- data.frame(direct_trait=character(),maternal_trait=character(),covariance=numeric(),covariance_SE=numeric(),correlation=numeric(),correlation_SE=numeric())
  if(is.null(C)) return(empty)
  results <- list()
  for(i in seq_len(job$n_traits)) for(j in seq_len(job$n_traits)) {
    value <- C[i,j]
    se <- read_native_function_se(lines,paste0("covam_",i,"_",j),value)
    if(G[i,i]<=0 || M[j,j]<=0) { cor <- cor_se <- NA_real_ } else {
      cor <- value/sqrt(G[i,i]*M[j,j])
      cor_se <- read_native_function_se(lines,paste0("ram_",i,"_",j),cor)
    }
    results[[length(results)+1L]] <- data.frame(direct_trait=job$traits[i],maternal_trait=job$traits[j],
      covariance=value,covariance_SE=se,correlation=cor,correlation_SE=cor_se)
  }
  do.call(rbind,results)
}
