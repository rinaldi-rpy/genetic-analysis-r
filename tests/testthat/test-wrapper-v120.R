test_that("selected EM predictors preserve observed cells and omit unselected predictors", {
  set.seed(71)
  d <- data.frame(ID=1:180, aux=rnorm(180), other=rnorm(180))
  d$x <- 2*d$aux+rnorm(180); d$y <- d$aux+rnorm(180)
  d$x[1:12] <- NA; d$y[13:24] <- NA
  r <- sequential_em(d,c("x","y"),predictor_traits="aux")
  expect_false(anyNA(r$data[c("x","y")]))
  expect_equal(r$data$x[13:180],d$x[13:180])
  expect_true(all(r$diagnostics$predictor_traits=="aux"))
  expect_equal(names(r$target_models),c("x","y"))
  d$other <- d$other * 100
  expect_equal(sequential_em(d,c("x","y"),predictor_traits="aux")$data,r$data |> transform(other=d$other))
  expect_error(sequential_em(d,c("x","y"),predictor_traits="missing"),"trait")
})
test_that("maternal and group effects generate separate native random groups", {
  p <- data.frame(ID=letters[1:8],Sire=c(NA,NA,rep("a",6)),Dam=c(NA,NA,rep("b",6)))
  d <- data.frame(ID=letters[3:8],w=c(10,12,9,15,13,16),Herd=rep(c("H1","H2"),3))
  job <- prepare_blupf90(p,d,"w",random=list(maternal=pedigree_effect("Dam"),herd=iid_effect("Herd")))
  expect_equal(length(job$random_effects),3L)
  expect_equal(job$random_effects$maternal$map$ID,job$id_map$ID)
  expect_equal(unname(job$data[,"maternal"]),rep(match("b",job$id_map$ID),6))
  expect_equal(sum(job$parameter_lines=="RANDOM_GROUP"),3L)
  expect_true(any(job$parameter_lines=="diagonal"))
  expect_error(prepare_blupf90(p,d,"w",random=list(bad=iid_effect("absent"))),"column")
})
test_that("platform specification chooses native executables", {
  expect_match(blup_platform("Linux","x86_64")$url,"Linux/64bit/blupf90")
  expect_match(blup_platform("Windows","AMD64")$filename,"exe$")
  expect_match(blup_platform("Darwin","x86_64")$url,"Mac")
  expect_error(blup_platform("Linux","aarch64"),"architecture")
  expect_error(blup_platform("Darwin","arm64"),"Intel|Rosetta")
})

test_that("selected predictor mappings handle disjoint targets and unconditional rows", {
  set.seed(222)
  d <- data.frame(x=rnorm(120),y=rnorm(120),z=rnorm(120),a=rnorm(120))
  d$x[1:8]<-NA; d$y[9:16]<-NA; d$a[1:3]<-NA
  r <- sequential_em(d,c("x","y","z"),predictor_traits=list(x="a",y="z",z=character()))
  expect_false(anyNA(r$data[c("x","y","z")]))
  expect_equal(r$data$x[1:3],rep(r$target_models$x$distribution_mean["x"],3),ignore_attr=TRUE)
  expect_equal(r$diagnostics$predictor_traits,c("a","z","intercept_only"))
  expect_error(sequential_em(d,c("x","y"),predictor_traits=list(x="a")),"named")
  expect_error(sequential_em(d,c("x","y"),predictor_traits="a",method="sequential"),"requires")
})

test_that("maternal and environmental predictions match exact mixed model equations", {
  exe <- Sys.getenv("BLUPF90_PLUS",unset="")
  skip_if(!nzchar(exe),"Set BLUPF90_PLUS for native mixed-model validation")
  p <- data.frame(ID=as.character(1:12),Sire=c(NA,NA,rep("1",10)),Dam=c(NA,NA,rep("2",4),rep("3",3),rep("4",3)))
  d <- data.frame(ID=as.character(5:12),w=c(10,12,9,15,13,16,8,11),Herd=rep(c("H1","H2"),4))
  job <- prepare_blupf90(p,d,"w",random=list(maternal=pedigree_effect("Dam"),herd=iid_effect("Herd")))
  path <- tempfile(); on.exit(unlink(path,recursive=TRUE),add=TRUE)
  write_blupf90_job(job,path)
  writeLines(blup_parameter_lines(job,matrix(2),matrix(3),FALSE,list(maternal=matrix(4),herd=matrix(1))),file.path(path,"prediction.par"))
  run_native_blup(resolve_blup_executable(exe),"prediction.par",path,"prediction",120,1L)
  A <- pedigree_relationship(p)$A
  Z <- lapply(job$random_effects,function(r) { M<-matrix(0,nrow(d),nrow(r$map)); M[cbind(seq_len(nrow(d)),r$codes)]<-1; M })
  C <- list(animal=2*A,maternal=4*A,herd=diag(2))
  V <- diag(3,nrow(d))
  for(name in names(Z)) V <- V+Z[[name]]%*%C[[name]]%*%t(Z[[name]])
  Vi <- solve(V); X <- matrix(1,nrow(d),1)
  P <- Vi - Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi
  for(name in names(Z)) {
    B <- C[[name]]%*%t(Z[[name]])
    expected <- drop(B%*%P%*%d$w)
    pev <- diag(C[[name]])-diag(B%*%P%*%t(B))
    native <- parse_effect_solutions(job,file.path(path,"solutions"),job$random_effects[[name]])
    expect_equal(native$EBV,unname(expected),tolerance=1e-6)
    expect_equal(native$SE,unname(sqrt(pev)),tolerance=1e-6)
  }
})

test_that("binary validation rejects mismatched operating systems", {
  path<-tempfile(); on.exit(unlink(path),add=TRUE)
  pe<-as.raw(rep(0,512));pe[1:2]<-charToRaw("MZ");pe[61]<-as.raw(128)
  pe[129:134]<-as.raw(c(80,69,0,0,100,134));writeBin(pe,path)
  expect_true(validate_native_binary(path,"Windows"))
  expect_error(validate_native_binary(path,"Linux"),"compatible")
  writeBin(as.raw(c(207,250,237,254,7,0,0,1)),path)
  expect_true(validate_native_binary(path,"Darwin"))
  expect_match(blup_platform("Darwin","arm64",TRUE)$url,"Mac_OSX")
  expect_error(blup_platform("Other","x86_64"),"operating")
})

test_that("real REML extracts all maternal components and uses total variance for h2", {
  exe<-Sys.getenv("BLUPF90_PLUS",unset="")
  skip_if(!nzchar(exe),"Set BLUPF90_PLUS for maternal REML integration")
  set.seed(512); n<-480
  p<-data.frame(ID=as.character(1:(n+50)),Sire=NA_character_,Dam=NA_character_)
  p$Sire[51:nrow(p)]<-sample(as.character(1:10),n,TRUE)
  p$Dam[51:nrow(p)]<-sample(as.character(11:50),n,TRUE)
  A<-pedigree_relationship(p,max_animals=1000)$A
  u<-drop(t(chol(A))%*%rnorm(nrow(A),sd=sqrt(3)))
  m<-drop(t(chol(A))%*%rnorm(nrow(A),sd=2))
  pe<-rnorm(40,sd=sqrt(3))
  d<-data.frame(ID=p$ID[51:nrow(p)])
  dam<-p$Dam[51:nrow(p)]
  d$w<-20+u[match(d$ID,rownames(A))]+m[match(dam,rownames(A))]+pe[match(dam,as.character(11:50))]+rnorm(n,sd=2)
  r<-run_genetic_analysis(p,d,"w",impute=FALSE,random=list(maternal=pedigree_effect("Dam"),mpe=iid_effect("Dam")),blupf90_control=list(executable=exe))
  total<-sum(vapply(r$random_covariances,function(M)M[1,1],numeric(1)))+r$residual_covariance[1,1]
  expect_equal(unname(r$heritability$h2),unname(r$genetic_covariance[1,1]/total))
  expect_equal(nrow(r$random_solutions$maternal),nrow(p))
  expect_equal(nrow(r$random_solutions$mpe),40L)
  expect_equal(nrow(r$variance_ratios),2L)
  expect_true(all(is.finite(r$variance_ratios$SE)))
  expect_true(all(is.finite(r$random_solutions$maternal$Accuracy)))
})

test_that("multi-trait maternal covariance and predictions follow trait order", {
  exe<-Sys.getenv("BLUPF90_PLUS",unset="")
  skip_if(!nzchar(exe),"Set BLUPF90_PLUS for multivariate maternal prediction")
  p<-data.frame(ID=as.character(1:12),Sire=c(NA,NA,rep("1",10)),Dam=c(NA,NA,rep("2",4),rep("3",3),rep("4",3)))
  d<-data.frame(ID=as.character(5:12),w=c(10,12,9,15,13,16,8,11),z=c(3,8,2,5,9,3,7,1))
  job<-prepare_blupf90(p,d,c("z","w"),random=list(maternal=pedigree_effect("Dam")))
  path<-tempfile();on.exit(unlink(path,recursive=TRUE),add=TRUE)
  write_blupf90_job(job,path)
  G<-matrix(c(2,.4,.4,3),2);M<-matrix(c(4,.3,.3,2),2);E<-diag(c(3,2))
  writeLines(blup_parameter_lines(job,G,E,FALSE,list(maternal=M)),file.path(path,"prediction.par"))
  run_native_blup(resolve_blup_executable(exe),"prediction.par",path,"prediction",120,1L)
  A<-pedigree_relationship(p)$A
  Z<-lapply(job$random_effects,function(r){v<-matrix(0,nrow(d),nrow(r$map));v[cbind(1:nrow(d),r$codes)]<-1;kronecker(diag(2),v)})
  C<-list(animal=kronecker(G,A),maternal=kronecker(M,A))
  V<-kronecker(E,diag(nrow(d)))
  for(name in names(Z)) V<-V+Z[[name]]%*%C[[name]]%*%t(Z[[name]])
  Vi<-solve(V);X<-kronecker(diag(2),matrix(1,nrow(d),1))
  P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi
  for(name in names(Z)) {
    B<-C[[name]]%*%t(Z[[name]])
    native<-parse_effect_solutions(job,file.path(path,"solutions"),job$random_effects[[name]])
    expect_equal(native$EBV,unname(drop(B%*%P%*%unlist(d[c("z","w")],use.names=FALSE))),tolerance=1e-6)
    expect_equal(native$SE,unname(sqrt(diag(C[[name]])-diag(B%*%P%*%t(B)))),tolerance=1e-6)
  }
})
