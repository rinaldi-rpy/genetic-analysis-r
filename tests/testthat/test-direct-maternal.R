test_that("correlated maternal effects share a full native covariance group", {
  p<-data.frame(ID=as.character(1:12),Sire=c(NA,NA,rep("1",10)),Dam=c(NA,NA,rep("2",4),rep("3",3),rep("4",3)))
  d<-data.frame(ID=as.character(5:12),w=c(10,12,9,15,13,16,8,11),z=c(3,8,2,5,9,3,7,1))
  job<-prepare_blupf90(p,d,c("z","w"),random=list(maternal=pedigree_effect("Dam",correlate_with="animal")))
  expect_equal(job$direct_maternal_effect,"maternal")
  expect_equal(sum(job$parameter_lines=="RANDOM_GROUP"),1L)
  expect_true(any(job$parameter_lines==paste(job$animal_effect,job$random_effects$maternal$effect)))
  expect_equal(dim(job$direct_maternal_start),c(2L,2L))
  expect_true(all(job$direct_maternal_start!=0))
  expect_true(any(grepl("OPTION se_covar_function ram_1_2",job$parameter_lines,fixed=TRUE)))
  expect_error(prepare_blupf90(p,d,"w",random=list(x=pedigree_effect("Sire",correlate_with="animal"))),"Dam")
})

test_that("correlated maternal native predictions and PEV match exact equations", {
  exe<-Sys.getenv("BLUPF90_PLUS",unset="")
  skip_if(!nzchar(exe),"Set BLUPF90_PLUS for correlated maternal prediction")
  for(k in 1:2) for(sign in c(-1,1)) {
    p<-data.frame(ID=as.character(1:12),Sire=c(NA,NA,rep("1",10)),Dam=c(NA,NA,rep("2",4),rep("3",3),rep("4",3)))
    d<-data.frame(ID=as.character(5:12),w=c(10,12,9,15,13,16,8,11),z=c(3,8,2,5,9,3,7,1))
    traits<-c("z","w")[seq_len(k)]
    job<-prepare_blupf90(p,d,traits,random=list(maternal=pedigree_effect("Dam",correlate_with="animal")))
    path<-tempfile();on.exit(unlink(path,recursive=TRUE),add=TRUE)
    write_blupf90_job(job,path)
    G<-matrix(c(2,.4,.4,3),2)[seq_len(k),seq_len(k),drop=FALSE]
    M<-matrix(c(4,.3,.3,2),2)[seq_len(k),seq_len(k),drop=FALSE]
    C<-sign*matrix(c(.4,.2,-.1,.3),2)[seq_len(k),seq_len(k),drop=FALSE]
    E<-diag(c(3,2))[seq_len(k),seq_len(k),drop=FALSE]
    writeLines(blup_parameter_lines(job,G,E,FALSE,list(maternal=M),C),file.path(path,"prediction.par"))
    run_native_blup(resolve_blup_executable(exe),"prediction.par",path,"prediction",120,1L)
    A<-pedigree_relationship(p)$A
    Z<-lapply(job$random_effects,function(r){v<-matrix(0,nrow(d),nrow(r$map));v[cbind(1:nrow(d),r$codes)]<-1;kronecker(diag(k),v)})
    Z<-do.call(cbind,Z)
    prior<-kronecker(direct_maternal_block(G,M,C),A)
    V<-Z%*%prior%*%t(Z)+kronecker(E,diag(nrow(d)))
    Vi<-solve(V);X<-kronecker(diag(k),matrix(1,nrow(d),1))
    P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi
    B<-prior%*%t(Z)
    expected<-drop(B%*%P%*%unlist(d[traits],use.names=FALSE))
    se<-sqrt(diag(prior)-diag(B%*%P%*%t(B)))
    native<-lapply(job$random_effects,function(r)parse_effect_solutions(job,file.path(path,"solutions"),r))
    expect_equal(unlist(lapply(native,`[[`,"EBV"),use.names=FALSE),unname(expected),tolerance=1e-6)
    expect_equal(unlist(lapply(native,`[[`,"SE"),use.names=FALSE),unname(se),tolerance=1e-6)
  }
})

test_that("real two-trait REML estimates and reports direct-maternal covariances", {
  exe<-Sys.getenv("BLUPF90_PLUS",unset="")
  skip_if(!nzchar(exe),"Set BLUPF90_PLUS for correlated maternal REML")
  set.seed(504); n<-600L
  p<-data.frame(ID=as.character(1:660),Sire=NA_character_,Dam=NA_character_)
  for(g in 1:6) {
    previous<-if(g==1) 1:60 else (61+(g-2)*100):(60+(g-1)*100)
    rows<-(61+(g-1)*100):(60+g*100)
    p$Sire[rows]<-sample(p$ID[previous[1:15]],100,TRUE)
    p$Dam[rows]<-rep(p$ID[previous[16:40]],each=4)
  }
  rel<-pedigree_relationship(p,max_animals=1000)
  G<-matrix(c(4,.7,.7,3),2);M<-matrix(c(3,.2,.2,2),2)
  C<-matrix(c(-.8,.2,-.1,-.5),2)
  u<-t(chol(direct_maternal_block(G,M,C)))%*%matrix(rnorm(4*nrow(p)),4)%*%chol(rel$A)
  rows<-61:nrow(p);d<-data.frame(ID=p$ID[rows])
  for(j in 1:2) d[[c("w","z")[j]]]<-20*j+u[j,match(d$ID,rownames(rel$A))]+u[j+2,match(p$Dam[rows],rownames(rel$A))]+rnorm(n,sd=sqrt(5))
  r<-run_genetic_analysis(p,d,c("w","z"),impute=FALSE,
    random=list(maternal=pedigree_effect("Dam",correlate_with="animal")),blupf90_control=list(executable=exe))
  expect_equal(dim(r$direct_maternal_covariance),c(2L,2L))
  expect_equal(dim(r$joint_genetic_covariance),c(4L,4L))
  expect_equal(nrow(r$direct_maternal_parameters),4L)
  expect_true(all(is.finite(r$direct_maternal_parameters$covariance_SE)))
  expect_true(all(is.finite(r$direct_maternal_parameters$correlation_SE)))
  expect_true(all(abs(r$direct_maternal_parameters$correlation)<=1))
  expect_true(any(abs(r$direct_maternal_covariance)>1e-4))
  total<-diag(r$genetic_covariance)+diag(r$random_covariances$maternal)+diag(r$residual_covariance)+diag(r$direct_maternal_covariance)
  expect_equal(r$heritability$h2,unname(diag(r$genetic_covariance)/total))
  expect_equal(nrow(r$random_solutions$maternal),nrow(p)*2)
  expect_true(r$diagnostics$direct_maternal_correlated)
})

test_that("joint covariance validation rejects inadmissible cross-effect blocks", {
  expect_error(direct_maternal_block(matrix(1),matrix(1),matrix(2)),"admissible")
  expect_error(direct_maternal_block(matrix(1),diag(2),matrix(1)),"dimensions")
  expect_error(pedigree_effect("Dam",correlate_with="other"),"correlate_with")
  p<-data.frame(ID=as.character(1:6),Sire=c(NA,NA,rep("1",4)),Dam=c(NA,NA,rep("2",4)))
  d<-data.frame(ID=as.character(3:6),w=c(2,5,3,7),Herd=c("a","a","b","b"))
  expect_error(prepare_blupf90(p,d,"w",random=list(m1=pedigree_effect("Dam","animal"),m2=pedigree_effect("Dam","animal"))),"one maternal")
  job<-prepare_blupf90(p,d,"w",random=list(herd=iid_effect("Herd"),maternal=pedigree_effect("Dam","animal")))
  expect_equal(job$random_effects$maternal$effect,job$animal_effect+1L)
  expect_error(blup_parameter_lines(job,matrix(1),matrix(1),FALSE,list(maternal=matrix(1),herd=matrix(1))),"converged")
})
