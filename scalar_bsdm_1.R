library(bayesplot)
scalar.bsdm <- function(bsdm.obj,perm=NULL){
  ## perm: permutation to list the terms
  names(bsdm.obj$value$theta[[1]]$coef)
  label_ <- sapply(bsdm.obj$.args$term2,paste0,collapse="=")
  eterms_ <- grep("E\\.",names(bsdm.obj$value$theta[[1]]$coef))
  mcmc_x <- do.call(rbind,lapply(
    bsdm.obj$value$theta,function(x){
      x$b[eterms_]
    }
  ))
  colnames(mcmc_x) <- label_
  mcmc_s <- t(apply(
    mcmc_x,2,quantile,prob=c(0.5,0.025,0.975)))
  colnames(mcmc_s) <- c("fit","lwr","upr")
  if(is.null(perm)){
    print(mcmc_dens(mcmc_x))
  }else{
    print(mcmc_dens(mcmc_x[,perm]))
  }
  if(is.null(perm)){
    return(mcmc_s)
  }else{
    return(mcmc_s[perm,])
  }
}