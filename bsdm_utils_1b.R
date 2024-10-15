ppm.resid <- function(x){
  model_ <- x$ppm$model
  names(model_)[ncol(model_)] <- "wt"
  e <- ppm.st.residuals(bsdm.obj$ppm,model_)
  z <- e/sqrt(model_$wt)
  cbind(model_,resid=e,rstand=z)
}

#x <- bsdm.obj
#data <- pa
#tcoord <- "amon"
# wcoord <- "wt_scale"
## cross-validate the fitted joint model using
## fishery component by month
ppm.cv <- function(x,data,tcoord="amon",wcoord = "wt_scale"){
  ## divide the data by month
  ## using half of the data for validation purposes
  amon_ <- sort(unique(data[,tcoord]))
  batchs_ <- amon_[seq(floor(length(amon_)/2),length(amon_))]
  
  ppm.obj <- x$ppm ## ppm object
  cat("cross-validate by time..,\n")
  #browser()
  storage <- foreach(i=1:length(batchs_)) %dopar%{
    cat("month ",batchs_[i],"\n")
    start_ <- amon_[i]
    sel_ <- (data[,tcoord]<batchs_[i]) & (data[,tcoord]>=start_)
    psel_ <- data[,tcoord]==batchs_[i]
    train_ <- ppm.obj$model[sel_,]
    test_ <- data[psel_,]
    names(train_)[ncol(train_)] <- "wt"
    test_$wt <- test_[,wcoord]
    #names(test_)[ncol(test_)] <- "wt"
    
    ## update the ppm model with only part of the data
    x$ppm <- gam(ppm.obj$formula,data=train_,
                      weights=wt,
                      family=quasi(link="log",variance="mu"))
    
    ## predict
    pred_ <- predict_sdm_bc_work(x,newd = test_)
    pred_$resid <- with(pred_,{
      wt*(z-fit)*1/sqrt(fit)
    })
    
    pred_
  }
  
  do.call(rbind,storage)
}

debug_ <- function(){
  rm(list=ls())
  library(mgcv)
  load(file="easi_fish/easifish.rda")
  load(file="easi_fish/bsdm5A3T203B1.rda")
  
  source("predict_sdm_bc_5d.R")
  source("predict_gam_term_2.R")
  source("term_utils.R")
  
  x <- bsdm.obj; data <- pa3;
  tcoord="imonth";wcoord = "wt_scale"
}