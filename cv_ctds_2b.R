require(doParallel)
## cross-validation of movement by individuals
cv.ctds <- function(x,data,zcoord="z", tcoord="t",
                    icoord="id",taucoord="tau",
                    type=1,fold=10,r=1,plotit=T){
  ## fitted object x: from joint modeling
  ## data: track data
  ## type: cross-validate  type 1 across individuals
  ##     : type 2: within individual running over days.
  ## fold: for type =2, divide the time periods into
  ##        multiple folds and move forward in time
  ##        use the first half as training set
  ## coord: names for movement, time, individual and residence
  ## binary movement response
  ## r: relative cost for false negative.
  ## plotit: default plot the ROC curve
  
  ## remove NAs
  tmp_ <- na.omit(data[,all.vars(x$formula)])
  na_action <- attr(tmp_,"na.action")
  if(length(na_action)>0){
    data <- data[-na_action,]
  }
  
  ## cross-validation storage of fitted values
  if(type==1){
    ## cross-validation by individuals
    batchid_ <-  unique(data[,icoord])
    cat("cross-validate by individuals...\n")
    parout <- foreach(i = 1:length(batchid_)) %dopar% {
      cat(i,",")
      sel_ <- data[,icoord] == batchid_[i] ## test set
      train_ <- data[!sel_,] ## training
      test_ <- data[sel_,]  ## testing
      ## train model
      ## build model using  data
      ctds.gam <- gam(x$formula,data=train_,
                      family=quasi(link="log",variance="mu"))
      pred <- predict(ctds.gam,newdata=test_,type="response")
      list(idx=sel_,val=pred)
    }
    cat("done.\n")
  }
  
  if(type==2){
    cat("cross-validation by time within individuals.\n")
    
    lst_ <- split(data,~id)
    for(i in 1:length(lst_)){
      nstep_ <- nrow(lst_[[i]])
      fold_ <- min(nstep_/2,fold)
      id_ <- sort(rep_len(1:(2*fold_),nstep_))
      id2_ <- pmax(0,id_-fold_)
      lst_[[i]]$batch <- id2_
    }
    
    data_ <- do.call(rbind,lst_)
    
    parout <- foreach(i =1:fold) %dopar% {
      cat(i,",")
      sel_ <- data_$batch < i ## training set
      psel_ <- data_$batch == i ## test set
      train_ <- data_[sel_,] ## training
      test_ <- data_[psel_,]  ## testing
      ## train model
      ## build model using  data
      ctds.gam <- gam(x$formula,data=train_,
                      family=quasi(link="log",variance="mu"))
      pred_ <- predict(ctds.gam,newdata=test_,type="response")
      list(idx=psel_,val=pred_)
    }
    #browser()
  }
  
  ## actual data
  z <- data[,all.vars(x$formula)[1]]
  
  ## assemble the prediction
  fit <- rep(NA,nrow(data))
  for(i in 1:length(parout)){
    fit[parout[[i]]$idx] <- parout[[i]]$val
  }
  ## CV results optimal threshold
  if(type==1){
    val <- cv.ctds.work(fit,z,r=r)
  }else{
    ## only extract the second half for validation
    val <- cv.ctds.work(fit[!is.na(fit)],z[!is.na(fit)],r=r)
  }
  
  if(plotit){
    plotCostROC(val,which="roc")
  }
  cat("Sensitivity=",val$T$sens,", Specificity=",val$T$spec,"\n")
  
  val
}

cv.ctds.work <- function(x,y,sample=T,r=1){
  ## y: binary status (healthy=0, disease=1)
  ## x: probability
  ## r: cost of mis-classification
  ## result classification threshold
  k1 <- x[y==0]+1e-8
  k2 <- x[y==1]+1e-8
  rho <- mean(y)
  
  if(sample){
    k1 <- k1[sample(length(k1),length(k2))]
  }
  costs <-  matrix(c(0, 0, r, (1 - rho)/rho), 
                   2, 2, byrow = TRUE)
  
  thres2(k1,k2,rho,costs=costs,
         method="empirical",ci=FALSE,extra.info = T)
  
}
