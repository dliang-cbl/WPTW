## use prediction to develop the smooth estimates
## ver 5d) fix non-focused variable at average condition,
## rather than zero them out in the following prediction
## script.
plot_sdm_bc <- function(object,select,n=100,rug=T,ylim=NULL,...){
  ## object: output from a fitted object
  ## select: specie distribution term to plot
  ## n: number of x values to output
  ## rug: add rug on x-axis
  ## ylim: override default limits
  ## output: approximate smooth plot through prediction
  term_ <- all.vars(object$ctds$formula)[-1][select]
  x0_ <- object$ctds$model[,term_]
  if(term_ %in% names(object$ppm$model)){
    tmp_ <- unique(c(all.vars(object$.args$formula),
                     all.vars(object$.args$bias)))
    X_ <- object$ppm$model[,c(tmp_,"(weights)")]
    x_ <- X_[,term_]
  }else{
    X_ <- object$.args$quad
    x_ <- X_[,term_]
  }
  
  ## develop average conditions
  allCovar <- all.vars(object$ctds$formula)[-1][-select]
  for(i in 1:length(allCovar)){
    
    y_ <- X_[,allCovar[i]]
    if(class(y_) %in% c("integer","numeric")){
      cat("averaging ",allCovar[i],"\n")
      X_[,allCovar[i]] <- mean(y_,na.rm=T)
    }else{
      cat("first level ",allCovar[i],"\n")
      X_[,allCovar[i]] <- levels(y_)[1]
    }
  }
  stopifnot(class(x_) %in% c("integer","numeric","factor"))
  
  if(is.numeric(x_)){
    ## discretize x to plot
    n_ <- min(n,length(x_))
    X_$z_ <- as.integer(n_*(x_-min(x_))/(max(x_)-min(x_)))
    X_ <- X_[order(x_),]
    X_ <- subset(X_,!duplicated(z_))
  }
  else{
    ## factor
    X_ <- subset(X_,!duplicated(x_))
  }
  
  ## make prediction
  pred_ <- predict_sdm_bc(object,newd=X_,type="link",
                          terms=term_)
  
  ## form confidence bands
  lwr_ <- pred_$fit - 2*pred_$sd
  fit_ <- pred_$fit
  upr_ <- pred_$fit + 2*pred_$sd
  
  if(is.numeric(x_)){
    ## linear interpolation
    xout <- seq(min(x_),max(x_),len=2*n)
    lwr <- approx(x=X_[,term_],y=lwr_,xout=xout)$y
    fit <- approx(x=X_[,term_],y=fit_,xout=xout)$y
    upr <- approx(x=X_[,term_],y=upr_,xout=xout)$y
    
    ## ecdf approximation and labeling
    label_ <- "Partial Effect"
    if(object$.args$grouped){
      sterms0 <- lapply(object$ctds$smooth,function(lst) lst$term)
      sterms <- sapply(sterms0,paste0,collapse="_")
    }else{
      sterms <- sapply(object$ctds$smooth,function(lst) lst$term)
    }
    if(!is.null(sterms)){
      if(term_ %in% sterms){
        which_ <- which(sterms == term_)
        first_ <- object$ctds$smooth[[which_]]$first.para
        last_ <- object$ctds$smooth[[which_]]$last.para
        edf_ <- sum(object$ctds$edf[first_:last_])
        label_ <- paste0("s(",term_,",",round(edf_,2),")")
      }
    }
    
    poly_ <- rbind(cbind(xout,lwr),cbind(rev(xout),rev(upr)))
    value <- list(x=xout,y=fit,band=poly_,
                  xlab=term_,
                  ylab=label_,
                  ylim=range(c(lwr,upr)))
    
    if(!is.null(ylim)){
      value$ylim <- ylim
    }
    ## plot
    with(value,plot(x,y,ylim=ylim,type="n",
                    xlab=xlab,ylab=ylab,...))
    if(rug){
      rug(unique(x0_))
      rug(X_[,term_],col=2)
    }
    with(value,polygon(band,col=grey(0.50)))
    with(value,lines(x,y,col=grey(1),lwd=2))
    
  }else{
    label_ <- "Partial Effect"
    
    value <- data.frame(X_[,term_],lwr_,fit_,upr_)
    names(value) <- c(term_,"lwr","fit","upr")
    
    ## plot
    plot_sdm_bc_draw(value[,1],value[,3],value[,2],value[,4],
                     xlab=xlab,ylab=ylab,...)
  }
  
  value
}

plot_sdm_bc_draw <- function(x,fit,lwr,upr,...){
  xl <- seq(1,length(x))
  
  yr <- range(c(fit,lwr,upr))
  plot(c(1,max(xl)+0.5*length(x)),yr,type="n",axes=F,...)
  currx <- 0
  for(i in 1:length(x)){
    poly_ <- cbind(currx+c(xl[i],xl[i]+1,xl[i]+1,xl[i]),
                   c(lwr[i],lwr[i],upr[i],upr[i]))
    polygon(poly_,col = grey(0.7))
    line_ <- poly_[1:2,]
    line_[,2] <- fit[i]
    lines(line_,col="white",lwd=2)
    currx <- currx + 0.5
  }
  axis(side=2,at=pretty(c(fit,lwr,upr)))
  axis(side=1,at=xl*1.5,labels = levels(x))  
}
dev_ <- function(){
  rm(list=ls())
  library(mgcv)
  source("term_utils.R")
  load(file="Scratch/sdm_bc_dev.rda")
  source("sim2jam_study_1.R")
  source("predict_gam_term_2.R")
  source("all_group_1.R")
  library(abind)
  ssh <- c(quad$ssh,glmobj$ssh)
  ssh <- cut(ssh,breaks=quantile(ssh,na.rm=T),labels = F)
  quad$ssh2 <- factor(ssh[seq(1,nrow(quad))])
  glmobj$ssh2 <- factor(ssh[nrow(quad)+seq(1,nrow(glmobj))])
  
  formula=z~bathy+s(sst,bs="cr",k=3)+fpi+ssh2
  term=c(1,2,3,3)
  bias=~s(long,lat,bs="gp",m=c(3,0.1),k=30)
  #bias=~s(bathy,bs="cr",k=3)+fpi+ssh
  timefact=30; verbose=T; iter=9; gcoord = "id"
  wcoord="wt"; tcoord="tau"; cterm="crw";trace=T
  source("sdm_bc_4b.R")
  if(FALSE){
    object <- sdm.bc(
      formula,term,bias,quad,glmobj,
      wcoord,gcoord,tcoord,cterm,timefact,
      iter,verbose,trace)
    save(object,file='Scratch/sdm_bc_dev_object_2d.rda')
  }else{
    load(file='Scratch/sdm_bc_dev_object_2d.rda')
  }
  
  source("predict_sdm_bc_4f.R")
  select <- 1;n=100;rug=T
  source("plot_sdm_bc_2d.R")
  summary(object$ctds)
  debug(plot_sdm_bc)
  dummy <- plot_sdm_bc(object,select=4)
  plot(object$ctds,select=1,residual=FALSE,scale=0)
  #debug(plot_sdm_bc)
  dummy <- plot_sdm_bc(object,select=2)
  #debug(plot_sdm_bc)
  dummy <- plot_sdm_bc(object,select=3)
  dummy <- plot_sdm_bc(object,select=1)
  
}