## bias corrected specie distribution models
## pooling fishery data and telemetry observations
## ver. 5d) allow informative priors on the scaling coefficients.
## ver. 4) multiple imputation of cut inference using mgcv
##         with joint samples stored
## sub ver. b) allow random effects to incorporate
##             within turtle correlation
sdm.bc <- function(formula,term,bias, quad, glmobj,
                   wcoord, gcoord, tcoord, cterm,
                   timefact,d=-1,S=0.125,niter=99,
                   verbose=TRUE,trace=FALSE){
  ## formula: a gam formula defining specie distribution model
  ## term   : a vector of integers defining the weighting
  ##          let nTerm number of terms in formula
  ##           J<nTerm the number of weight parameters
  ##          from vector[i] = j, j=1,..,J
  ##          Default j=i, i.e. all terms weighted before
  ##          entering specie distribution modeling
  ## bias   : a one-sided gam formula defining 
  ##          observation bias
  ## quad   : a data frame of fishery observations
  ##          with quadrature points added
  ## glmobj : a data frame of telemetry data converted
  ##          into a glm formulation (Hanks et al. 2015)
  ##
  ## wcoord : variable name for weight in quad
  ##
  ## gcoord : variable name for group in glmobj
  ## tcoord : variable name for residence time in glmobj
  ## ccoord : variable name for crw term in glmobj
  
  ## d      : Normal prior mean for scaling coefficients
  ##          default (-1)
  ## S     : Normal prior standard deviation for scaling coefficients
  ##          default 0.125 such that scale within 25% of -1 with 95% prior probability
  ## timefact: number of time units for the motility coefficient
  ##           in each ppm coefficient
  ##           default is 1,  that is the motility
  ##           coefficients were estimated using the same
  ##           time units as the ppm coefficients
  ##           
  ##           this can change to 30 when ppm coefficients on a 
  ##           different time unit (month) than motility 
  ##           coefficient (day)
  
  ## niter   : number of imputation to approximate the 
  ##          "cut" joint posterior
  ## return  : 
  ##  value: a list of samples from 
  ##      beta: resource selection coefficients
  ##      theta: a list of bias correction coefficients given beta
  ##  ctds : a gam(lm) object of quasi-likelihood fits for CTDS
  ##  ppm  : a gam(lm) object of a quasi-likelihood fit for PPM
  ##         mainly used to structure the predictive results
  ##  dic  : deviance information criterion
  ##  arg  : a list of arguments
  ##        quad : input quadrature data
  
  grouped <- !is.null(gcoord)
  if(grouped) glmobj$id <- glmobj[,gcoord]
  glmobj$tau <- glmobj[,tcoord]
  glmobj$crw <- glmobj[,cterm]
  glmobj$tau <- glmobj$tau / timefact ## report residence time in month

  sdm.terms <- all.vars(formula)
  if(!all(sdm.terms %in% names(glmobj))){
    stop("Some terms in formula not found in glmobj.\n")
  }
  
  ## define specie distribution terms
  sdm.name_ <- c(term.rhs(formula),"crw","offset(log(tau))")
  if(grouped){
    ## formula must contain habitat term only
    ## no bivariate smooth allowed
    re.terms <- paste0("s(id,",var.rhs(formula),",bs='re')")
    for0 <- reformulate(term=c(sdm.name_,"s(id,bs='re')",
                               re.terms),
                        response=var.lhs(formula))
  }else{
    for0 <- reformulate(term=sdm.name_,
                        response=var.lhs(formula))
  }
  ## fit specie distribution model using tracking data
  ctds.gam <- gam(formula=for0,data=glmobj,
                  family=quasi(link="log",variance="mu"),
                  #method="ML",
                  control = gam.control(trace=trace))
  
  ## sample from posterior of resource selection coefficients
  bmu <- coef(ctds.gam)
  bsigma <- vcov(ctds.gam,unconditional=TRUE)
  bsam <- rmvn(niter,bmu,bsigma)
  ## specie distribution terms can't include id
  if(grouped){
    ## pick out the individual specific random effects
    ri.pattern <- "^s\\(id\\)\\.[0-9]+$"
    #bsam[,grep(ri.pattern,names(bmu))] <- 0
    ## pick out the individual specific random coefficients
    rc.pattern <- "^s\\(id\\,[A-Za-z]+\\)\\.[0-9]+$"
    #rc.ii <- grep("^s\\(id\\,[A-Za-z]+\\)\\.[0-9]+$",names(bmu))
    #bsam[,grep(rc.pattern,names(bmu))] <- 0
    group.col <- c(grep(ri.pattern,names(bmu)),
                   grep(rc.pattern,names(bmu)))
  }else{
    group.col <- NULL
  }
  
  ## define multiplicative scaling
  definition_ <- split(var.rhs(formula),term)

  ## pre-processing fishery data
  quad$wt <- quad[,wcoord]
  bias.terms <- all.vars(bias)
  if(!all(c(sdm.terms[-1],bias.terms) %in% names(quad))){
    stop("Some terms in bias not found in quad.\n")
  }
  bias.name_ <- term.rhs(bias)
  
  ## prepare linear predictor of unbiased fishery
  quad$crw <- 0
  quad$tau <- 1
  quad$id <- ctds.gam$model$id[1]
  
  ## remove missing environmental covariates
  l.quad <- na.omit(quad[,c(bias.terms,"wt",all.vars(for0))])
  X <- predict(ctds.gam,l.quad,type="lpmatrix")

  alphaLst <- vector("list",niter)
  for(i in 1:niter){
    ## unbiased and un-scaled prediction of log-intensity
    if(grouped){
      E.lst <- predict.gam.term(
        formula=reformulate(term=sdm.name_,
                            response=var.lhs(formula)),
        x = definition_,
        beta = bsam[i,-group.col],
        X = X[,-group.col,drop=FALSE])
    }else{
      E.lst <- predict.gam.term(
        formula=ctds.gam$formula,
        x = definition_,
        beta = bsam[i,],
        X = X)
    }
    
    ## estimate bias and scales
    for1 <- reformulate(c(names(E.lst),bias.name_),
                        response="z")
    gam.bsdm <- gam(formula=for1,data=cbind(l.quad,E.lst),
                    weights=wt,
                    family=quasi(link="log",variance="mu"))
    
    #browser()
    
    ## simulate from compute posterior distribution of
    # conjugate Normal family
    ## mu0: prior mean
    ## sigma0: prior standard deviation
    ## xbar: observed mean
    ## sigmasqhat: observed covariance matrix
    ## idx: the index to which prior to affect
    ## value: an output
    simPostNorm <- function(mu0,sigma0,xbar,sigmasqhat,idx=NULL){
      m0 <- rep(mu0,length(xbar))
      q0 <- diag(1/sigma0^2,length(xbar))
      if(!is.null(idx)){
        m0[-idx] <- 0
        diag(q0)[-idx] <- 0
      }
      qPost <- solve(sigmasqhat) + q0
      bPost <- solve(sigmasqhat,xbar) + q0 %*% m0
      mPost <- as.vector(solve(qPost,bPost))
      names(mPost) <- names(xbar)
      list(mean=mPost,Cov=solve(qPost),
           value=rmvn(1,mPost,solve(qPost)))
    }
    ## a sample from Gaussian approximation
    tSigma <- vcov(gam.bsdm,unconditional = T)
    #asam <- rmvn(1,coef(gam.bsdm),tSigma)
    #browser()
    tSim <- simPostNorm(
      mu0 = d,sigma0=S,xbar=coef(gam.bsdm),
      sigmasqhat = tSigma,idx=1+1:ncol(E.lst))
    asam <- tSim$value

    if(verbose){
      cat("i=",i,",",asam[1+1:ncol(E.lst)],"\n")
    }
    
    ## fitted log likelihood
    Xp <- predict(gam.bsdm,type="lpmatrix")
    etahat <- Xp %*% tSim$mean
    muhat <- as.numeric(gam.bsdm$family$linkinv(etahat))
    lhat <- with(gam.bsdm, sum(
      model[,"(weights)"]*(y*log(muhat)-muhat)
      #model[,"(weights)"]*(y*log(fitted.values)-fitted.values)
    ))
    
    ## current deviance
    eta <- Xp %*% asam
    mu <- as.numeric(gam.bsdm$family$linkinv(eta))
    llik <- with(gam.bsdm, sum(
      model[,"(weights)"]*(y*log(mu)-mu)
    ))
    
    ## store sample from the current iteration
    sam_ <- list(
      b=asam, #coef=coef(gam.bsdm),Vp=tSigma,
      coef =tSim$mean, vp=tSim$Cov,
      scale=gam.bsdm$scale)
    
    ## return
    alphaLst[[i]] <- list(value=sam_,lhat=lhat,llik=llik)
  }

  helper <- function(object,x){
    arr_ <- lapply(object,function(lst){lst[[x]]})
    dim_ <- dim(arr_[[1]])
    dim_ <- ifelse(is.null(dim_),1,length(dim_))
    arr_ <- abind(arr_,along=0)
    apply(arr_,1+seq(1,dim_),mean)
  }

  args <- list(
    eDefinition =all_group_label(term),
    term=term,
    term2=definition_,
    formula=formula,
    bias = bias,
    bformula=for1,
    grouped=grouped,
    group.col=group.col)
  
  ## dic calculation
  lhat <- helper(alphaLst,"lhat")
  llik <- helper(alphaLst,"llik")
  dic <- data.frame(Deviance=-2*llik,
                    Dhat = -2*lhat)
  dic$pD <- with(dic,Deviance - Dhat)
  dic$DIC <- with(dic,Deviance+pD)
  
  ## samples from theta
  theta_ <- lapply(alphaLst,function(lst) lst$value)
  
  ## combine data from quad to gam.bsdm
  sdm.only.terms__ <- all.vars(formula)[-1]
  sdm.only.terms__ <- sdm.only.terms__[!(sdm.only.terms__ %in% bias.terms)]
  args$quad <- NULL
  if(length(sdm.only.terms__)>0){
    args$quad <- l.quad
  }
  
  list(value=list(theta=theta_,beta=bsam),
       dic=dic,.args=args,ctds=ctds.gam,
       ppm=gam.bsdm)
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
  formula=z~bathy+s(sst,bs="cr",k=3)+fpi+ssh
  term=c(1,2,3,3)
  bias=~s(long,lat,bs="gp",m=c(3,0.1),k=10)
  #bias=~s(bathy,bs="cr",k=3)+fpi+ssh
  timefact=30; verbose=T; niter=9; gcoord = "id"
  wcoord="wt"; tcoord="tau"; cterm="crw";trace=T
  source("sdm_bc_5d.R")
  object <- sdm.bc(
    formula,term,bias,quad,glmobj,
    wcoord,gcoord,tcoord,cterm,timefact,
    niter=iter,verbose=verbose,trace=trace)
  save(object,file='Scratch/sdm_bc_dev_object.rda')
}

