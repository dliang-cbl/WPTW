## ctmc2glm.study: helper function to query dynamic environmental covariates
## query assuming a static environmental covariate, then use the resulting
## coordinates to query to list of stack to build dynamic covariates
## D. Liang

###
## NOTE: assumed that tmap , temporal resolution is the SAME for multiple
## environmental covariates. 
## if the environmental covariate is static, use a list of raster layer instead
## then the tmap argument will be ignored.
## ASSUME:
## assume all individual rasters must have the same resolution.
## assume the environmental covariate does not vary rapidly 
##        as the animal moves
## or the raster resolution is not too coarse with no movement possible.

## the gradient support from dynamic covariate is disabled for NOW!
## version 1: start this model with a PPM model

ctmc2glm.study <- function(xy,t,tmap,stack.static, stack.grad, verbose=0,...) 
{
  #ARGUMENTS:
  #xy	A matrix of x,y locations at T time points.
  #t	  A vector of T times associated with the T locations in "xy".
  #tmap A vector of T periods associated with the T locations in "xy" and stacks of environmental data
  #stack.static A list of rasterStack objects, 
  #  where each stack is a location-based covariate and
  #  where each layer in the stack is a snap shot of the location-based covariate.
  #stack.grad A list of rasterStack objects,
  #  where each stack is a directional gradient-based covariate
  #  where each layer in the stack is a snap shot of the gradient-based covariate
  # ... : additional arguments to ctmc2glm
  #VALUES:
  # z: Response variable (either zero or 1) for analysis using GLM software
  # X: matrix of predictor variables for analysis using GLM software. Created
  #    from the location-based and gradient-based covariates
  # tau: Offset for each row in a Poisson GLM with log link.
  # t: Vector of the time each raster grid cell was entered.
  
  ## check raster input of same dimension
  n1 <- "null"
  if(length(stack.static)>0){
    n1 <- unique(sapply(stack.static,class))
    if(length(n1)>1){
      stop("stack.static must be of same temporal dimension.\n")
    }
    if(n1=="RasterLayer"){
      warning("stack.static assumed static, tmap ignored\n")
    }
    if(n1=="RasterStack"){
      v1 <- sapply(stack.static,function(elmt) dim(elmt)[3])
      if(length(table(v1))>1){
        stop("stack.static must contain stack of same temporal dimension.\n")
      }
    }
  }
  if(n1=="null"){
    stop("statck.static must not be NULL\n.")
  }
  n2 <- "null"
  if(length(stack.grad)>0){
    n2 <- unique(sapply(stack.grad,class))
    if(length(n2)>1){
      stop("stack.grad must be of same temporal resolution\n")
    }
    if(n2=="RasterLayer"){
      warning("stack.grad static, tmap ignored\n")
    }
    if(n2=="RasterStack"){
      v2 <- sapply(stack.grad,function(elmt) dim(elmt)[3])
      if(length(table(v2))>1){
        stop("stack.grad must contain stack of same temporal dimension.\n")
      }
    }
  }
  
  ## only static environmental covariates
  if(n1 !="RasterStack" & n2 !="RasterStack"){
    # query static environmental covariates
    if(n1=="RasterLayer"){
      if(length(stack.static)>1){ ## more than one covariates
        l.stack.static <- stack(stack.static) ## return a stack
        examplerast <- raster(l.stack.static,layer=1)
      }
      else{ ## only one covariates
        l.stack.static <- stack.static[[1]] ## return a raster
        examplerast <- l.stack.static
      }
    }
    # query gradient environmental covariates
    if(n2=="null"){
      l.stack.grad <- NULL
    }
    else if(n2=="RasterLayer"){
      if(length(stack.grad)>1){
        l.stack.grad <- stack(stack.grad)
        examplerast <- raster(l.stack.grad,layer=1)
      }
      else{
        l.stack.grad <- stack.grad[[1]]
        examplerast <- l.stack.grad
      }
    }
    ctmc <- path2ctmc(xy=xy,t=t,rast=examplerast,print.iter=FALSE)
    glm.data <- ctmc2glm.study.work(ctmc,l.stack.static,l.stack.grad,...)
  }
  if(n1 =="RasterStack" | n2 =="RasterStack"){
    #browser()
    ## query as if the environmental covariates are static;
    # query static environmental covariates
    static.lst <- ctmc2glm.st.build(n1,stack.static)
    # query gradient environmental covariates
    grad.lst <- ctmc2glm.st.build(n2,stack.grad)
    
    ## exit if track too short
    if(nrow(xy)<3 || length(unique(t))<3){
      return(NULL)
    }    
    ## build glm.data assuming static covariates
    ctmc <- path2ctmc(xy=xy,t=t,rast=static.lst$example,print.iter=FALSE)
    
    ## exit if track too short
    if(length(ctmc$ec)<2){
      return(NULL)
    }
    
    glm.data <- ctmc2glm(ctmc,static.lst$stack,grad.lst$stack,...)
    ##
    
    
    ## identify the layer of environmental covariates at each stay/move
    glm.data.tidx <- ctmc2glm.knn(dataTime=t,queryTime=glm.data[,"t"],tmax=NULL)
    glm.data$tmap__ <- tmap[glm.data.tidx]
    ##
    
    ## split the glm data according to the time animal entering the grid cell
    glm.data.list <- split(x=glm.data,f=glm.data[,"tmap__"])
    ## each entry corresponds to a different layer in the list
    
    ## identify the static and gradient covariates if any
    static.id <- 1+seq(1,static.lst$nvar)
    if(grad.lst$nvar>0) {
      grad.id <- max(static.id)+seq(1,grad.lst$nvar)
    }else{
      grad.id <- NULL
    }
    if(verbose>3){
      cat("static covariates =",static.id," ")
      if(grad.lst$nvar>0){
        cat(" gradient covariate(s) = ",grad.id," ")
      }
    }
    ## extract environmental covariates at each group of environmental condition
    out <- vector("list",length(glm.data.list))
    for(i in 1:length(out)){
      out[[i]] <- glm.data.list[[i]] ## initialize output
      layerid <- glm.data.list[[i]][1,"tmap__"] ## define layerid
      if(verbose>3){
        cat("querying layer ",layerid," ")
      }
      #xy <- as.matrix(glm.data.list[[i]][,c("x.current","y.current")])
      xy <- as.matrix(glm.data.list[[i]][,c("x.adj","y.adj")])
      
      ## define static and dynamic layers
      l.static.lst <- ctmc2glm.st.build(n1,stack.static,layerid=layerid) 
      l.grad.lst <- ctmc2glm.st.build(n2,stack.grad,layerid=layerid)
      #browser()
      
      out[[i]][,static.id] <- extract(l.static.lst$stack,xy)
      if(!is.null(grad.id)){
        stop("gradient not supported\n")
        out[[i]][,grad.id] <- extract(l.grad.lst$stack,xy)
      }
      if(verbose>3){
        cat("done\n")
      }
    }
    
    ## combine the output and generate output
    glm.data <- do.call(rbind,out)
  }
  glm.data
}

ctmc2glm.study.work <- function (ctmc, stack.static, stack.grad, crw = TRUE, normalize.gradients = FALSE, 
          grad.point.decreasing = TRUE, include.cell.locations = TRUE, 
          directions = 4, zero.idx = integer()) 
{
  p.static = nlayers(stack.static)
  p.crw = 0
  if (crw) {
    p.crw = 1
  }
  if (class(stack.grad) == "RasterLayer" | class(stack.grad) == 
      "RasterStack") {
    p.grad = nlayers(stack.grad)
    stack.gradient = rast.grad(stack.grad)
    if (normalize.gradients) {
      lengths = sqrt(stack.gradient$grad.x^2 + stack.gradient$grad.y^2)
      stack.gradient$grad.x <- stack.gradient$grad.x/lengths
      stack.gradient$grad.y <- stack.gradient$grad.y/lengths
    }
  }
  else {
    p.grad = 0
  }
  p = p.static + p.crw + p.grad
  if (class(stack.static) == "RasterStack") {
    examplerast = stack.static[[1]]
  }
  if (class(stack.static) == "RasterLayer") {
    examplerast = stack.static
  }
  locs = ctmc$ec
  wait.times = ctmc$rt
  notzero.idx = 1:ncell(examplerast)
  if (length(zero.idx) > 0) {
    notzero.idx = notzero.idx[-zero.idx]
  }
  adj = adjacent(examplerast, locs, pairs = TRUE, sorted = TRUE, 
                 id = TRUE, directions = directions, target = notzero.idx)
  adj.cells = adj[, 3]
  rr = rle(adj[, 1])
  time.idx = rep(rr$values, times = rr$lengths)
  start.cells = adj[, 2]
  z = rep(0, length(start.cells))
  idx.move = rep(0, length(z))
  diag.move = rep(0, length(locs))
  for (i in 1:(length(locs))) {
    idx.t = which(time.idx == i)
    idx.m = which(adj.cells[idx.t] == locs[i + 1])
    z[idx.t[idx.m]] <- 1
    if (length(idx.m) == 0) {
      diag.move[i] = 1
    }
  }
  tau = rep(wait.times, times = rr$lengths)
  t = rep(ctmc$trans.times, times = rr$lengths)
  if (nlayers(stack.static) > 1) {
    #X.static = values(stack.static)[start.cells, ]
    X.static = values(stack.static)[adj.cells,]
  }
  else {
    #X.static = matrix(values(stack.static)[start.cells], 
    #                  ncol = 1)
    X.static = matrix(values(stack.static)[adj.cells],ncol=1)
  }
  colnames(X.static) <- names(stack.static)
  xy.cell = xyFromCell(examplerast, start.cells)
  xy.adj = xyFromCell(examplerast, adj.cells)
  v.adj = (xy.adj - xy.cell)/sqrt(apply((xy.cell - xy.adj)^2, 
                                        1, sum))
  if (p.grad > 0) {
    X.grad = v.adj[, 1] * stack.gradient$grad.x[start.cells, 
    ] + v.adj[, 2] * stack.gradient$grad.y[start.cells, 
    ]
    if (grad.point.decreasing == TRUE) {
      X.grad = -X.grad
    }
    colnames(X.grad) <- colnames(stack.gradient$grad.x)
  }
  idx.move = which(z == 1)
  idx.move = c(idx.move, length(z))
  v.moves = v.adj[rep(idx.move[1:(length(rr$lengths) - 1)], 
                      times = rr$lengths[-1]), ]
  v.moves = rbind(matrix(0, ncol = 2, nrow = rr$lengths[1]), 
                  v.moves)
  X.crw = apply(v.moves * v.adj, 1, sum)
  if (crw == FALSE & p.grad > 0) {
    X = cbind(X.static, X.grad)
  }
  if (crw == TRUE & p.grad > 0) {
    X = cbind(X.static, X.grad, X.crw)
    colnames(X)[ncol(X)] = "crw"
  }
  if (crw == FALSE & p.grad == 0) {
    X = cbind(X.static)
  }
  if (crw == TRUE & p.grad == 0) {
    X = cbind(X.static, X.crw)
    colnames(X)[ncol(X)] = "crw"
  }
  if (include.cell.locations) {
    xys = cbind(xy.cell, xy.adj)
    colnames(xys) = c("x.current", "y.current", 
                      "x.adj", "y.adj")
    X = cbind(X, xys)
  }
  T = length(wait.times)
  p = ncol(X)
  out = data.frame(z = z, X, tau = tau, t = t)
  T = nrow(out)
  out = out[-((T - 3):T), ]
  out
}
debug1 <- function(){
  rm(list=ls())
  load(file="simu_4.rda")
  load(file="ctmc2glm_study_debug2.rda")
  source("ctmc2glm_study_3.R")
  debug(ctmc2glm.study)
  with(tmp,ctmc2glm.study(
    xy=cbind(x,y),t=t,tmap=p,cov__,NULL))
}
test <- function()
{
  rm(list=ls())
  library(ctmcmove)
  #setwd("C:/Users/dliang/Google Drive/Summer_2016/Aimee")
  load("sim_test_1.rda")
  args(ctmc2glm.study)
  #function (xy, t, tmap, stack.static, stack.grad, ...) 
  dim(cov__[[1]])
  #[1] 25 25 60
  range(dat$t)
  #[1]   1 180
  
  ## argument
  xy <- as.matrix(dat[,1:2])
  t <- as.vector(dat[,3])
  tmap <- dat[,4]  ## three records per image
  
  source("ctmc2glm_study_2.R")
  ## static only
  glmdata0 <- ctmc2glm.study(xy,t,tmap,list(raster(cov__[[1]],layer=1)),NULL)
  str(glmdata0)
  # 'data.frame':	334 obs. of  9 variables:
  # $ z        : num  0 1 0 0 0 0 1 0 0 0 ...
  # $ layer.1  : num  0.929 0.929 0.929 0.929 -0.722 ...
  # $ crw      : num  0 0 0 0 -1 1 0 0 0 0 ...
  # $ x.current: num  12 12 12 12 13 ...
  # $ y.current: num  12 12 12 12 12 ...
  # $ x.adj    : num  11.1 13 12 12 12 ...
  # $ y.adj    : num  12 12 13 11.1 12 ...
  # $ tau      : num  0.657 0.657 0.657 0.657 0.186 ...
  # $ t        : num  1.66 1.66 1.66 1.66 1.84 ...
  extract(raster(cov__[[1]],layer=1),glmdata0[1:10,6:7])- glmdata0[1:10,2]
  # [1]  0  0  0  0  0 NA NA  0  0 NA
  
  glmdata1 <- ctmc2glm.study(xy,t,tmap,NULL,list(raster(cov__[[1]],layer=1)))
  # Error in ctmc2glm.st(xy, t, tmap, NULL, list(raster(cov1, layer = 1))) : 
  #   statck.static must not be NULL
  glmdata2 <- ctmc2glm.study(xy,t,tmap,list(raster(cov__[[1]],layer=1)),
                          list(raster(cov__[[1]],layer=1)))
  str(glmdata2)
  # 'data.frame':	334 obs. of  10 variables:
  # $ z        : num  0 1 0 0 0 0 1 0 0 0 ...
  # $ layer.1  : num  0.929 0.929 0.929 0.929 -0.722 ...
  # $ X.grad   : num  -7.75e-06 7.75e-06 -8.58e-07 8.58e-07 7.01e-07 ...
  # $ crw      : num  0 0 0 0 -1 1 0 0 0 0 ...
  # $ x.current: num  12 12 12 12 13 ...
  # $ y.current: num  12 12 12 12 12 ...
  # $ x.adj    : num  11.1 13 12 12 12 ...
  # $ y.adj    : num  12 12 13 11.1 12 ...
  # $ tau      : num  0.657 0.657 0.657 0.657 0.186 ...
  # $ t        : num  1.66 1.66 1.66 1.66 1.84 ...
  
  X__ <- list(x1=stack(cov__[[1]]),x2=stack(cov__[[2]]))
  ## dynamic for location based covariate
  ## no gradient
  debug(ctmc2glm.study)
  glmdata3 <- ctmc2glm.study(xy,t,tmap,X__,NULL,verbose=999)
  str(glmdata3)
  # 'data.frame':	334 obs. of  10 variables:
  # $ z        : num  0 1 0 0 0 0 1 0 0 0 ...
  # $ layer.1  : num  0.929 0.929 0.929 0.929 -0.722 ...
  # $ crw      : num  0 0 0 0 -1 1 0 0 0 0 ...
  # $ x.current: num  12 12 12 12 13 ...
  # $ y.current: num  12 12 12 12 12 ...
  # $ x.adj    : num  11.1 13 12 12 12 ...
  # $ y.adj    : num  12 12 13 11.1 12 ...
  # $ tau      : num  0.657 0.657 0.657 0.657 0.186 ...
  # $ t        : num  1.66 1.66 1.66 1.66 1.84 ...
  # $ tmap__   : int  1 1 1 1 1 1 1 1 1 1 ...
  unique(glmdata3$tmap__)
  line <- subset(glmdata3, tmap__==2)
  eps <- extract(raster(cov__[[1]],layer=2),line[,7:8])- line[,2]
  summary(eps)
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  #    0       0       0       0       0       0       2 
  
  ## static gradient
  #debug(ctmc2glm.st)
  glmdata4 <- ctmc2glm.study(xy,t,tmap,X__,list(raster(cov__[[1]],layer=1)),verbose=99)
  ## Error in gradient not supported
  ## ignore below for now.
  
  ## dynamic gradient
  glmdata5 <- ctmc2glm.st(xy,t,tmap,list(cov1),list(cov1),verbose=99)
  
  ## static stack
  summary(glmdata5[,2]-glmdata4[,2])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0       0       0       0       0       0 
  summary(glmdata5[,2]-glmdata3[,2])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0       0       0       0       0       0 
  
  ## gradient layer
  summary(glmdata4[,3]-glmdata2[,3])
}
#require(raster)
ctmc2glm.st.stack <- function(lst,i)
{
  ## helper function to build stack from a list of raster stack:
  ## lst: list of raster stack
  ## i : layer within each stack to extract
  ## value: a raster stack
  out <- vector("list",length(lst))
  for(idx_ in 1:length(lst)){
    out[[idx_]] <- raster(lst[[idx_]],layer=i)
  }
  if(length(out)>1){
    r <- stack(out)
  }
  else{
    r <- out[[1]]
  }
  r
}
ctmc2glm.st.crw.impute <- function(x)
{
  ## input missing values of crw with last known value
  ## if not found, with the most adjacent no missing value
  if(all(is.na(x))){
    warning("nothing to impute\n")
    return(x)
  }
  if(any(is.na(x))){
    na.idx <- which(is.na(x))
    nna.idx <- which(!is.na(x))
    for( i in 1:length(na.idx)){
      if(na.idx[i]<min(nna.idx)){
        ## no previous non-missing
        x[na.idx[i]] <- x[min(nna.idx)]
      }
      else{
        ## with previous non-missing
        sel <- max(nna.idx[nna.idx<na.idx[i]]) ## first most recent non-missing id
        x[na.idx[i]] <- x[sel]
      }
    }
  }
  x
}


#
# helper function to build static and gradient stack given a list arguments
# build a rasterlayer or stack based on the list input
#
# arugments: n2: the class name
# stack.grad: the input list of raster objects
# layerid: the layer to extract
# return the raster object for glm.data modeling;
ctmc2glm.st.build <- function(n2,stack.grad,layerid=1)
{
  if(n2=="null"){
    l.stack.grad <- NULL
    nvar <- 0
    examplerast <- NULL
  }
  else if(n2=="RasterLayer"){ ## static covariates
    nvar <- length(stack.grad)
    if(length(stack.grad)>1){
      l.stack.grad <- stack(stack.grad)
      examplerast <- raster(l.stack.grad,layer=layerid)
    }
    else{
      l.stack.grad <- stack.grad[[1]]
      examplerast <- l.stack.grad
    }
  }                           ## end static covariates
  else if(n2=="RasterStack"){ ## dynamic covariates
    nvar <- length(stack.grad)
    l.stack.grad <- ctmc2glm.st.stack(stack.grad,layerid)
    if(class(l.stack.grad)=="RasterLayer"){
      examplerast <- l.stack.grad
    }
    else{
      examplerast <- raster(l.stack.grad,layer=1)
    }
  }                           ## end dynamic covariates
  list(stack=l.stack.grad,example=examplerast,nvar=nvar)
}

#
# helper function to return the nearest time in data to query
require(FNN)
ctmc2glm.knn <- function(dataTime,queryTime,tmax=NULL)
{
  if(is.null(tmax)){
    tmax <- diff(range(dataTime))
  }
  knn.temp <- get.knnx(data=matrix(dataTime,ncol=1),
                       query=matrix(queryTime,ncol=1),
                       k=1)
  if(any(knn.temp$nn.dist>tmax)){
    stop("No time found in data")
  }
  knn.temp$nn.index
}

test1 <- function()
{
  x <- c(1,NA,2,NA,3)
  ctmc2glm.st.crw.impute(x)
  #[1] 1 1 2 3 3
  ctmc2glm.st.crw.impute(1:4)
  #[1] 1 2 3 4
  x <- c(NA,NA,2,NA,3)
  #debug(ctmc2glm.st.crw.impute)
  ctmc2glm.st.crw.impute(x)
  #[1]  2  2  2 NA  3
  ctmc2glm.st.crw.impute(rep(NA,4))
}
