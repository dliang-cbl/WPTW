## ppm_st_glm_6.R
## convert the spatio-temporal point pattern and covariates into a poisson regression data
## with weights

## ver 6) allow predictive analysis from future raster* object
## sub ver d) make total weight equaling the raster area
## sub ver e) fix presence environmental data
## sub ver f) modulize code to generate data for prediction
##            and build stack from data

## BUG:  the environmental condition of the data should not change
##       only the resolution of the model change
## BUG:  for varying resolution, the total weights do not equal the raster area
## BUG:  can not change resolution of raster directly without maintaining the masks
## BUG:  can not subset time period to observed.
## BUG:  when data are scaled up, some near shore observations were scaled outside
## the study area, leading to problem with updating the weights of the pseudo-absence
## points, these points were now removed altogether.


## arguments
## x : a R data.frame containing the geographic coorindates and time stamp of the point process
## window: a R raster* object with NA denoting masked area and 1 the study area
## covariates: a R list of R::raster stack object defining environmental variables
##  ( Note: these stack objects must be re-sampled to the same resolution)
##  ( And: the same resolution was used to discretize the point process by default)
## trend: a R formula object defining the covariate name
## coordinate: a R formula object defining the x/y and t coordinates in data.frame x
## eps : Optional. A positive number, or a vector of two positive numbers, giving
##       the horizontal and vertical spacing respectively, of the grid of 
##       dummy points. (default being the environmental raster resolution)
## newdata   : a R list of R::raster stack object defining environmental variables
##             for prediction, they should be the same format as covariates
##             and have the same named arguments
## newdata.time: a R vector denoting the time-stamp for each layer of newdata
##              default to -1
## nfac: factor to search for nearest neighbor in the raster for observation point
## values: a R data.frame with variables
## z: auxiliary variable (1/w) if observed, zero otherwise
## w: weights equal area each point covers
## covariates listed as input
## x: x coordinate
## y: y coordinate
## t: t coordinate

ppm.st.glm <- function(
  x,window,covariates,trend=~1,coordinate=~x+y+t,eps=NULL,
  newdata=NULL,newdata.ctime=NULL,nfac=2,verbose=TRUE)
{
  x <- as.data.frame(x)
  xcoord <- all.vars(coordinate)[1]
  ycoord <- all.vars(coordinate)[2]
  tcoord <- all.vars(coordinate)[3]
  
  ## determine the quadrature grid of the study domain
  eraster <- window
  
  ## develop quadrature  grid
  if(!is.null(eps)){
    kfac <- max(ceiling(eps/res(eraster)))
    eraster1 <- eraster
    res(eraster1) <- eps
    if(kfac>1){
      ## if grid resolution larger than input resolution,
      ## then aggregate the input data and resample to the grid using nearest neighbor
      eraster0 <- aggregate(eraster,fact=kfac)
      eraster <- resample(eraster0,eraster1,method="ngb")
    }
    else{
      ## if grid resolution smaller than input resolution,
      ## then disaggregate the input data using nearest neighbor algorithm
      eraster <- resample(eraster,eraster1,method="ngb")
    }
  }
  
  #browser()  
  area <- prod(res(eraster))
  if(!is.null(eps)){
    ## adjust the area so the total area matches the input window
    orig.area <- prod(res(window))*prod(dim(window))
    new.area <- area*prod(dim(eraster))
    area.fac <- orig.area / new.area 
    area <- area * area.fac
  }
  
  ## quadrature points locations
  xyz <- rasterToPoints(eraster)
  id <- cellFromXY(eraster,xyz[,1:2])
  colnames(xyz)[1:2] <- c(xcoord,ycoord)
  
  ## build psedo-absence points by time periods 
  ## (assuming time periods are equal)
  ## pseduo-absence points from window
  nt.range <- range(x[,tcoord],na.rm=T)
  nt.data <- seq(nt.range[1],nt.range[2])
  nt <- length(nt.data)
  
  grid <- vector("list",nt)
  #browser()
  for(ii in 1:nt){
    i <- nt.data[ii]
    if(verbose) cat("pseduo-absence at time= ",i,"\n") ## for each time period,
    
    ## x/y coordinates
    ## define time coordinates 
    ## assume time stamp in x matches with the layer id in covariates
    ## therefore could use the time-stamp in x directly to query the layers in grid
    ##
    ## weight as area per cell per unit time
    
    ## extract environmental data from quadrature points
    l.stack <- extract.helper(object = covariates,layer = i)
    l.tmp <- extract(l.stack,xyz[,1:2])
    
    l.cell <- id
    if(length(all.vars(trend))==0){
      ## no trend
      grid[[i]] <- cbind(cell=l.cell,xyz[,1:2],t=i,wt=area)
    }
    else{
      ## with trend
      grid[[i]] <- cbind(cell=l.cell,xyz[,1:2],l.tmp,t=i,wt=area)
    }
  }
  
  ## assign the observed locations to the discretization grid
  obs.lst <- split(x,x[,tcoord])
  res__ <- max(mean(res(window)),mean(res(eraster))) ## resolution to search for nearest neighbor
  for(j in 1:length(obs.lst)){
    if(verbose) cat("j= ",j,"\n") ## for each time period,
    
    ## the following section identifies the observed location in the grid
    ## count the number of grids with observed data,
    ## define the weight for each observed location as the raster area -
    ## divided by the number of observed location plus one (ie the grid point)
    
    cells <- cellFromXY.nbr(eraster,obs.lst[[j]][,c(xcoord,ycoord)],
                            eps=nfac*res__)
    stopifnot(all(!is.na(cells))) ## stop here if data not in study area
    
    tbl <- as.data.frame(table(cells))
    cells.unique <- as.numeric(as.character(tbl$cells))
    cells.freq <- tbl$Freq
    
    idx.t <- obs.lst[[j]][1,tcoord]
    l.cell.na <- match(cells,grid[[idx.t]][,"cell"])
    obs.lst[[j]]$cell <- cells
    obs.lst[[j]]$freq_ <- cells.freq[match(cells,cells.unique)] + 1
    obs.lst[[j]]$wt <- area / obs.lst[[j]]$freq_
    
    ## extract the covariates from the corresponding layer
    l.stack <- extract.helper(object = covariates,layer = idx.t)
    tmp <- extract(l.stack,obs.lst[[j]][,c(xcoord,ycoord)])
    
    
    if(length(all.vars(trend))>0){
      obs.lst[[j]] <- cbind(obs.lst[[j]],tmp)
    }
    
    ## remove out-of-range observations
    obs.lst[[j]] <- subset(obs.lst[[j]],!is.na(l.cell.na))
    
    ##  adjusts the weights in the discretization grid
    ## so that the total weights equals the area
    l.grid <- match(cells.unique,grid[[idx.t]][,"cell"])
    l.grid.nna <- which(!is.na(l.grid))
    grid[[idx.t]][l.grid[l.grid.nna],"wt"] <- area/(cells.freq[l.grid.nna]+1)
    
  }
  
  ## combine the discretization grid and the observed data
  tmp1 <- do.call(rbind,grid)
  colnames(tmp1)[ncol(tmp1)-1] <- tcoord
  varnames <- colnames(tmp1)
  tmp2 <- do.call(rbind,obs.lst)[,varnames]
  #browser()
  
  value <- rbind(data.frame(tmp1,z=0),data.frame(tmp2,z=1/tmp2[,"wt"]))
  
  #browser()
  ## build prediction data - follow the same routine as for pseudo-absence data
  if(!is.null(newdata)){
    np <- dim(newdata[[1]])[3]
    l.pred <- vector("list",np)
    if(is.null(newdata.ctime)){
      newdata.ctime <- rep(-1,np)
    }
    for(k in 1:np){
      if(verbose) cat("prediction at time= ",k,"\n") ## for each time period
      
      ## extract environmental predictions
      l.stack <- extract.helper(object=newdata,layer=k)
      l.tmp.2 <- extract(l.stack,xyz[,1:2])
      
      l.cell <- id
      if(length(all.vars(trend))==0){
        ## no trend
        l.pred[[k]] <- cbind(cell=l.cell,xyz[,1:2],t=newdata.ctime[k],wt=area)
      }
      else{
        ## with trend
        l.pred[[k]] <- cbind(cell=l.cell,xyz[,1:2],l.tmp.2,t=newdata.ctime[k],wt=area)
      }
    }
    tmp3 <- do.call(rbind,l.pred)
    colnames(tmp3)[ncol(tmp3)-1] <- tcoord
  }
  
  
  if(!is.null(newdata)){
    value <- rbind(value,data.frame(tmp3,z=NA))
  }
  
  value
}
ppm.st.glm.window <- function(window,eps=NULL){
  ## determine the quadrature grid of the study domain
  eraster <- window
  
  ## develop quadrature  grid
  if(!is.null(eps)){
    kfac <- max(ceiling(eps/res(eraster)))
    eraster1 <- eraster
    res(eraster1) <- eps
    if(kfac>1){
      ## if grid resolution larger than input resolution,
      ## then aggregate the input data and resample to the grid using nearest neighbor
      eraster0 <- aggregate(eraster,fact=kfac)
      eraster <- resample(eraster0,eraster1,method="ngb")
    }
    else{
      ## if grid resolution smaller than input resolution,
      ## then disaggregate the input data using nearest neighbor algorithm
      eraster <- resample(eraster,eraster1,method="ngb")
    }
  }
  
  area <- prod(res(eraster))
  if(!is.null(eps)){
    ## adjust the area so the total area matches the input window
    orig.area <- prod(res(window))*prod(dim(window))
    new.area <- area*prod(dim(eraster))
    area.fac <- orig.area / new.area 
    area <- area * area.fac
  }
  
  ## quadrature points locations
  xyz <- rasterToPoints(eraster)
  id <- cellFromXY(eraster,xyz[,1:2])
  list(value=eraster,id=id,xyz=xyz,area=area)
}
ppm.st.glm.fill <- function(
  window,newdata,formula,eps=NULL,verbose=FALSE){
  ##window: window to store data frame objects
  ##newdata: prediction on the same grid
  ##formula: LHS: value to fill
  ##         RHS: x+y+t coordinates in new data
  ##eps    : resolution to fill the prediction
  ##verbose: show intermediate results
  ## value: a raster Stack of resolution eps
  var_ <- all.vars(formula)
  stopifnot(all(var_ %in% names(newdata)))
  
  fill.lst <- split(newdata,newdata[,var_[4]])
  fill.times <- names(fill.lst)
  l.xyz2 <- vector("list",length(fill.lst))
  names(l.xyz2) <- fill.times
  if(!is.null(eps)){
    if(verbose){
      cat("downscale resolution=",eps,".\n")
    }
    for(i in 1:length(fill.lst)){
      if(verbose){
        cat("processing ",fill.times[i],"\n")
      }
      l.xyz2[[i]] <- intensity(
        object=window,
        xyz = fill.lst[[i]][,var_[c(2,3,1)]],
        eps = eps
      )
    }
  }else{
    for(i in 1:length(fill.lst)){
      if(verbose){
        cat("processing ",fill.times[i],"\n")
      }
      xyz <- fill.lst[[i]][,var_[c(2,3,1)]]
      cell <- cellFromXY(window,xyz[,-3])
      line_ <- window
      line_[] <- NA
      line_[cell] <- xyz[,3]
      l.xyz2[[i]] <- line_
    }
  }
  stack(l.xyz2)
}
ppm.st.glm.new <- function(window,newdata,coord=~x+y+t,
                           newdata.ctime=NULL,eps=NULL,
                           verbose=FALSE){
  ## window to data frame objects
  window_ <- ppm.st.glm.window(window,eps)
  
  colnames(window_$xyz)[1:2] <- all.vars(coord)[1:2]
  ## build prediction data - follow the same routine as for pseudo-absence data
  np <- dim(newdata[[1]])[3]
  l.pred <- vector("list",np)
  if(is.null(newdata.ctime)){
    if(!is.null(names(newdata[[1]]))){
      newdata.ctime <- names(newdata[[1]])  
    }else{
      newdata.ctime <- seq(1,dim(newdata[[1]])[3])
    }
  }
  for(k in 1:np){
    if(verbose) cat("prediction at time= ",k,"\n") ## for each time period
    ## extract environmental predictions
    l.stack <- extract.helper(object=newdata,layer=k)
    l.tmp.2 <- extract(l.stack,window_$xyz[,1:2])
    l.pred[[k]] <- data.frame(
      cell=window_$id,
      window_$xyz[,1:2],
      l.tmp.2,
      wt=window_$area,
      t=newdata.ctime[k])
  }
  pred <- do.call(rbind,l.pred)
  colnames(pred)[ncol(pred)] <- all.vars(coord)[3]
  pred
}
test <- function()
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  library(FNN)
  library(raster)
  library(mgcv)
  source("cellFromXY_nbr.R")
  source("ppm_st_utils.R")
  source("extract_1.R")
  load(file="simu_6.rda")
  source("ppm_st_glm_6f.R")
  base <- raster(cov__[[1]],layer=1)
  eraster <- ppm.st.glm.window(base,eps=1)
  newd <- ppm.st.glm.new(base,
                         newdata = cov__,
                         newdata.ctime = names(cov__$bathy),
                         coord = ~long+lat+amonth,
                         eps = 0.5)
  ## select a starting model to evaluate the resolutions
  loc <- newd[sample(10000,100),]
  loc$t <- match(loc$amonth,names(cov__[[1]]))
  quad <- ppm.st.glm(loc,window=base,
                     eps=0.5,
                     covariates = cov__, 
                     trend=~bathy+sst+sst2+fpi+ssh+chl,
                     coordinate = ~long+lat+t,
                     verbose = FALSE)
}
