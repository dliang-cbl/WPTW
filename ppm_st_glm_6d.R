## ppm_st_glm_6.R
## convert the spatio-temporal point pattern and covariates into a poisson regression data
## with weights

## ver 6) allow predictive analysis from future raster* object
## sub ver d) make total weight equaling the raster area

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
require(FNN)
#require(INLA)
#inla.setOption(enable.inla.argument.weights=TRUE)

## minor update 5/30) change simulation to include raster area.

ppm.st.glm <- function(x,window,covariates,trend=~1,coordinate=~x+y+t,eps=NULL,
                       newdata=NULL,newdata.ctime=NULL,nfac=2,verbose=TRUE)
{
  x <- as.data.frame(x)
  xcoord <- all.vars(coordinate)[1]
  ycoord <- all.vars(coordinate)[2]
  tcoord <- all.vars(coordinate)[3]
  
  ## determine the discretization grid of the study domain
  eraster <- window
  
  ## resample discretization grid
  if(!is.null(eps)){
    kfac <- max(ceiling(eps/res(eraster)))
    #res(eraster) <- eps
    eraster1 <- eraster
    res(eraster1) <- eps
    if(kfac>1){
      ## if grid resolution larger than input resolution,
      ## then aggregate the input data and resample to the grid using nearest neighbor
      eraster0 <- aggregate(eraster,fact=kfac)
      eraster <- resample(eraster0,eraster1,method="ngb")
      
      covariates <- lapply(covariates,function(elmt){
        r2 <- aggregate(elmt,fact=kfac)
        resample(r2,eraster,method="ngb")
      })
     
      ## also process prediction data
      if(!is.null(newdata)){
        newdata <- lapply(newdata,function(elmt){
          r2 <- aggregate(elmt,fact=kfac)
          resample(r2,eraster,method="ngb")
        })
      } 
    }
    else{
      ## if grid resolution smaller than input resolution,
      ## then disaggregate the input data using nearest neighbor algorithm
      eraster <- resample(eraster,eraster1,method="ngb")
      
      covariates <- lapply(covariates,function(elmt){
        resample(elmt,eraster,method="ngb")
      })
      
      ## if new data, also process it
      if(!is.null(newdata)){
        newdata <- lapply(newdata,function(elmt){
          resample(newdata,eraster,method="ngb")
        })
      }
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
  
  ## build psedo-absence points by time periods 
  ## (assuming time periods are equal)
  ## pseduo-absence points from window
  xyz <- rasterToPoints(eraster)
  id <- cellFromXY(eraster,xyz[,1:2])
  colnames(xyz)[1:2] <- c(xcoord,ycoord)
  
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
    
    l.tmp <- sapply(covariates,function(elmt){
      elmt[[i]][id]
    })
    colnames(l.tmp) <- names(covariates)
    
    l.cell <- id
    
    if(length(all.vars(trend))==0){
      ## no trend
      grid[[i]] <- cbind(cell=l.cell,xyz[,1:2],t=i,wt=area)
    }
    else{
      ## with trend
      grid[[i]] <- cbind(cell=l.cell,xyz[,1:2],l.tmp,t=i,wt=area)
    }
    #colnames(grid[[i]])[2:4] <- c(xcoord,ycoord,tcoord)
  }

  ## assign the observed locations to the discretization grid
  obs.lst <- split(x,x[,tcoord])
  res__ <- max(mean(res(window)),mean(res(eraster))) ## resolution to search for nearest neighbor
  for(j in 1:length(obs.lst)){
    if(verbose) cat("j= ",j,"\n") ## for each time period,
    #if(j==7){
    #  browser()
    #}
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
    tmp <- sapply(covariates,function(elmt){
      elmt[[idx.t]][cells]
    })

    if(is.null(dim(tmp))){
      tmp <- matrix(tmp,1,length(tmp))
    }
    colnames(tmp) <- names(covariates)
    if(length(all.vars(trend))>0){
      obs.lst[[j]] <- cbind(obs.lst[[j]],tmp)
    }
    
    ## remove out-of-range observations
    obs.lst[[j]] <- subset(obs.lst[[j]],!is.na(l.cell.na))
    
    ##  adjusts the weights in the discretization grid
    ## so that the total weights equals the area
    l.grid <- match(cells.unique,grid[[idx.t]][,"cell"])
    l.grid.nna <- which(!is.na(l.grid))
    #grid[[idx.t]][l.grid,"wt"] <- area/(cells.freq+1)
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
      l.tmp.2 <- sapply(newdata,function(elmt){
        elmt[[k]][id]
      })
      colnames(l.tmp.2) <- names(newdata)
      
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

test <- function()
{
  rm(list=ls())
  sptw <- read.csv("leatherback_4.csv",head=T)
  
  library(raster)
  ## Load the extended raster
  load("./sp_raster.rda")
  load("./SST_aggregate_data.Rdata")
  load("./bathyMask_median_nn.Rdata")
  plot(bathyMask_median_nn)
  
  bathy1 <- vector('list',dim(allSST_median_nn)[3])
  for(i in 1:length(bathy1)){
    bathy1[[i]] <- bathyMask_median_nn/1000
  }
  bathy <- stack(bathy1)
  
  ## build pseudo-absence data
  tmp <- with(sptw,as.Date(paste(ayear,amonth,1,sep="/"),"%Y/%m/%d"))
  sptw$cmonth <- match(format(tmp,"%b%y"),names(allSST_median_nn))
  with(subset(sptw,is.na(cmonth)),table(ayear))
  #ayear
  #1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 
  #3    9   12    9   10    6    8    7    8    7   10   54  113 
  
  sptw1 <- subset(sptw,!is.na(cmonth))
  
  bathy.pred <- subset(bathy,c(162,168))
  sst.pred <- subset(allSST_median_nn,c(162,168))
  #save(sptw1,sp,bathy,allSST_median_nn,file="debug_1.rda")
  source("ppm_st_glm_6c.R")
  options(warn=1)
  glm1 <- ppm.st.glm(
    x=sptw1,window=sp,covariates=list(bathy=bathy,sst=allSST_median_nn),
    trend=~bathy+sst,coordinate = ~long+lat+cmonth,
    newdata=list(bathy=bathy.pred,sst=sst.pred),newdata.ctime=c(162,168),
    eps=3.0)
}
