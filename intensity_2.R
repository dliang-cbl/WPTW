## intensity_2.R
## ver 2) allow more resolution to be shown
require(raster)
require(MBA)
intensity.rescale <- function(object,eps=NULL)
{
  if(!is.null(eps)){
    kfac <- max(ceiling(eps/res(object)))
    eraster1 <- object
    res(eraster1) <- eps
    if(kfac>1){
      object <- aggregate(object,fact=kfac)
    }
    output <- resample(object,eraster1,method="ngb")
  }
  else{
    output <- object
  }
  output
}
intensity <- function(object,xyz,eps=NULL,n=1,m=1,h=8,extend=FALSE,sp=FALSE,...)
{
  ## helper script to show irregular data on a raster by MBA package.
  ## object: a raster object
  ## xyz   : xyz matrix to visualize on raster
  ## eps   : resolution of the raster object
  ## Value: a raster object for plotting

  object <- intensity.rescale(object,eps=eps)  ## rescale the raster* object
  
  l.xyz <- rasterToPoints(object)      ## extract all non-missing cells
  l.id <- cellFromXY(object,l.xyz[,1:2]) ## the cell number in the raster 
  
  l.xyz.est <- mba.points(as.matrix(xyz),l.xyz[,1:2])[[1]] ## multilevel B-splines interpolation on raster
  
  r.object <- object
  r.object[l.id] <- l.xyz.est[,3]
  
  r.object
}

test2 <- function()
{
  rm(list=ls())
  .libPaths(c(.libPaths(),"~/Rpacks"))
  library(raster)
  ## Load the extended raster
  load("./sp_raster.rda")
  ## Load model results
  load("Scratch/inla3c_data.rda")
  load("Scratch//inla3c.rda")
  
  source("intensity_1.R")
  l.fit <- subset(cbind(glm0b,test1c$summary.fitted.value),z==0)
  l.obs <- subset(cbind(glm0b,test1c$summary.fitted.value),z>0)
  
  zlim <- range(l.fit[,10])
  
  ## June prediction
  l.xyz <- subset(l.fit,z==0 & amonth==6,select=c(2,3,10))
  l.xyz.obs <- subset(l.obs, amonth==6, select=c(2,3,10))
  source("choropleth_raster_1.R")
  fit.june <- intensity(sp,l.xyz)
  pal <- tim.colors(8)
  brks <- c(0,0.001,0.025,0.05,0.1,0.5,1,2,5.359)
  png("Scratch/fit_legend.png",height=1000,width=1000,res=400)
  fit.june.chor <- choropleth.raster(
    fit.june,pal,"Intensity",legend=TRUE,
    style="fixed",fixedBreaks=brks,dataPrecision=4)
  dev.off()
  
  png("Scratch/fit_june_1.png",height=4000,width=6000,res=500)
  op <- par(mar=c(3.1,3.1,1.1,1.1))
  #plot(fit.june,zlim=zlim,col=tim.colors())
  plot(fit.june.chor,col=pal,legend=FALSE)
  points(l.xyz.obs[,1],l.xyz.obs[,2])
  par(op)
  dev.off()
  
  ## December prediction
  l.xyz <- subset(l.fit,z==0 & amonth==12,select=c(2,3,10))
  l.xyz.obs <- subset(l.obs, amonth==12, select=c(2,3,10))
  fit.dec <- intensity(sp,l.xyz)
  fit.dec.chor <- choropleth.raster(
    fit.dec,pal,"intensity",style="fixed",fixedBreaks=brks,dataPrecision=4)
  png("Scratch/fit_december_1.png",height=4000,width=6000,res=500)
  op <- par(mar=c(3.1,3.1,1.1,1.1))
  #plot(fit.dec,zlim=zlim,col=tim.colors())
  plot(fit.dec.chor,col=pal,legend=FALSE)
  points(l.xyz.obs[,1],l.xyz.obs[,2],col="green")
  par(op)
  dev.off()
  
  for(i in 1:12){
    l.xyz <- subset(l.fit,z==0 & amonth==i,select=c(2,3,10))
    l.xyz.obs <- subset(l.obs, amonth==i, select=c(2,3,10))
    l.month.fit <- intensity(sp,l.xyz)
    l.month.fit.chor <- choropleth.raster(
      l.month.fit,pal,"intensity",
      style="fixed",fixedBreaks=brks,dataPrecision=4
    )
    png(paste("Scratch/fit_amonth_",i,"_1.png",sep=""),height=4000,width=6000,res=500)
    op <- par(mar=c(3.1,3.1,1.1,1.1))
    #plot(l.month.fit,zlim=zlim)
    plot(l.month.fit.chor,col=pal,legend=FALSE)
    points(l.xyz.obs[,1],l.xyz.obs[,2],col="black")
    par(op)
    dev.off()
  }
  
}

test <- function()
{
  rm(list=ls())
  library(raster)
  r <- raster(nrows=10, ncols=10)
  xy <- xyFromCell(r,1:100)
  mu <- apply((xy/200)^2,1,sum)
  values(r) <- rnorm(100,mu,0.1)
  set.seed(112300)
  r[sample(100,40)] <- NA
  plot(r)
  
  r1 <- r
  res(r1) <- res(r)/50
  dim(r1)
  obs <- cbind(runif(50,-250,250),runif(50,-80,80),rnorm(50))
  
  xyz <- obs
  object <- r1
  
  library(fields)
  quilt.plot(xyz[,1],xyz[,2],xyz[,3],nrow=10,ncol=10)
  
  l.xyz <- rasterToPoints(object)      ## extract all non-missing cells
  l.id <- cellFromXY(object,l.xyz[,1:2]) ## the cell number in the raster 
  
  require(MBA)
  l.xyz.est <- mba.points(xyz,l.xyz[,1:2])[[1]] ## multilevel B-splines interpolation on raster
  
  r.object <- object
  r.object[l.id] <- l.xyz.est[,3]
  plot(r.object)
}