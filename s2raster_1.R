#object: basis raster
#x: plotted smooth from mgcv
s2raster <- function(object,x,se=FALSE,...){
  xy <- expand.grid(x=x[[1]]$x,y=x[[1]]$y)
  if(se){
    xyz <- cbind(xy,x[[1]]$se)
  }else{
    xyz <- cbind(xy,x[[1]]$fit)
  }
  xyz <- na.omit(xyz)
  intensity(object,xyz,...)
}
dev_ <- function(){
  rm(list=ls())
  library(mgcv)
  library(raster)
  library(MBA)
  library(automap)
  source("intensity_2.R")
  load(file="Scratch/s2raster_dev.rda")
  load(file="../leatherback_4_Dong_new_noIATTC_pa_3b.rda")
  rx <- range(pa3$point_x,na.rm=T)
  ry <- range(pa3$point_y,na.rm=T)
  x[[1]]$x <- x[[1]]$x*abs(diff(rx)) + rx[1]
  x[[1]]$y <- x[[1]]$y*abs(diff(ry)) + ry[1]
  object <- raster("../../Data/sp_raster_proj.nc")
  source("s2raster_1.R")
  out <- s2raster(object,x)
  plot(out,xlim=rx,ylim=ry,asp=F)
}