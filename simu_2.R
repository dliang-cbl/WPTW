rm(list=ls())
options(warn = 1)
library(raster); library(ctmcmove); library(mgcv);
options("rgdal_show_exportToProj4_warnings"="none")

list.files(pattern="tif")

## define seascape variables as a list
seascape <- list(bathy=stack("bathy.tif"),
              sst=stack("sst.tif"),
              chl=stack("chl.tif"),
              ssh=stack("ssh.tif"))


## load tracking data codes
source("nmiss_2.R")
source("s2inla_6.R")
source("ctmc2glm_study_3.R")

# load simulated tracks
tracks0 <- read.csv("tracks.csv")
str(tracks0)
## aday: daily time stamp
## amon: monthly time stamp (continuous month to link
##                           with layers of seascape)
tracks <- split(tracks0,~id)

# Process tracks into a data structure
# that allows generalized linear modeling of telemetry data

tmplist <- vector("list",length(tracks))
for(i in 1:length(tracks)){
  cat("track",i,"/",length(tracks),".\n")
  line <-  with(tracks[[i]],ctmc2glm.study(
    xy=cbind(long,lat),t=aday,tmap=amon,
    seascape,NULL))
  if(!is.null(line)){
    line$id <- i  
  }
  tmplist[[i]] <- line
}

## estimating population level "resource selection"
track.data <- do.call(rbind,tmplist)
track.data$id <- factor(track.data$id)

## rename the environmental covariates
names(track.data)[1+1:length(seascape)] <- names(seascape)



cat("processing fishery observations.\n")
loc <- read.csv("loc.csv")
str(loc)
## amon: continuous month to match seascape layers

## load point process codes
source("cellFromXY_nbr.R")
source("ppm_st_glm_6d.R")
source("ppm_st_utils.R")


## study area raster layer (NA=outside area, non-NA=in the area)
sp <- raster(seascape[["bathy"]],layer=1)

# generate pseudo-absence data 
args(ppm.st.glm)

pa <- ppm.st.glm(
  loc, 
  window = sp, #a R raster* object with NA denoting masked area and 1 the study area,
  covariates = seascape,
  trend=reformulate(names(seascape)),
  coordinate=~long+lat+amon)


## Load ISDM modeling codes
library(abind)
source("all_group_1.R")
source("predict_gam_term_2.R")
source("sdm_bc_5d.R")
source("sim2jam_study_1.R")
source("term_utils.R")


## load prediction codes
library(doParallel)
source("predict_sdm_bc_5d.R")
source("ppm_st_glm_6f.R")
source("extract_1.R")
source("intensity_2.R")

## telemetry model: linking motility to seascape variables
term0 <- names(seascape)
for0 <- reformulate(term0,response="z")
## bias correction model
#fpr1 <- ~bathy+fpi+ssh
for1 <- ~s(long,lat,bs="ds",m=c(1,0.5),k=10)

## build scaling model list all subset for sdm terms
scale_terms <- all_group(
  reformulate(term0,response = "z"),verbose=F)
head(scale_terms)
## each row denotes a possible combination of terms
## value in the columns denote whether the scaling terms are equal
## for example 1,2,2,2 indicates that the last three terms
##  share the same scaling constant 

term1 <- scale_terms[1,]

ts0 <- proc.time()
bsdm.obj <- sdm.bc(
  ## see sdm_bc_5d.R for argument definitions
  formula = for0,
  term = term1,
  bias = for1,
  quad = pa,
  glmobj = track.data,
  wcoord = "wt", ## weight for pseudo-absence, 
  gcoord=NULL,#"id", ## terms defining individual
  tcoord = "tau",cterm="crw", ## terms defining residence time and last movement angle
  d=-1,S=300/8, ## lognormal prior mean and standard deviation for scalar
  timefact = 30,  ## 30 days in a month assumed
  niter =9,trace=T
)
ts1 <- proc.time()
cat("Joint analyses took",ts1-ts0,"\n")
save(bsdm.obj,file="bsdm_obj.rData")

## build predictive seascape
load(file="bsdm_obj.rData")
newdata <- list(bathy=stack("bathy.tif"),
                sst=stack("sst.tif"),
                chl=stack("chl.tif"),
                ssh=stack("ssh.tif"))

## build predictive data frame
newd <- ppm.st.glm.new(
  window=sp,newdata=newdata,
  coord=~long+lat+amon)

## amon is only a place holder for months

cat("predicting monthly intensity...\n")
(ts1 <- proc.time())
raw.pred <- predict_sdm_bc(bsdm.obj,newd,all=F,
                           verbose=F,batch=10)
(ts2 <- proc.time())


predMonth <- ppm.st.glm.fill(
  window = sp,
  newdata = raw.pred,
  formula = fit~long+lat+amon,
  verbose = T,
)

res(sp)
output <- raster(predMonth,layer=1)

library(raster)
library(sf)
library(tmap)
library(tmaptools)
tmap_mode("plot")
#library(OpenStreetMap)
library(RColorBrewer)
library(viridis)

summary(values(output))
col <-brewer.pal(n=5,"RdPu")
tm_shape(output)+tm_raster(
  palette=col,breaks = c(-999,0.043,0.044,0.045,0.046,0.048,999),
  labels = c("0 to 0.043","0.043 to 0.044","0.044 to 0.045",
             "0.045 to 0.046","0.046 to 0.048", "0.046=8 above"),
  title=expression(lambda))+
  tm_layout(legend.bg.color = "white")
#save.image(file=fn)

block2_build_data <- function(){
  rm(list = ls())
  library(terra)
  #library(raster)
  library(ctmcmove)
  
  args__ <- c("10","50")
  
  ntracks <- as.integer(args__[1])  ## number of tracks to simulate
  ndays <- as.integer(args__[2]) ## number of days to simulate
  
  ## output file name
  fn <- paste0("simu2/simuT_",ntracks,"L_",ndays,".rda")
  cat("output to ",fn,"\n")
  
  ## define seascape variables
  
  cov__ <- list(bathy=rast("simu2/bathy.tif"),
                sst=rast("simu2/sst.tif"),
                chl=rast("simu2/chl.tif"),
                ssh=rast("simu2/ssh.tif"))
  
  ## define resource selection coefficients
  load(file="simu2/coef.rda")
  
  
  ## load  intensity function (population level) codes
  source("residence_5.R")
  library(abind)
  
  ## load ctmc simulate codes
  source("ctmc_sim_study_5.R");
  source("ctmc_data_2.R")
  
  cat("simulating movement data each with ",
      ntracks," tracks of ", ndays,"days.\n")
  
  
  ## population level resource selection motility coefficients
  beta0 <- beta[-length(beta)]
  beta0[1] <- 2.273766
  crw0 <- crw
  
  ## simulate tracks
  tracks <- vector("list",ntracks)
  for(i in 1:ntracks){
    cat("simulating track ",i,'.\n')
    ## change 7/28/2022 DL
    ## no individual level variation
    #sel.i <- sample(nrow(mu),1)
    ## change 11/2/2020 DL
    ## coefficients were linked to residence
    ## not motility (i.e. negated)
    tracks[[i]] <- ctmc.sim.study(
      stack=cov__,coef=-1*beta0,start = c(165,-20),
      crw = crw0)
    tracks[[i]]$value <- ctmc.data(tracks[[i]])
  }
  
  ## subsample long tracks
  for(i in 1:length(tracks)){
    all_t <- sort(unique(tracks[[i]]$value$t))
    start_ <- sample(all_t[1:floor(length(all_t)/2)],1)
    len_ <- 1+rpois(n=1,lambda = ndays)
    tracks[[i]]$sub <- subset(
      tracks[[i]]$value,t %in% seq(start_,start_+len_))  
  }
  
  i <- 2
  sapply(tracks,function(x) length(unique(x$sub$t)))
  
  ## simulation fishery observations
  
  ## scalar from motility to intensity
  #delta <- c(1,1.5,-1.14,-1.14,1.31,-1.61)
  delta <- c(1,-1,-1,-1,-1)
  
  ## population level resource selection intensity coefficients
  beta <- beta0*delta
  
  ## additive bias in fishery observations
  e <- cov__$bathy[[1]]
  pts <- as.data.frame(e,xy=T)
  e[cellFromXY(e,pts[,1:2])] <- with(pts,0.7*(3-((0.5*(x-160))^2+(6*(y+25))^2)/2500))
  #plot(e)
  
  ## scale down the number of fishery observations
  #beta[1] <- -7
  beta[1] <- -10.5
  
  ##  daily intensity
  lambda0 <- residence(beta,cov__)$fit
  
  ## apply additive bias
  lambda <- vector("list",dim(lambda0)[3])
  for(i in 1:length(lambda)){
    lambda[[i]] <- exp(log(lambda0[[i]])+e)
  }
  lambda <- rast(lambda)
  
  ## convert to monthly units
  lambda <- 30*lambda
  
  ## simulate observations
  source("rppm_3.R")
  
  ## expected number of events
  if(FALSE){
    tmp <- as.array(lambda)
    ss <- apply(tmp,3,function(x) sum(x,na.rm = T))
    area <- prod(res(lambda))
    barplot(ss*area,xlab="month",ylab="#",
            main="expected observations")
  }
  
  cat("simulating biased fishery observations.\n")
  loc <- rppm.st(lambda)
  
  save(loc,tracks,beta,cov__,lambda0,beta0,delta,file=fn)
  
  tmp <- vector("list",length(tracks))
  for(i in 1:length(tracks)){
    sub_ <- tracks[[i]]$sub
    names(sub_) <- c("long","lat","aday","amon")
    sub_$id <- i
    tmp[[i]] <- sub_
  }
  track.df <- do.call(rbind,tmp)
  write.csv(track.df,file="simu2/tracks.csv",row.names = F)
  write.csv(loc,file="simu2/loc.csv",row.names = F)
  
  library(ggplot2)
  ggplot(track.df,aes(x=long,y=lat))+geom_line()+
    facet_wrap(~id)+xlim(c(173,175))
  
  #174.875	-35.875
  
  
  ggplot(loc,aes(x=long,y=lat))+geom_point()+
    facet_wrap(~amon)
}

block1_build_seascape <- function(){
  rm(list = ls())
  ## load effect sizes
  load(file="../Analyses/easi_fish/bsdm5A3T203B1.rda")
  beta <- coef(bsdm.obj$ctds)
  crw <- coef(bsdm.obj$ctds)["crw"]
  beta <- beta[c(1:3,5,9,10)]
  beta <- beta[c(1,3,6,2,4,5)]
  beta[3] <- beta[3]/6.5 ## normalize the results
  beta[5] <- beta[5]/0.24
  
  ## load sea scape variables
  library(terra)
  bathy0 <- rast("../GIS/bathy_025deg.tif")
  bathy <- -1*bathy0[[1]]/1000
  
  sst <- rast("../GIS/SST_2000-01_2024-01.tif")
  chl <- rast("../GIS/CHL_2000-01_2024-01.tif")
  ssh <- rast("../GIS/SSH_2000-01_2024-01.tif")

  library(ctmcmove)
  
  
  ## crop to smaller study areas
  ext <- extent(150,175,-40,-10)
  sst.c <- crop(sst,ext)
  bathy.c <- crop(bathy,ext)
  ssh.c <- crop(ssh,ext)
  chl.c <- crop(chl,ext)
  
  ## subset of data
  sel <- 1:10 ## for first ten months
  bathy.c <- bathy.c[[sel]]
  sst.c <- sst.c[[sel]]
  ssh.c <- ssh.c[[sel]]
  chl.c <- chl.c[[sel]]
  
  ## pre-process to fill in missing values
  library(gstat)
  library(sf)
  source("idw_rast_3.R")
  base <- bathy.c
  sst.c <- idw.rast.stack(sst.c,base)
  ssh.c <- idw.rast.stack(ssh.c,base)
  chl.c <- idw.rast.stack(chl.c,base)
  
  lst0 <- vector("list",nlyr(sst.c))
  for(i in 1:nlyr(sst.c)){
    lst0[[i]] <- bathy.c
  }
  
  ## organize into a list
  cov__ <- list(
    bathy=rast(lst0),
    sst=sst.c,
    chl=chl.c,
    ssh=ssh.c)
  
  writeRaster(rast(lst0),file="simu2/bathy.tif",overwrite=TRUE)
  writeRaster(sst.c,file="simu2/sst.tif",overwrite=TRUE)
  writeRaster(chl.c,file="simu2/chl.tif",overwrite=TRUE)
  writeRaster(ssh.c,file="simu2/ssh.tif",overwrite=TRUE)
  
  ## collinearity
  cor(as.vector(cov__[[1]]),as.vector(cov__[[2]]),use ="complete.obs")
  

  ## scale down the number of fishery observations
  beta[1] <- -4
  
  ## scratch
  ## determine bias coefficients so that bias 
  ## between -100% and 100% with and without coefficients
  
  source("residence_5.R")
  library(abind)
  
  ## without bias
  sptw.lambda <- residence(beta[-length(beta)],cov__)
  
  
  look <- function(x,y){
    l.x <- as.array(x$fit/y$fit)
    l.x <- as.vector(l.x) - 1
    l.x <- na.omit(l.x)
    print(summary(l.x))
    hist(l.x,xlab="ratio",
         main="relative bias due to fishery observation")
  }
  
  
  save(beta,crw,file="simu2/coef.rda")
  
}