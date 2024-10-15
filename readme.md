Monthly prediction using both tracking and fishery observations
================
Dong Liang @ CBL UMCES
2024-10-15

This file documents the estimation of a joint model described in the
manuscript and the prediction of monthly intensity based on the joint
model. Simulated tracks and observations were used to demonstrate the
codes.

``` r
rm(list=ls())
options(warn = 1)
library(raster); library(ctmcmove); library(mgcv);
options("rgdal_show_exportToProj4_warnings"="none")
```

## Environmental data

This section assumes that environmental variables have been downloaded
manually as raster files onto the current work directory of R. In
addition, these files have been re-sampled to the same geographic extent
and same resolution. The longitude and latitude coordinate system is
assumed on all input rasters. The geographic extent used was between 150
and 175 degrees in longitude and between -40 and -10 degrees in
latitude.

``` r
list.files(pattern="tif")
```

    ## [1] "bathy.tif"     "chl.tif"       "intensity.tif" "ssh.tif"       "sst.tif"

The input rasters are then loaded into a <i> list </i> format required
by modeling codes. Currently, the monthly resolution is assumed for the
environmental rasters.

``` r
seascape <- list(bathy=stack("bathy.tif"),
              sst=stack("sst.tif"),
              chl=stack("chl.tif"),
              ssh=stack("ssh.tif"))
dim(seascape[[1]])
```

    ## [1] 120 100  10

## Track Data Processing

This section takes the simulated tracks and processes the data into a
structure needed for modeling. The simulated tracks are at the daily
resolution, also with a variable linking with the environmental raster
layer <i> amon </i>.

``` r
## load tracking data codes
source("nmiss_2.R")
source("s2inla_6.R")
source("ctmc2glm_study_3.R")

# load simulated tracks
tracks0 <- read.csv("tracks.csv")
str(tracks0)
```

    ## 'data.frame':    511 obs. of  5 variables:
    ##  $ long: num  175 175 175 175 175 ...
    ##  $ lat : num  -35.9 -35.6 -35.9 -36.4 -36.4 ...
    ##  $ aday: int  141 142 143 144 145 146 147 148 149 150 ...
    ##  $ amon: int  5 5 5 5 5 5 5 5 5 6 ...
    ##  $ id  : int  1 1 1 1 1 1 1 1 1 1 ...

``` r
tracks <- split(tracks0,~id)
```

For each individual, the daily tracks were converted into a continuous
time discrete space path, consisting of the residential cells according
to the environmental raster, the residence time in terms of days, and
the transition between cells.

``` r
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
```

To include random effects at an individual level, the output also
includes the individual <i> id </i>. The seascape variables were renamed
to facilitate modeling.

``` r
track.data <- do.call(rbind,tmplist)
track.data$id <- factor(track.data$id)
names(track.data)[1+1:length(seascape)] <- names(seascape)
```

## Fishery Observation Processing

This section processes the monthly observations into a structure needed
for joint models. The simulated observations are at monthly resolution.

``` r
loc <- read.csv("loc.csv")
str(loc)
```

    ## 'data.frame':    434 obs. of  3 variables:
    ##  $ long: num  165 160 160 160 157 ...
    ##  $ lat : num  -24.1 -26.1 -15.6 -16.4 -13.6 ...
    ##  $ amon: int  1 1 1 1 1 1 1 1 1 1 ...

Here <i> amon </i> should correspond to the layers in the seascape <i>
stack </i>.

Pseudo-absence data are generated for each month to estimate the
intensity of fishery observations.

``` r
## load point process codes
source("cellFromXY_nbr.R")
source("ppm_st_glm_6d.R")
source("ppm_st_utils.R")


## study area raster layer (NA=outside area, non-NA=in the area)
sp <- raster(seascape[["bathy"]],layer=1)

# generate pseudo-absence data 
args(ppm.st.glm)
```

    ## function (x, window, covariates, trend = ~1, coordinate = ~x + 
    ##     y + t, eps = NULL, newdata = NULL, newdata.ctime = NULL, 
    ##     nfac = 2, verbose = TRUE) 
    ## NULL

``` r
pa <- ppm.st.glm(
  loc, 
  window = sp, #a R raster* object with NA denoting masked area and 1 the study area,
  covariates = seascape,
  trend=reformulate(names(seascape)),
  coordinate=~long+lat+amon)
```

    ## pseduo-absence at time=  1 
    ## pseduo-absence at time=  2 
    ## pseduo-absence at time=  3 
    ## pseduo-absence at time=  4 
    ## pseduo-absence at time=  5 
    ## pseduo-absence at time=  6 
    ## pseduo-absence at time=  7 
    ## pseduo-absence at time=  8 
    ## pseduo-absence at time=  9 
    ## pseduo-absence at time=  10 
    ## j=  1 
    ## j=  2 
    ## j=  3 
    ## j=  4 
    ## j=  5 
    ## j=  6 
    ## j=  7 
    ## j=  8 
    ## j=  9 
    ## j=  10

## Joint Modeling

The joint model involves linking the motility from tracks with the
environmental condition, here using the seascape conditions as
covariates. This model component matches that of Hoover et al. 2019.

``` r
library(abind)
source("all_group_1.R")
source("predict_gam_term_2.R")
source("sdm_bc_5d_par.R")
source("sim2jam_study_1.R")
source("term_utils.R")

## load parallel code for Monte Carlo
library(doParallel)
nodes <- 7 ## define number of nodes
if(F){
  cl <- makeCluster(nodes)
  dump <- clusterEvalQ(cl,{
    library(mgcv)
    source("predict_gam_term_2.R")
    source("term_utils.R")
  })
  registerDoParallel(cl)
}
```

Change the <i> F </i> to <i> T </i> above if modeling is conducted on a
high-performance computer, and replace the <i> nodes </i> from 7 with
the actual number of nodes available.

``` r
## load prediction codes
library(doParallel)
source("predict_sdm_bc_5d.R")
source("ppm_st_glm_6f.R")
source("extract_1.R")
source("intensity_2.R")

## telemetry model: linking motility to seascape variables
term0 <- names(seascape)
for0 <- reformulate(term0,response="z")
for0
```

    ## z ~ bathy + sst + chl + ssh

Here motility is linked with four seascape variables as linear terms.
Note smooth terms are allowable here as well.

The bias in fishery observation is assumed to be a spatial term,
representing the observational process of the fisheries, this bias, plus
a scaled version of the motility from the movement data, is used to
jointly model the intensity of the animal occurrences, without possible
observational biases.

``` r
for1 <- ~s(long,lat,bs="ds",m=c(1,0.5),k=10)
```

The scaling of the motility term to intensity was defined based on a set
of integers of length equal to the number of predictors. All predictors
with the same integer share the same scales. In below, only one scaling
term was assumed, i.e. the motility coefficients were multiplied by a
single scalar to model intensity. This single scalar parameter has a
prior estimate of -1, thus higher morality is associated with lower
intensity.

``` r
scale_terms <- all_group(
  reformulate(term0,response = "z"),verbose=F)
head(scale_terms)
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    1    1    1
    ## [2,]    1    2    2    2
    ## [3,]    1    2    1    1
    ## [4,]    1    1    2    2
    ## [5,]    1    2    3    3
    ## [6,]    1    1    2    1

``` r
term1 <- scale_terms[1,]
```

The joint model is estimated using Monte Carlo methods below.

``` r
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
save(bsdm.obj,file="bsdm_obj.rData")
```

The arguments <i> d </i> and <i> S </i> control the prior scales linking
motility and intensity, <i> niter </i> indicates the number of Monte
Carlo simulations to estimate the joint model, this number should be
large (at least 199) to reduce Monte Carlo errors. <i> timefact </i>
denotes the ratio of temporal units in the tracks (daily) and the
observations (monthly), here a 30-day month was assumed. Due to the
Monte Carlo simulation, this script might take up to 20 hours to finish,
thus the results were saved as a <i> rData </i> file for offline
processing.

If using parallel package, be sure to stop the cluster.

``` r
if(F){
  stopCluster(cl)
}
```

## Cross validation

The movement component and the point process components of the joint
modeling were evaluated using cross-validation. Cross-validation for the
movement components was done in two ways. First, individual tracks were
iteratively held out and predicted using a movement model estimated from
all the other tracks. Secondly, the first half of all tracks were used
to train a movement model, which was used to predict the second half of
the tracks. Given that movement data were represented using a sequence
of binary actions among residence cells along the path, we quantified
the predictive performance of the model using sensitivity, specificity,
and area under the receiver operating characteristic (ROC) curve. The
point process component was cross-validated in time: the model was
trained using the first half of the observations and used to predict the
fishery observation during the second half of the times series.
Prediction was conducted iteratively forward by months. Perason
residuals were calculated based on the fitted point process model to
quantify the predictive capability of the point process component.

Cross-validation of the movement model involves many repeated runs of
the modeling code. If possible, these runs should be conducted in
parallel. We first define a cluster object and load the necessary
scripts.

``` r
source('cv_ctds_2b.R')
if(F){
  cl <- makeCluster(nodes)
  dump <- clusterEvalQ(cl,{
    library(mgcv)
    source("predict_sdm_bc_5d.R")
    source("predict_gam_term_2.R")
    source("term_utils.R")
  })
  registerDoParallel(cl)
}
```

Next, the cross-validation can be conducted among tracks, along with the
ROC curve.

``` r
library(ThresholdROC)
## cross-validate movement component
## by individual
cv1 <- cv.ctds(bsdm.obj$ctds,track.data,type=1,plotit = F)
```

    ## cross-validate by individuals...
    ## 1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,done.
    ## Sensitivity= 0.784 , Specificity= 0.585

``` r
plotCostROC(cv1,which.plot = "roc")
```

<img src="readme_files/figure-gfm/unnamed-chunk-18-1.png" width=".5\linewidth" style="display: block; margin: auto;" />

Sensitivity indicates the proportion of movements correctly predicted by
the model. Specificity indicates the proportion of non-movements
correctly predicted by the model. The area under the ROC curve provides
a single index of the classification capability of the movement model
based on the environmental data.

Cross-validation can also be conducted within tracks but forward in
time.

``` r
## by time within each individuals
cv2 <- cv.ctds(bsdm.obj$ctds,track.data,type=2,plotit = F)
```

    ## cross-validation by time within individuals.
    ## 1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,Sensitivity= 0.841 , Specificity= 0.573

``` r
plotCostROC(cv2,which.plot = "roc")
```

<img src="readme_files/figure-gfm/unnamed-chunk-19-1.png" width=".5\linewidth" style="display: block; margin: auto;" />

Likewise, cross-validation of the point process model can make use of
the parallel nodes.

``` r
source("bsdm_utils_1b.R")
source("ppm_st_utils.R")
edat2 <- ppm.cv(bsdm.obj,pa,tcoord="amon",wcoord = "wt")
names(edat2)
```

    ##  [1] "cell"  "long"  "lat"   "bathy" "sst"   "chl"   "ssh"   "amon"  "wt"    "z"     "crw"   "tau"  
    ## [13] "fit"   "sd"    "resid"

The output contains the actual fishery observations, the fitted
intensity of the bias-corrected observations, along with the Pearson
residuals. These residuals should remain random to known predictor
variables included in the modeling. A residual plot versus these
predictors thus can graphically diagnose the point process model.

``` r
plot.ppm.st.residuals(cbind(wt,resid)~amon,edat2,
                      cex=3,pch=21,
                      bg="steelblue1",xlab="month")
```

    ##   amon  wt  resid
    ## 1    5 713  -5.54
    ## 2    6 713   6.91
    ## 3    7 713  25.86
    ## 4    8 713 -29.51
    ## 5    9 713 -26.37
    ## 6   10 713   5.48

<img src="readme_files/figure-gfm/unnamed-chunk-22-1.png" width=".5\linewidth" style="display: block; margin: auto;" />

The residuals are around zero according to month, with the 95%
confidence interval containing the observed residuals, this suggests no
lack of fit of the point process model.

The following code generates a similar plot for another predictor.

``` r
plot.ppm.st.residuals(cbind(wt,resid)~bathy,edat2,
                      cex=3,pch=21,
                      bg="steelblue1",xlab="bathy")
```

    ##    bathy    wt   resid
    ## 1  -5.73 143.2 -13.768
    ## 2  -5.20 140.6  20.111
    ## 3  -4.89 324.8   2.118
    ## 4  -4.83 103.9  -4.537
    ## 5  -4.71 145.5  -3.944
    ## 6  -4.42 139.9 -18.472
    ## 7  -4.21 142.1  -9.802
    ## 8  -3.98 142.5  -4.700
    ## 9  -3.68 183.0  -4.663
    ## 10 -3.60 102.0   0.586
    ## 11 -3.53 149.6  14.606
    ## 12 -3.34 247.9  32.116
    ## 13 -3.22  28.9  -0.856
    ## 14 -3.20 142.5  12.115
    ## 15 -3.10 143.2  -0.370
    ## 16 -2.92 141.8 -10.978
    ## 17 -2.85 142.1  30.971
    ## 18 -2.66 142.5  -0.390
    ## 19 -2.53 142.5  -1.035
    ## 20 -2.33 143.2  -6.301
    ## 21 -2.22 141.4   7.007
    ## 22 -2.00 142.9 -15.361
    ## 23 -1.85 143.2 -24.451
    ## 24 -1.66 142.5   0.155
    ## 25 -1.46 142.9 -12.402
    ## 26 -1.30 146.6  -0.823
    ## 27 -1.06 139.1   5.615
    ## 28 -0.66 142.5 -14.479
    ## 29 -0.18 142.5  -0.559

<img src="readme_files/figure-gfm/unnamed-chunk-23-1.png" width=".5\linewidth" style="display: block; margin: auto;" />

## Statistical Inferences

Population-level coefficients from the movement model can be summarized
below.

``` r
summary(bsdm.obj$ctds)
```

    ## 
    ## Family: quasi 
    ## Link function: log 
    ## 
    ## Formula:
    ## z ~ bathy + sst + chl + ssh + crw + offset(log(tau))
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.5427     0.5054    3.05   0.0023 ** 
    ## bathy        -0.0619     0.0412   -1.50   0.1332    
    ## sst           0.0981     0.0308    3.18   0.0015 ** 
    ## chl           0.3669     0.1514    2.42   0.0155 *  
    ## ssh          -1.7674     0.6570   -2.69   0.0072 ** 
    ## crw           0.5906     0.0395   14.93   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## 
    ## R-sq.(adj) =  -0.556   Deviance explained = 10.6%
    ## GCV = 0.78558  Scale est. = 0.74286   n = 2047

The tracks show strong directional dependence, and resource selection in
temperature, chlorophyll-a and sea surface height. The aggregated
effects of each seascape variable to the partial effect of the intensity
can be shown below.

``` r
source("plot_sdm_bc_5d.R")
source("predict_gam_term_2.R")
sm1 <- plot_sdm_bc(bsdm.obj,select=1)
```

    ## averaging  sst 
    ## averaging  chl 
    ## averaging  ssh 
    ## averaging  crw 
    ## averaging  tau

``` r
sm2 <- plot_sdm_bc(bsdm.obj,select=2)
```

    ## averaging  bathy 
    ## averaging  chl 
    ## averaging  ssh 
    ## averaging  crw 
    ## averaging  tau

``` r
sm3 <- plot_sdm_bc(bsdm.obj,select=3)
```

    ## averaging  bathy 
    ## averaging  sst 
    ## averaging  ssh 
    ## averaging  crw 
    ## averaging  tau

``` r
sm4 <- plot_sdm_bc(bsdm.obj,select=4)
```

    ## averaging  bathy 
    ## averaging  sst 
    ## averaging  chl 
    ## averaging  crw 
    ## averaging  tau

<img src="readme_files/figure-gfm/unnamed-chunk-25-1.png" width=".5\linewidth" style="display: block; margin: auto;" /><img src="readme_files/figure-gfm/unnamed-chunk-25-2.png" width=".5\linewidth" style="display: block; margin: auto;" /><img src="readme_files/figure-gfm/unnamed-chunk-25-3.png" width=".5\linewidth" style="display: block; margin: auto;" /><img src="readme_files/figure-gfm/unnamed-chunk-25-4.png" width=".5\linewidth" style="display: block; margin: auto;" />
Scaling parameters linking estimated motility effects to the intensity
of the thinned Point Process Model can be summarized below.

``` r
source("scalar_bsdm_1.R")
scalar.bsdm(bsdm.obj)
```

    ##                       fit    lwr   upr
    ## bathy=sst=chl=ssh -0.0428 -0.703 0.508

<img src="readme_files/figure-gfm/unnamed-chunk-26-1.png" width=".5\linewidth" style="display: block; margin: auto;" />

The common scalar links each movement term to the intensity component of
the model, the above figure shows the posterior estimates of this scalar
to be most negative. Thus motility is negatively associated with
intensity.

The spatial random effects are mapped to capture the potential
observational process.

``` r
source("s2raster_1.R")
source("intensity_2.R")
library(RColorBrewer)
x <- plot(bsdm.obj$ppm,select = 0,scheme=2,
          contour.col=NA,rug = FALSE)
base <- raster(seascape[[1]],layer=1)
sm_term <- s2raster(base,x)
smPal <- brewer.pal(6,"BrBG")
plot(sm_term,col=smPal)
```

<img src="readme_files/figure-gfm/unnamed-chunk-27-1.png" width=".65\linewidth" style="display: block; margin: auto;" />

## Prediction

The above training procedure is connected with the raw track and
observation data and should be conducted whenever new raw data becomes
available. The trained model can be used to process incoming
environmental data and generate predictive maps of interaction hot spot.

The fitted model is dumped into a <i> RData </i> file. We next load the
fitted object.

``` r
load(file="bsdm_obj.rData")
```

The object name is <i> bsdm.obj </i>, this object is trained from all
available fishery and tracking observations.

Depending on the resolution of the raster, prediction could be
computationally intensive. The prediction code can be run in parallel if
necessary. If parallel is needed, change the following block from <i>
FALSE </i> to <i> TRUE </i> and define <i> nodes </i> as the number of
nodes available on the server.

``` r
## prepare parallel nodes
library(doParallel)
if(FALSE){
  cl <- makeCluster(nodes)
  dump <- clusterEvalQ(cl,{
    source("ppm_st_glm_6f.R")
    source("extract_1.R")
    source("intensity_2.R")
    source("predict_gam_term_2.R")
    source("term_utils.R")
    source("predict_sdm_bc_5d.R")})
  registerDoParallel(cl)
}
```

Now we define the area of prediction. For demonstration, the same
environmental data are used for estimating the model as well as
prediction. But new environmental data should be used here for
prediction.

``` r
newdata <- list(bathy=stack("bathy.tif"),
                sst=stack("sst.tif"),
                chl=stack("chl.tif"),
                ssh=stack("ssh.tif"))
```

The incoming seascape condition is converted into a <i> data.frame </i>
with corresponding environmental conditions, the variable <i> amon </i>
defines the layer (month) name, but it was not used in the month by
month prediction.

``` r
newd <- ppm.st.glm.new(
  window=sp,newdata=newdata,
  coord=~long+lat+amon)
```

The data frame <i> newd </i> contains the seascape variables and the
geographical coordinates, as well as additional variables necessary to
fit the joint Poisson Process Model.

Now we can generate prediction.

``` r
cat("predicting monthly intensity...\n")
(ts1 <- proc.time())
raw.pred <- predict_sdm_bc(bsdm.obj,newd,all=F,
                           verbose=F,batch=10)
(ts2 <- proc.time())
```

The <i> all </i> argument of <i> predict_sdm_bc </i> indicates whether
the bias field is used in the prediction. It was turned off to remove
the observation process in the predictive mapping of intensity. And <i>
batch </i> denotes processing predictions in 10 batches, possibly in
parallel.

## Visualization

The output can be converted into an <i> raster </i> object.

``` r
predMonth <- ppm.st.glm.fill(
  window = sp,
  newdata = raw.pred,
  formula = fit~long+lat+amon,
  verbose = T,
)
```

    ## processing  deptho.1 
    ## processing  deptho.10 
    ## processing  deptho.2 
    ## processing  deptho.3 
    ## processing  deptho.4 
    ## processing  deptho.5 
    ## processing  deptho.6 
    ## processing  deptho.7 
    ## processing  deptho.8 
    ## processing  deptho.9

The output in <i> raw.pred </i> was converted into a <i> raster stack
</i>. Specifically the posterior mean in <i> fit </i> was geo-coded onto
the raster <i> base </i> according to the geographical coordinates <i>
point_x </i>, and <i> point_y </i>, with layer names in <i> amon </i>.

The unit in the output is originally in 0.25 squared decimal degrees.

``` r
res(sp)
```

    ## [1] 0.25 0.25

``` r
output <- raster(predMonth,layer=1)
```

The following code visualize the first layer.

``` r
library(raster)
library(sf)
library(tmap)
library(tmaptools)
tmap_mode("plot")
#library(OpenStreetMap)
library(RColorBrewer)
library(viridis)

summary(values(output))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##       0       0       0       0       0       0     592

``` r
col <-brewer.pal(n=5,"RdPu")
tm_shape(output)+tm_raster(
  palette=col,
  breaks = c(-999,0.043,0.044,0.045,0.046,0.048,999),
  labels = c("0 to 0.043","0.043 to 0.044",
             "0.044 to 0.045","0.045 to 0.046",
             "0.046 to 0.048", "0.048 above"),
  title=expression(lambda))+
  tm_layout(legend.bg.color = "white")
#save.image(file=fn)
```

<img src="readme_files/figure-gfm/unnamed-chunk-36-1.png" width=".5\linewidth" style="display: block; margin: auto;" />

The output can be saved as a <i>tif </i> file <i>intensity.tif </i>
under the current directory and visualized manually.

``` r
writeRaster(output,file=paste0("intensity.tif"),overwrite=T)
```
