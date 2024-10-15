## loglik function
## x: the quasi-likelihood fit
## data: the observed and quadrature points
## value: log-likelihood of the fitted poisson process
## adjust the total area so it match
ppm.st.lik <- function(x,data){
  lambda <- fitted(x)
  l.lik <- with(subset(data,!is.na(z)),{
    wt * (z*log(lambda)-lambda)
  })
  sum(l.lik)
}

## value: residual measure of each quadrature point
ppm.st.residuals <- function(x,data){
  l.data <- subset(data,!is.na(z))
  lambda <- predict(x,newdata=l.data,type="response")
  l.res <- with(l.data,{
    wt*(z-lambda)*1/sqrt(lambda)
  })
  l.res
}
test <- function()
{
  ppm.st.lik(gam0,glm0)
  ## -200.7945
}

## x: predicted data with residuals
## formula: define residual on left and area on right and time period on right
## value:
## a residual plot
## plus the aggregated time
plot.ppm.st.residuals <- function(formula,x,nclass=30,cex=3,pch=21,bg="steelblue1",...)
{
  x_ <- x[,all.vars(formula)[3]]
  ux_ <- unique(x_)
  if(length(ux_)>nclass){
    #browser()
    ## discretize if the predictor is continous
    br_ <- quantile(x_,probs=seq(0,1,len=nclass+1))
    cx_ <- cut(x_,breaks=unique(br_),label=FALSE)
    cx_ <- unique(br_)[cx_]
    x[,all.vars(formula)[3]] <- cx_
  }
  tmp1 <- aggregate(formula,data=x,FUN=sum)
  vars <- all.vars(formula)
  wt <- tmp1[,vars[1]]
  upr <- 2*sqrt(wt)
  lwr <- -2*sqrt(wt)
  e <- tmp1[,vars[2]]
  t_ <- tmp1[,vars[3]]
  plot(t_,e,ylim=range(e,upr,lwr),
       ylab="Pearson Residual",...)
  polygon(x=c(t_,rev(t_)),y=c(lwr,rev(upr)),col=grey(0.9))
  points(t_,e,cex=cex,pch=pch,bg=bg)
  lines(lowess(t_,e),col="red")
  tmp1
}

