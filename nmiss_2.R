## helper script to subset to data where the formula is not missing.
nmiss <- function(formula,data)
{
  if(is.null(data)){
    return(data)
  }
  tmp <- na.omit(data[,all.vars(formula)])
  #browser()
  if(!is.null(attr(tmp,"na.action"))){
    data <- data[-attr(tmp,"na.action"),]
  }
  data
}

test <- function()
{
  setwd("C:/Users/dliang/Google Drive/Spring_2017/Turtle_Watch/Analysis")
  df <- data.frame(x=rnorm(100),y=rnorm(100))
  df$x[1:2] <- NA
  source("nmiss_1.R")
  df2 <- nmiss(y~x,df)
  df3 <- nmiss(y~x,df2)
}