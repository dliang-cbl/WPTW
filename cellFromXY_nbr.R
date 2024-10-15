## extract cell number from xy given a distance threshold
## object/xy: see documentation from cellFromXY
## eps : distance threshold
## value: cell ID of identified
require(FNN)
cellFromXY.nbr <- function(object,xy,eps=NULL)
{
  cell <- cellFromXY(object,xy)
  if(!is.null(eps)){
    ## identify any points outside the area
    l.na <- which(is.na(cell))
    ## identify the nearest cell
    if(length(l.na)){
      l.xy <- xy[l.na,]
      l.knn <- get.knnx(query=l.xy,data=coordinates(object),k=1)
      l.r <- l.knn$nn.index[,1]
      ## remove if the nearest cell beyond
      l.r[l.knn$nn.dist[,1]>eps] <- NA
      cell[l.na] <- l.r
    }
  }
  cell
}