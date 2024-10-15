extract.helper <- function(object,layer){
  ## object: a list of raster stacks
  ## layer: layer id from each stack
  ## value: a raster stack objects with each environmental layer
  
  out <- vector("list",length(object))
  for(i in 1:length(out)){
    out[[i]] <- raster(object[[i]],layer=layer)
  }
  names(out) <- names(object)
  stack(out)
}