all_group <- function(formula,verbose=T){
  var_ <- all.vars(formula)[-1]
  id_ <- seq(1,length(var_))
  lst_ <- vector("list",length(var_))
  for(i in 1:length(var_)){
    lst_[[i]] <- seq(1,length(var_))
  }
  grid_ <- expand.grid(lst_)
  grid2_ <- apply(grid_,1,function(x){
    as.numeric(factor(x,levels = unique(x)))
  })
  key_ <- apply(grid2_,2,paste0,collapse="")
  grid3_ <- grid2_[,!duplicated(key_)]
  if(verbose){
    return(apply(grid3_,2,function(x){
      split(var_,x)
    }))
  }else{
    return(t(grid3_))
  }
}

all_group_label <- function(x){
  lst <- split(seq(1,length(x)),x)
  label <- sapply(lst,paste0,collapse="-")
  paste(label,collapse = ", ")
}