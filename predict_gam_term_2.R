## helper function to predict from fitted gam object
## according to terms
## formula: a gam formula of model
## x: a list of character vector defining the "variables" to combine
## beta: current estimate of coefficients
## X: lpmatrix output from gam fitted object
## value: constructed terms according to x
predict.gam.term <- function(formula,x,beta,X){
  ## name of terms
  name_ <- term.name(formula)
  terms_ <- term.c(name_,arg=FALSE)
  all_terms <- name_$term[attr(name_$arg,"init"):length(name_$term)]
  
  lst <- vector("list",length(x))
  for(i in 1:length(x)){
    curr_terms <- terms_[match(x[[i]],all_terms)]
    curr_col <- lapply(curr_terms,grep,x=colnames(X))    
    curr_col <- do.call(c,curr_col)
    lst[[i]] <- as.numeric(X[,curr_col,drop=FALSE] %*% beta[curr_col])
  }
  
  out <- as.data.frame(do.call(cbind,lst))
  names(out) <- paste0("E.",1:ncol(out))
  attr(out,"definition") <- x
  
  out
}
dev_ <- function(){
 load(file="Scratch/predict_gam_term_1.rda")
  X <-predict(object,quad,type="lpmatrix") 
  e <- predict.gam.term(
    formula=object$formula,
    x=list("bathy","sst",c("fpi","ssh")),
    beta=bmu,
    X=X)
}