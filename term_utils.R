term.c <- function(x,arg=T){
  ## x: output from term.name
  ## value: character vector of terms
  is.f_ <- nchar(x$func)>0
  if(!arg){
    tmp_ <- with(x,paste0(func,"\\(",term,"\\)"))
    term_ <- with(x,ifelse(nchar(func)>0,tmp_,term))
    sel <- attr(x$arg,"init")+seq(0,length(x$arg)-1)
    return(term_[sel])
  }
  o_ <- rep("",length(x$term))
  o_[!is.f_] <- x$term[!is.f_]
  
  ## convert argument to character vectors
  a_ <- sapply(x$arg,function(lst){
    r <- ""
    if(length(lst)>0){
      u <- rep("",length(lst))
      for(i in 1:length(lst)){
        aa_ <- lst[[i]]
        if(is.name(aa_)){
          aa_ <- as.character(aa_)
        }
        else if(is.character(aa_)){
          aa_ <- paste0('"',aa_,'"')
        }
        u[i] <- paste0(names(lst)[i],"=",aa_)
      }
      r <-  paste0(u,collapse=",")
    }
    r
  })
  a1_ <- rep("",length(x$term))
  a1_[attr(x$arg,"init")+seq(0,length(x$arg)-1)] <- a_
  
  is.a_ <- is.f_ & nchar(a1_)>0
  is.na_ <- is.f_ & nchar(a1_)==0
  o_[is.a_] <- paste0(
    x$func[is.a_],"(",x$term[is.a_],",",a1_[is.a_],")")
  o_[is.na_] <- paste0(
    x$func[is.na_],"(",x$term[is.na_],")"
  )
  o_
}
test.c <- function(){
  debug(term.c)
  term.c(term.name(z~s(x)+y+offset(log(tau))))
  term.c(term.name(~s(x)+y))
  term.c(term.name(~s(x)+s(y)))
  term.c(term.name(~s(x,bs="cr",k=5)+s(y)))
  term.c(term.name(z~s(x)+y))
  term.c(term.name(z~s(x)+s(y)))
  term.c(term.name(z~s(x,bs="cr",k=5)+s(y)))
  term.c(term.name(z~s(x,bs="cr",k=5)+s(y,k=10)))
}
term.name <- function(formula,ops=c("+","~")){
  ## ops: a list of binary operators
  
  ## extract the terms and functions/transformations
  ## from a formula
  
  n_ <- all.names(formula)
  v_ <- all.vars(formula)
  
  ## check for duplicate terms
  v0_ <- all.vars(formula, unique = FALSE)
  if(length(v0_)>length(v_)){
    stop("duplicate term found, please check.")
  }
  ## index  of all variables
  i_ <- match(v_,n_)
  ## index of all operators
  p_ <- which(n_%in% ops) 
  ## index of all function names
  fidx_ <- rep(T,length(n_))
  fidx_[c(i_,p_)] <- FALSE
  f_ <- which(fidx_)
  
  ## Now to match function names to the terms
  ## Start with the last operator
  ## Now to match function names to the terms
  ## Start with the last operator
  last_ops <- p_[length(p_)]
  ## find the following two variables if available
  if(length(i_)>0){
    if(sum(i_>last_ops)>1){
      pair_last <- i_[i_>last_ops][1:2]
    }else{
      pair_last <- i_[i_>last_ops]
    }
  }else{
    return(list(term="",func=""))
  }
  
  ## helper function to find preceding names of variables
  helper <- function(x,y,z){
    ## x: location of operator
    ## y: locations of variables
    ## z: all locations of function names
    ## value: find the preceding function of each variable
    ##        if available
    o <- rep(NA,length(y))
    for(i in 1:length(y)){
      prev_ <- y[i]-1
      if(prev_>x){## is not an operator
        if(prev_ %in% z){ ## and is a function name
          o[i] <- prev_
        }
      }
    }
    o
  }
  ## index of function names to each variable
  o <- rep(NA,length(i_)) 
  o[match(pair_last,i_)] <- helper(last_ops,pair_last,f_)
  ## index of unprocessed names
  uname_ <- i_[!(i_%in%pair_last)]
  
  for(i in (length(p_)-1):1){ 
    ## find the previous operator
    last_ops <- p_[i]
    ## and the following un-processed variable
    next_ <- uname_[uname_>last_ops][1]
    if(!is.na(next_)){
      ## find the corresponding function
      o[match(next_,i_)] <- helper(last_ops,next_,f_)
      ## remove the variable as processed
      uname_ <- uname_[!(uname_==next_)]
    }
  }
  
  ## output function
  r <- rep("",length(o))
  r[!is.na(o)] <- n_[o[!is.na(o)]]
  
  #browser()
  ## output arguments
  rhs_ <- term.rhs(formula)
  removechar <- function(x,y){
    ## remove y from a list of characters
    ux <- unlist(x)
    ii <- !(ux %in% y)
    x[ii]
  }
  
  i0 <- length(r) - length(rhs_) ## last response term
  a_ <- vector("list",length(rhs_))
  for(i in 1:length(rhs_)){
    lang_ <- as.list(str2lang(rhs_[i]))
    fn_ <- c(r[i0+i],v_[i0+i])
    a_[[i]] <- removechar(lang_,fn_)
  }
  attr(a_,"init") <- i0+1
  list(term=v_,func=r,arg=a_)
}
test.name <- function(){
  ## need to make sure no duplicate
  term.name(y~s(x,bs="cr",k=5))
  term.name(y+z~s(a)+s(b)+d) 
  term.name(y~s(x)+s(z)+s(w))
  term.name(y~x+z+w)
  term.name(cbind(y)~x+z+w)
  term.name(~s(x)+y) 
  term.name(~x) 
  term.name(y~1) ## this does not work.
  term.name(~1) ## this does not work.
}
term.lhs <- function(formula){
  l <- terms(formula)
  vars_ <- as.character(attr(l,"variable"))[-1]
  resp_ <- attr(l,"response")
  if(resp_>0){
    return(vars_[resp_])
  }
  return(NULL)
}
var.lhs <- function(formula){
  l <- terms(formula)
  vars_ <- as.character(attr(l,"variable"))[-1]
  resp_ <- attr(l,"response")
  if(resp_>0){
    tmp <- reformulate(vars_[resp_])
    return(all.vars(tmp))
  }
  return(NULL)
}
test.lhs <- function(){
  formula <- cbind(x,y)~s(z)
}
var.rhs <- function(formula){
  all.vars(reformulate(term.rhs(formula)))
}
term.rhs <- function(formula,arg=T){
  l <- terms(formula)
  vars_ <- as.character(attr(l,"variable"))[-1]
  resp_ <- attr(l,"response")
  if(resp_>0){
    vars_ <- vars_[-resp_]
  }
  vars_
}
test.rhs <- function(){
  term.rhs(y+z~s(x)) 
  term.rhs(y~s(x))
  term.rhs(~s(x))
  term.rhs(y~x+s(y))
  term.rhs(y~x+y)
  term.rhs(y~1)
  term.rhs(y~0)
}
term.lhs <- function(formula){
  ## extract response from univariate or one
  ## sided formula
  l <- terms(formula)
  vars_ <- as.character(attr(l,"variable"))[-1]
  vars_[attr(l,"response")]
}
test.lhs <- function(){
  all.vars(cbind(n,y)~s(x)) ## need a way to figure out all.var
  term.lhs(cbind(n,y)~s(x)) ## this again does not work
  term.lhs(y+z~s(x)) ## this does not work
  term.lhs(y~s(x))
  term.lhs(~s(x))
}
