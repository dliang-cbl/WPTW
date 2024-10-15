sim2jam_study <- function (sam, pregam, edf.type = 2) 
{
  if (is.null(sam$b)) 
    stop("coefficient simulation data is missing")
  #pregam$Vp <- cov(t(sam$b[, , 1,drop=FALSE]))
  pregam$Vp <- sam$Vp
  #pregam$coefficients <- rowMeans(sam$b[, , 1,drop=FALSE])
  pregam$coefficients <- sam$b
  #pregam$sig2 <- if (is.null(sam$scale)) 
  #else mean(sam$scale)
  pregam$sig2 <- sam$scale
  if (edf.type > 0) {
    if (is.null(sam$mu)) {
      eta <- pregam$X %*% pregam$coefficients
      mu <- pregam$family$linkinv(eta)
    }
    else {
      #mu <- rowMeans(sam$mu)
      mu <- sam$mu
      eta <- pregam$family$linkfun(mu)
    }
    w <- as.numeric(pregam$w * pregam$family$mu.eta(eta)^2/pregam$family$variance(mu))
    XWX <- t(pregam$X) %*% (w * pregam$X)
  }
  else XWX <- t(pregam$X) %*% (pregam$X)
  if (edf.type < 2) {
    rho <- rowMeans(sam$rho)
    lambda <- exp(rho)
    XWXS <- XWX
    for (i in 1:length(lambda)) {
      ind <- pregam$off[i]:(pregam$off[i] + ncol(pregam$S[[i]]) - 
                              1)
      XWXS[ind, ind] <- XWXS[ind, ind] + pregam$S[[i]] * 
        lambda[i]
    }
    pregam$edf <- diag(solve(XWXS, XWX))
  }
  else pregam$edf <- rowSums(pregam$Vp * t(XWX))/pregam$sig2
  class(pregam) <- "jam"
  pregam
}