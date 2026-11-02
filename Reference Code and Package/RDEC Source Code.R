# RDEC package source code

nrdec.fit <- function( x, y, id, S, omegainit=c(.0,.0), ltol=.01, omega.low = c(0.001,0), omega.high=c(.95,.95), verbose=FALSE
)
{
  if (length(omegainit)!=2) stop("omegainit must have length 2")
  if (length(y) != nrow(x)) stop("x and y non-conforming")
  if (length(y) != length(id)) stop("id and y non-conforming")
  if (length(y) != length(S)) stop("id and s non-conforming")
  
  omega <- omega <- omegainit
  
  beta.sig2 <- gls.beta.sig2( x, y, id, S, omega )
  # what is gls.beta.sig2 function?
  
  beta <- beta.sig2$beta
  sig2 <- beta.sig2$sig2
  newl <- rdec.m2llik( omega, x, y, id, S, beta, sig2 )
  
  if (verbose) cat("Initial -2 ln L = ",newl,"\n")
  oldl <- newl + 10*ltol
  while ( (oldl - newl) > ltol )
  {
    oldl <- newl
    optstr <-        nlminb( start=omega, objective = rdec.m2llik,
                             gradient=rdec.prof.grad, hessian=rdec.prof.hess, lower=omega.low, 
                             upper=omega.high, x=x, y=y, id=id, S=S, beta=beta, sig2=sig2
    )
    
    msgs <- c("X CONVERGENCE", "RELATIVE FUNCTION CONVERGENCE", 
              "BOTH X AND RELATIVE FUNCTION CONVERGENCE", 
              "ABSOLUTE FUNCTION CONVERGENCE")
    if (is.na(match(optstr$message,msgs))) 
      warning(paste("Convergence suspect",optstr$message))
    
    omega <- optstr$parameter
    beta.sig2 <- gls.beta.sig2( x, y, id, S, omega )
    
    beta <- beta.sig2$beta
    sig2 <- beta.sig2$sig2
    newl <- rdec.m2llik( omega, x, y, id, S, beta, sig2 )
    
    if (verbose) cat("Updated -2 ln L = ",newl,"\n")
  }
  
  optstr$aux <- NULL
  list( beta=beta, sig2=sig2, omega=omega, m2ll=newl, opt=optstr )
  
}
rdec.m2llik <- function( omega, x, y, id, S, beta, sig2 )
{
  N <- length(y)
  z <- (y - x %*% beta)/sqrt(sig2)
  (N*log(2*pi*sig2) + log.det.bd.dec( id, S, omega )
    + quadform.bd.dec( z, id, S, omega )
  )
}
quadform.bd.dec <- function(x, id, S, omega)
{
  n <- nrow(x)
  p <- ncol(x)
  if(length(omega) != 2)
    stop("omega must have 2 elements")
  qfp <- double(p * p)
  
  ans <-         ans <- .C("eval_quadform_bdmat_dec_intf",
                           x,
                           as.integer(n),
                           as.integer(p),
                           id,
                           S,
                           omega,
                           qfp = qfp, PACKAGE="RDEC")
  matrix(ans$qfp,p,p)
}

gls.beta.sig2 <-
  function(x, y, id, S, omega)
  {
    xpvix <- quadform.bd.dec(x, id, S, omega)
    xpviy <- wtd.inrprod.bd.dec(x, y, id, S, omega)
    beta <- solve(xpvix) %*% xpviy
    e <- y - x %*% beta
    sig2 <- quadform.bd.dec(e, id, S, omega)/length(y)
    list(beta = beta, sig2 = as.vector(sig2))
  }

log.det.bd.dec <- function( id, S, omega )
{
  N <- length(id)
  if (length(omega) != 2)stop("omega must be length 2")
  ans <- .C("eval_det_bdmat_dec_intf",
            as.integer(N),
            id,
            S,
            omega,
            ldetout=double(1), PACKAGE="RDEC")
  ans$ldetout
}


log.det.bd.cs <- function( id, rho )
{
  N <- length(id)
  ans <- .C("eval_det_bdmat_cs_intf",
            as.integer(N),
            id,
            rho,
            ldetout=double(1), PACKAGE="RDEC")
  ans$ldetout
}


wtd.inrprod.bd.dec <-
  function(x, y, id, S, omega)
  {
    if(length(omega) != 2) stop("omega must be length 2")
    if(is.matrix(x)) {
      n <- nrow(x)
      p <- ncol(x)
    }
    else {
      n <- length(x)
      p <- 1
    }
    if(is.matrix(y)) {
      q <- ncol(y)
    }
    else q <- 1
    ipp <- double(p * q)
    ans <- .C("eval_wtd_inrprod_dec_intf",
              x,
              as.integer(n),
              as.integer(p),
              y,
              as.integer(q),
              id,
              omega,
              S,
              ipp = ipp, PACKAGE="RDEC")
    matrix(ans$ipp, p, q)
  }


rdec.omega.score <- function( omega, x, y, id, S, beta, sig2 )
{
  n <- length(y)
  p <- ncol(x)	
  # void score_mvg_dec_intf( double* omega, double* beta, double* sig2, double* x, 
  #  double* y, double* id, double* s, int* n, int* p, double* scomega );
  ans <- .C("score_mvg_dec_intf",
            omega,
            beta,
            sig2,
            x,
            y,
            id,
            S,
            as.integer(n),
            as.integer(p),
            sc.omega = double(2), PACKAGE="RDEC")
  ans$sc.omega
}

rdec.prof.grad <- function( omega, x, y, id, S, beta, sig2 )
  -2*rdec.omega.score( omega, x, y, id, S, beta, sig2 )

rdec.full.hess <- function( omega, x, y, id, S, beta, sig2 )
{
  #void hess_mvg_dec_intf( double* omega, double* beta, double* sig2, double* x, 
  #                                 double* y, double* id, double* s, int* n, 
  #                                 int* p, double* hessbeta_g, double* hessbeta_t, 
  #                                 double* hesssig2, double* hsg, double* hst, 
  #					double* hessomega )
  p <- ncol(x)
  n <- nrow(x)
  hbb <- -1/sig2 * quadform.bd.dec(x, id, S, omega)
  hbg <- double(p)
  hbt <- double(p)
  hss <- double(1)
  hsg <- double(1)
  hst <- double(1)
  homega <- double(4)
  ans <- .C("hess_mvg_dec_intf",
            as.double(omega),
            as.double(beta),
            as.double(sig2),
            as.double(x),
            as.double(y),
            as.double(id),
            as.double(S),
            as.integer(n),
            as.integer(p),
            Hbg = hbg,
            Hbt = hbt,
            Hss = hss,
            Hsg = hsg,
            Hst = hst,
            Homega = homega, PACKAGE="RDEC")
  ans$Homega <- matrix(ans$Homega, 2, 2)
  tmp1 <- cbind(hbb, 0, ans$Hbg, ans$Hbt)
  tmp2 <- rbind(tmp1, c(rep(0, p), ans$Hss, ans$Hsg, ans$Hst))
  tmp3 <- rbind(tmp2, c(ans$Hbg, ans$Hsg, ans$Homega[1,  ]))
  tmp4 <- rbind(tmp3, c(ans$Hbt, ans$Hst, ans$Homega[2,  ]))
  list(hbb = hbb, homega = ans$Homega, mat = tmp4)
}

rdec.prof.hess <- function( omega, x, y, id, S, beta, sig2 )
{
  ans <- rdec.full.hess( omega, x, y, id, S, beta, sig2 )
  as.vector(-2*ans$homega)[c(1,2,4)]
}


rdec <- function(formula, id, S, data = sys.parent(), subset, na.action=na.fail, omega.init = c(0, 0), 
                 omega.low = c(0.001, 0), omega.high = c(0.95, 0.95), ltol = 0.01, 
                 contrasts = NULL, verbose=FALSE)
{
  call <- match.call()
  m <- match.call(expand = FALSE)
  if(missing(data) & !missing(na.action))
    stop("supply data frame (data=) if na.action is to be used")
  if(!missing(data))
    attach(na.action(data))
  m$id <- m$S <- m$omega.init <- m$omega.low <- m$omega.high <- m$ltol <- 
    m$contrasts <- m$verbose <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  X <- model.matrix(Terms, m, "contrasts")
  id <- as.double(id)
  S <- as.double(S)
  fit <- rdec.fit(X, Y, id, S, omegainit = omega.init, omega.high = 
                    omega.high, omega.low = omega.low, ltol = ltol, verbose=verbose)
  fit$terms <- Terms
  fit$call <- call
  fit$N <- length(Y)
  fit$Nclust <- length(unique(id))
  fit$fitted.values <- X %*% fit$coefficients
  fit$residuals <- Y - fit$fitted.values
  attr(fit, "na.message") <- attr(m, "na.message")
  class(fit) <- c("rdec") #, "lm")
  fit
}


rdec.fit <- function (x, y, id, S, omegainit = c(0, 0), ltol = 0.01, omega.low = c(0.001, 
                                                                                   0), omega.high = c(0.95, 0.95), verbose=FALSE) 
{
  if (length(omegainit) != 2) 
    stop("omegainit must have length 2")
  if (length(y) != nrow(x)) 
    stop("x and y non-conforming")
  if (length(y) != length(id)) 
    stop("id and y non-conforming")
  if (length(y) != length(S)) 
    stop("id and s non-conforming")
  omega <- omega <- omegainit
  beta.sig2 <- gls.beta.sig2(x, y, id, S, omega)
  beta <- beta.sig2$beta
  sig2 <- beta.sig2$sig2
  newl <- rdec.m2llik(omega, x, y, id, S, beta, sig2)
  if (verbose) cat("Initial -2 ln L = ", newl, "\n")
  oldl <- newl + 10 * ltol
  while ((oldl - newl) > ltol) {
    oldl <- newl
    optstr <- optim(par = omega, fn = rdec.m2llik, method = "L-BFGS-B", 
                    gr = rdec.prof.grad, lower = omega.low, upper = omega.high, 
                    x = x, y = y, id = id, S = S, beta = beta, sig2 = sig2)
    msgs <- c("X CONVERGENCE", "RELATIVE FUNCTION CONVERGENCE", 
              "BOTH X AND RELATIVE FUNCTION CONVERGENCE", "ABSOLUTE FUNCTION CONVERGENCE")
    if (is.na(match(optstr$message, msgs))) 
      warning(paste("Convergence suspect", optstr$message))
    omega <- optstr$par
    beta.sig2 <- gls.beta.sig2(x, y, id, S, omega)
    beta <- beta.sig2$beta
    sig2 <- beta.sig2$sig2
    newl <- rdec.m2llik(omega, x, y, id, S, beta, sig2)
    if (verbose) cat("Updated -2 ln L = ", newl, "\n")
  }
  optstr$aux <- NULL
  hess <- rdec.full.hess(omega, x, y, id, S, beta, sig2)$mat
  xnames <- dimnames(x)[[2]]
  beta <- as.vector(beta)
  names(beta) <- xnames
  allnames <- c(xnames, "sigma2", "gamma", "theta")
  dimnames(hess) <- list(allnames, allnames)
  names(omega) <- c("gamma", "theta")
  theta <- omega["theta"]
  hessd <- nrow(hess)
  if (theta < 0.001 | (theta > 0.999 & theta < 1.001)) 
    hess <- hess[-hessd, -hessd]
  list(coefficients = beta, sig2 = sig2, omega = omega, m2ll = newl, 
       opt = optstr, rank = ncol(x), hessian = hess)
}

old.rdec.fit <- function(x, y, id, S, omegainit = c(0, 0), ltol = 0.01, omega.low = c(0.001, 0),
                         omega.high = c(0.95, 0.95), verbose=FALSE)
{
  if(length(omegainit) != 2)
    stop("omegainit must have length 2")
  if(length(y) != nrow(x))
    stop("x and y non-conforming")
  if(length(y) != length(id))
    stop("id and y non-conforming")
  if(length(y) != length(S))
    stop("id and s non-conforming")
  omega <- omega <- omegainit
  beta.sig2 <- gls.beta.sig2(x, y, id, S, omega)
  beta <- beta.sig2$beta
  sig2 <- beta.sig2$sig2
  newl <- rdec.m2llik(omega, x, y, id, S, beta, sig2)
  if (verbose) cat("Initial -2 ln L = ", newl, "\n")
  oldl <- newl + 10 * ltol
  while((oldl - newl) > ltol) {
    oldl <- newl
    optstr <- nlminb(start = omega, objective = rdec.m2llik, 
                     gradient = rdec.prof.grad, hessian = rdec.prof.hess, 
                     lower = omega.low, upper = omega.high, x = x, y = y, id
                     = id, S = S, beta = beta, sig2 = sig2)
    msgs <- c("X CONVERGENCE", "RELATIVE FUNCTION CONVERGENCE", 
              "BOTH X AND RELATIVE FUNCTION CONVERGENCE", 
              "ABSOLUTE FUNCTION CONVERGENCE")
    if(is.na(match(optstr$message, msgs)))
      warning(paste("Convergence suspect", optstr$message))
    omega <- optstr$parameter
    beta.sig2 <- gls.beta.sig2(x, y, id, S, omega)
    beta <- beta.sig2$beta
    sig2 <- beta.sig2$sig2
    newl <- rdec.m2llik(omega, x, y, id, S, beta, sig2)
    if (verbose) cat("Updated -2 ln L = ", newl, "\n")
  }
  optstr$aux <- NULL
  hess <- rdec.full.hess(omega, x, y, id, S, beta, sig2)$mat
  xnames <- dimnames(x)[[2]]
  beta <- as.vector(beta)
  names(beta) <- xnames
  allnames <- c(xnames, "sigma2", "gamma", "theta")
  dimnames(hess) <- list(allnames, allnames)
  names(omega) <- c("gamma", "theta")
  theta <- omega["theta"]
  hessd <- nrow(hess)
  if (theta < .001 | (theta > .999 & theta < 1.001))
    hess <- hess[-hessd,-hessd]
  list(coefficients = beta, sig2 = sig2, omega = omega, m2ll = newl, opt
       = optstr, rank = ncol(x), hessian = hess)
}

print.rdec <- function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  coef <- coefficients(x)
  cat("\nCoefficients:\n")
  print(coef, ...)
  rank <- x$rank
  if(is.null(rank))
    rank <- sum(!nas)
  nobs <- x$N
  rdf <- nobs - x$rank
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
  if(!is.null(attr(x, "na.message")))
    cat(attr(x, "na.message"), "\n")
  invisible(x)
}

summary.rdec <- function(object, wt = NULL)
{
  #this method is designed on the assumption that the coef method
  # returns only the estimated coefficients.  It will (it's asserted)
  # also work, however, with fitting methods that don't follow this
  # style, but instead put NA's into the unestimated coefficients
  coef <- coefficients(object)
  cnames <- labels(coef)
  ctotal <- object$coef
  ptotal <- length(ctotal)
  resid <- object$residuals
  fv <- fitted(object)
  n <- length(resid)
  p <- object$rank
  if(is.null(p))
    p <- sum(!is.na(ctotal))
  if(any(na <- is.na(coef))) {
    coef <- coef[!na]
    p <- length(coef)
  }
  rdf <- object$df.resid
  if(is.null(rdf))
    rdf <- n - p
  if(!is.null(wt)) {
    wt <- wt^0.5
    resid <- resid * wt
    fv <- fv * wt
    excl <- wt == 0
    if(any(excl)) {
      warning(paste(sum(excl), 
                    "rows with zero weights not counted"))
      resid <- resid[!excl]
      fv <- fv[!excl]
      if(is.null(object$df.residual))
        rdf <- rdf - sum(excl)
    }
  }
  #	rinv <- diag(p)
  #	dimnames(rinv) <- list(cnames, cnames)
  stddev <- sqrt(object$sig2)
  #	R <- object$R
  #	if(p < ptotal)
  #		R <- R[1:p, 1:p, drop = F]
  #	rinv <- solve(R, rinv)
  #	rowlen <- (rinv^2 %*% rep.int(1, p))^0.5
  #	names(rowlen) <- cnames
  #	if(correlation) {
  #		correl <- rinv * array(1/rowlen, c(p, p))
  #		correl <- correl %*% t(correl)
  #	}
  #	else correl <- NULL
  coef <- array(coef, c(p, 4))
  dimnames(coef) <- list(cnames, c("Value", "Std. Error", "Z value", 
                                   "Pr(>|Z|)"))
  InvInfo <- solve( - object$hessian)
  coef[, 2] <- sqrt(diag(InvInfo))[1:p]
  coef[, 3] <- coef[, 1]/coef[, 2]
  coef[, 4] <- if(rdf > 0) 2 * (1 - pnorm(abs(coef[, 3]))) else NA
  yy <- fv + resid
  int <- attr(object$terms, "intercept")
  if(is.null(int)) {
    r2 <- NA
    fstat <- NA
  }
  else if(int) {
    if(is.null(wt)) {
      # intercept and no weights
      mn <- mean(fv)	# same as mean(yy)
      r2 <- sum((fv - mn)^2)/sum((yy - mn)^2)
      fstat <- c(value = (sum((fv - mn) * (yy - mn))/(p - 1) * 
                            rdf)/sum(resid^2), numdf = p - 1, dendf = rdf)
    }
    else {
      # intercept and weights
      w <- wt/sum(wt^2)
      r2 <- sum((fv - sum(fv * wt) * w)^2)/sum((yy - sum(yy * 
                                                           wt) * w)^2)
      fstat <- c(value = (sum((fv - sum(fv * wt) * w) * (yy - 
                                                           sum(yy * wt) * w))/(p - 1) * rdf)/sum(resid^2), 
                 numdf = p - 1, dendf = rdf)
    }
  }
  else {
    # no intercept
    r2 <- sum((fv)^2)/sum((yy)^2)
    fstat <- c(value = (sum(fv * yy)/p * rdf)/sum(resid^2), numdf
               = p, dendf = rdf)
  }
  initobj <- object
  object <- object[c("call", "terms")]
  object$residuals <- resid
  object$coefficients <- coef
  object$sigma <- stddev # sqrt MLE
  object$omega <- initobj$omega
  theta <- object$omega[2]
  if ( theta < .001 | ( theta > .999 & theta < 1.001 )) object$se.omega <- c(sqrt(InvInfo[p + (2), p + (2)]),NA)
  else object$se.omega <- sqrt(diag(InvInfo[p + (2:3), p + (2:3)]))
  object$df <- c(p, rdf, ptotal)
  object$cov.unscaled <- InvInfo[1:p, 1:p]	#	object$correlation <- correl
  object$m2ll <- initobj$m2ll
  class(object) <- "summary.rdec"
  object
}

print.summary.rdec <- function(x, digits = max(3, .Options$digits - 3), ...)
{
  cat("\nCall: ")
  dput(x$call)
  resid <- as.vector(x$residuals)
  df <- x$df
  rdf <- df[2]
  if(rdf > 5) {
    cat("Residuals:\n")
    if(length(dim(resid)) == 2) {
      rq <- apply(t(resid), 1, quantile)
      dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", 
                             "Max"), dimnames(resid)[[2]])
    }
    else {
      rq <- quantile(resid)
      names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print(rq, digits = digits, ...)
  }
  else if(rdf > 0) {
    cat("Residuals:\n")
    print(resid, digits = digits, ...)
  }
  cat("\nCorrelation parameter estimates:\n")
  corm <- rbind(x$omega, x$se.omega)
  dimnames(corm)[[1]] <- c("par. est.", "s.e.")
  print(corm)
  cat("\nRegression coefficients:\n")
  print(format(round(x$coef, digits = digits)), quote = FALSE, ...)
  cat("\nMLE of resid. variance:", format(signif(x$sigma^2, digits)), "\n")
  cat("-2 log likelihood: ", format(x$m2ll, digits = log10(abs(x$m2ll)) + 4), 
      "\n")
  invisible(x)
}


final_data_beforeretire$CONST<-rep(1,nrow(final_data_beforeretire))
model.dec<-rdec(YREXP~CONST, id=ID, S=YR_TV, omega.init = c(.5,.5), omega.low=c(.01,.01),
              omega.high=c(.95, .95), data=final_data_beforeretire)

valid_data<-valid_data %>% group_by(ID) %>%
  mutate(timing=YEAR-YEAR[1]+1)
error.model.dec<-rdec(EVALUE~VALUE, id=ID, S=timing, omega.init = c(.5,.5), omega.low=c(.01,.01),
                      omega.high=c(.95, .95), data=valid_data)



