nullMixedModel <- function(scanAnnot,
                           VC,
                           outcome,
                           covar.vec = NULL,
                           scan.include = NULL,
                           group.var = NULL,
                           verbose = TRUE){

  cholSigmaInv <- VC$cholSigmaInv
  varComp <- VC$varComp

  # read in outcome and covariate data
  if(verbose) message("Reading in Phenotype and Covariate Data...")
  dat <- .createDesignMatrix(scanAnnot = scanAnnot, outcome = outcome, covar.vec = covar.vec, ivar.vec = NULL, scan.include = scan.include)
  # outcome vector
  Y <- dat$Y
  # design matrix
  W <- dat$W
  k <- ncol(W)
  rm(dat)
  # scans to include
  scan.include <- rownames(W)

  # which samples to remove from cholSigmaInv
  if(!all(scan.include %in% colnames(cholSigmaInv))){
    stop("All of the included Samples must be in the cholSigmaInv matrix")
  }
  chol.idx <- which(!(colnames(cholSigmaInv) %in% scan.include))
  cholSigmaInv <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)
  

  # sample size 
  n <- length(scan.include)
  if(verbose) message("Fitting Model with ", n, " Samples")

  # calculate fixed effects estimates (beta)
  CW <- crossprod(cholSigmaInv,W)
  CY <- crossprod(cholSigmaInv,Y)
  XtSigInvX <- crossprod(CW)
  XtSigInvXInv <- chol2inv(chol(XtSigInvX))
  beta <- crossprod(XtSigInvXInv, crossprod(CW,CY))


  # marginal residuals
  fits <- tcrossprod(W, t(beta))
  residM <- as.vector(Y - fits)

  # conditional residuals  
  m <- length(varComp)
  residtmp <- crossprod(residM,cholSigmaInv)

  if(VC$hetResid){
    if(is.null(group.var)){
      stop("group.var must be specified when heterogeneous residual variances were used to estimate variance components")
    }
    # vector of estimated residual variances
    sigma2epsilon <- rep(NA, n)
    # read in group variable and subset to included samples
    group <- getVariable(scanAnnot, group.var)[getScanID(scanAnnot) %in% scan.include]
    # determine unique groups
    group.names <- as.character(unique(group))
    g <- length(group.names)
    for(i in 1:g){
      sigma2epsilon[group == group.names[i]] <- varComp[m-g+i]
    }
    residC <- as.vector(sigma2epsilon*tcrossprod(cholSigmaInv, residtmp))

  }else{
    residC <- as.vector(varComp[m]*tcrossprod(cholSigmaInv, residtmp))
  }

  # Variance Covariance of betas
  RSS <- sum(residtmp^2)/(n-k)
  Vbeta <- RSS*XtSigInvXInv

  # test statistics and p-values
  SE <- sqrt(diag(Vbeta))
  Stat <- (beta/SE)^2
  pval <- pchisq(Stat, df=1, lower.tail=FALSE)

  # find likelihood and AIC
  logLik <- as.numeric(-0.5*n*log(2*pi*RSS) + sum(log(diag(cholSigmaInv))) - 0.5*tcrossprod(residtmp)/RSS)
  AIC <- 2*(k+m) - 2*logLik
  logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - 0.5*log(det(XtSigInvX)) )  

  # prepare results
  dimnames(Vbeta) <- list(colnames(W), colnames(W))

  fixef <- data.frame(Est = beta, SE = SE, Stat = Stat, pval = pval)
  rownames(fixef) <- colnames(W)
  
  return(list(fixef = fixef,
              varComp = varComp,
              resid.marginal = residM,
              resid.conditional = residC,
              logLik = logLik,
              AIC = AIC,
              logLikR = logLikR,
              model.matrix = W,
              Vbeta = Vbeta,
              RSS =RSS))
}



.subsetCholSigmaInv <- function(cholSigmaInv, chol.idx) {
  if(length(chol.idx) > 0){
    # subset cholSigmaInv
    SigmaInv <- tcrossprod(cholSigmaInv)
    for(i in sort(chol.idx, decreasing=TRUE)){
      SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
    }
    cholSigmaInv <- t(chol(SigmaInv))
  }
  
  cholSigmaInv
}
