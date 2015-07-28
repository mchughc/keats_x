.createDesignMatrix <- function(scanAnnot, outcome, covar.vec, ivar.vec, scan.include){
  # set which samples to keep
  scanID <- getScanID(scanAnnot)
  # samples to include for the analysis
  if(!is.null(scan.include)){
    if(!all(scan.include %in% scanID)){
      stop("Not all of the scan IDs in scan.include are in scanAnnot")
    }
    keep <- scanID %in% scan.include
  }else{
    keep <- rep(TRUE, length(scanID))
  }  

  if(!is.null(covar.vec)){
    cvnames <- unique(unlist(strsplit(covar.vec,"[*:]")))
    # read in data
    dat <- as.data.frame(getVariable(scanAnnot, c(outcome,cvnames)))
    # identify samples with any missing data
    keep <- keep & apply(dat,1,function(x){ all(!is.na(x)) })
    
    # read in interaction variable data
    if(!is.null(ivar.vec)){
      ivnames <- unique(unlist(strsplit(ivar.vec,"[*:]")))
      # read in data
      idat <- as.data.frame(getVariable(scanAnnot, c(outcome,ivnames)))
      # identify samples with any missing data
      keep <- keep & apply(idat,1,function(x){ all(!is.na(x)) })
      # remove samples with any missing data
      idat <- as.data.frame(idat[keep,])
      # create design matrix for interaction variables
      iformula <- as.formula(paste(paste(outcome,"~"), paste(ivar.vec,collapse="+")))
      V <- model.matrix(iformula, data=idat)
      rownames(V) <- scanID[keep]
    }    
    
    # remove samples with any missing data
    dat <- as.data.frame(dat[keep,])
    # outcome vector
    Y <- dat[,outcome]
    # create design matrix
    model.formula <- as.formula(paste(paste(outcome,"~"), paste(covar.vec,collapse="+")))
    W <- model.matrix(model.formula, data=dat)
    rownames(W) <- scanID[keep]
    
  }else{
    # read in data
    dat <- getVariable(scanAnnot,outcome)
    # identify samples with any missing data
    keep <- keep & !is.na(dat)
    
    # read in interaction variable data
    if(!is.null(ivar.vec)){
      ivnames <- unique(unlist(strsplit(ivar.vec,"[*:]")))
      # read in data
      idat <- as.data.frame(getVariable(scanAnnot, c(outcome,ivnames)))
      # identify samples with any missing data
      keep <- keep & apply(idat,1,function(x){ all(!is.na(x)) })
      # remove samples with any missing data
      idat <- as.data.frame(idat[keep,])
      # create matrix to store interaction variables
      iformula <- as.formula(paste(paste(outcome,"~"), paste(ivar.vec,collapse="+")))
      V <- model.matrix(iformula, data=idat)
      rownames(V) <- scanID[keep]
    }
    
    # outcome
    Y <- dat[keep]
    # design matrix
    W <- matrix(1,nrow=length(Y),ncol=1)
    rownames(W) <- scanID[keep]
  }

  # return output
  if(!is.null(ivar.vec)){
    return(list(Y = Y, W = W, V = V))
  }else{
    return(list(Y = Y, W = W))
  }

}