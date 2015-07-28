
findrho <- function(x,y,c){
  a=-1; b=1, precision=0.00001
  while( (b-a)>precision){
    m <- (a+b)/2
    if(pmvnorm(x,y,corr=m)==c) break
    if(pmvnorm(x,y,corr=m)<c) a <- m
    else b <- m
  }
  return(m)
}

covmat <- function(nr,cormat,probs,quants,covmat){
  covmat <- matrix(NA,nrow=nrow(cormat),ncol=ncol(cormat))
  for(i in 0:nr){
    for(j in i:nr){
      c <- cormat[i,j]*sqrt(probs[i]*(1-probs[i])*probs[j]*(1-probs[j])) + probs[i]*probs[j]
      covmat[i,j] <- findrho(quants[i],quants[j],c)
    }
  }
  return(covmat)
}