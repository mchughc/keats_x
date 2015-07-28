
library(MASS)
library(GWASTools)
library(QCpipeline)
#library(OLGApipeline)
#library(OLGAanalysis)


setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
#setwd("~/Documents/keats-x")

source("nullMixedModel_v3.R")
source("createDesignMatrix.R")
source("estVarComp_combined_v4.R")

install.packages("kinship_1.1.0-23.tar.gz", repos = NULL, type = "source")
library(kinship)
wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
getSNP <- function(){
  set.seed(rpois(1,lambda=100))
  return(sample(0:2,1000,replace=TRUE,prob=c(0.81,0.18,0.01)))
}
         
#set.seed(12345)
fullped<-data.frame(famid=rep(1:250,each=4),id=10001:11000,fa=rep(0,1000),mo=rep(0,1000))
fullped$fa[(1:250)*4-1]<-fullped$fa[(1:250)*4]<-(1:250)*4+9997
fullped$mo[(1:250)*4-1]<-fullped$mo[(1:250)*4]<-(1:250)*4+9998

fullkins=makekinship(fullped$famid, fullped$id, fullped$fa, fullped$mo)
sqrtweights<- wuweights

id <- 10001:11000
tmpidx<-!is.na(match(dimnames(fullkins)[[1]], id))
tmpkins<-fullkins[tmpidx, tmpidx]
tmpidx<-match(id, dimnames(tmpkins)[[1]])
if(any(is.na(tmpidx))) stop("Some id not exist in fullkins. Check your full pedigree data...")
kins<-as.matrix(tmpkins)[tmpidx, tmpidx]

snpStart <- 1
snpEnd <- 1:1000

kinAuto <- as.matrix(fullkins)
colnames(kinAuto) <- id
covMatList <- list(kinAuto*2) # this is the kinship matrix for all samples
names(covMatList) <- c("kinshipAutos")

res <- data.frame(matrix(NA,nrow=1000,ncol=3))
colnames(res) <- c("mean","max","min")


for(j in 1:1000){

  geno <- data.frame(matrix(NA,nrow=1000,ncol=501))
  colnames(geno) <- c("id",paste0("snp",2:500))

  geno$id <- 10001:11000
  for(i in 2:ncol(geno)){
    geno[,i] <- getSNP()
  }

  set.seed(j)
  noise <- rnorm(1000,0,1)
  pheno<-data.frame(id=10001:11000,y=noise+mvrnorm(1,mu=rep(0,1000),Sigma=0.5*kinAuto),sex=rep(1:2,500),bmi=rnorm(1000,25,2))

#  sigmaMat <- mvrnorm(1000,mu=rep(0,8),Sigma=s2X*kinX[1:8,1:8])

  phenotype=pheno$y
  id=pheno$id
  covariates=as.matrix(pheno[,c("sex","bmi")])
  n<-length(phenotype)

  X<-cbind(rep(1,n),as.matrix(covariates))
  tmpdata<-data.frame(phenotype,id,covariates)
  nc<-ncol(covariates)
  exprs<-paste("phenotype ~", paste(names(tmpdata)[-c(1,2)],collapse=" + "))

  sex=tmpdata$sex
  sex[sex==1] <- "M"
  sex[sex==2] <- "F"
  scan1 <- data.frame(scanID=tmpdata$id,sex=sex,pheno=tmpdata$phenotype,bmi=tmpdata$bmi)
  scan1 <- ScanAnnotationDataFrame(scan1)


#geno<-data.frame(id=10001:11000,
#                 snp1=sample(0:2,1000,replace=T,prob=c(0.81,0.18,0.01)),
#                 snp2=sample(0:2,1000,replace=T,prob=c(0.81,0.18,0.01)),
#                 snp3=sample(0:2,1000,replace=T,prob=c(0.81,0.18,0.01)),
#                 snp4=getSNP(),
#                 snp5=getSNP(),snp6=getSNP()
#                 )


  genotypes=as.matrix(geno[,-1])
  
  Z<-genotypes
  nZ<-ncol(Z)
  MAF<-colMeans(Z, na.rm = TRUE)/2

  weights<-sqrtweights(MAF)
  G<-t(t(Z)*weights)

  nullmod<-lmekin(as.formula(exprs),tmpdata,random=~1|id,varlist=list(kinAuto))
  
  
  genoMt <- MatrixGenotypeReader(genotype=t(genotypes),snpID=as.integer(1:ncol(genotypes)),
                                 chromosome=as.integer(rep(23,ncol(genotypes))),
                                 position=as.integer(1:ncol(genotypes)), scanID=scan1$scanID)
  
  genoData <- GenotypeData(genoMt,scanAnnot=scan1)

  #varCompBoth <- estVarComp(scan1,covMatList=covMatList,"pheno",dropZeros=FALSE,covar.vec=c("sex","bmi"))#,AIREML.tol=1e-10)
varCompBoth <- estVarComp(scan1,covMatList=covMatList,"pheno",covar.vec=c("sex","bmi"))#,AIREML.tol=1e-10)
  varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
  estVarCompCI(varCompBoth,prop=TRUE)


  
  SIGMA<-nullmod$theta[1]*2*kinAuto+nullmod$theta[2]*diag(n)
  SIGMAi <- solve(SIGMA)
#  diffV <- SIGMAi-varCompBoth$cholSigmaInv
  
  resids <- nullmod$res
  Q<-t(resids) %*% SIGMAi %*% G %*% t(G) %*% SIGMAi %*% resids
  P<-SIGMA - X %*% solve(t(X) %*% SIGMAi %*% X) %*% t(X)
  eig<-eigen(t(G) %*% SIGMAi %*% P %*% SIGMAi %*% G, symmetric=T, only.values=T)
  evals<-eig$values[eig$values>1e-6*eig$values[1]]

  ###
  
   ## run nullMixedModel so we have the output somewhere
  mod <- nullMixedModel(scanAnnot = scan1,
                        VC = varCompBoth,
                        outcome = "pheno",
                        covar.vec = colnames(covariates),
                        verbose=TRUE)
  
  resids.Matt <- mod$resid.marginal

  ourSig <- varCompBoth$varComp[1]*2*kinAuto + varCompBoth$varComp[2]*diag(n)
  ourSigI <- solve(ourSig) # differs from the varCompBoth$cholSigmaInv
  ourSigI <- varCompBoth$cholSigmaInv %*% varCompBoth$cholSigmaInv
  cholDec <- chol(ourSig)
  
  QMine <- t(resids.Matt) %*% ourSigI %*% G %*% t(G) %*% ourSigI %*% resids.Matt
  P <- ourSig - X %*% solve(t(X) %*% ourSigI %*% X) %*% t(X)

  tmpEs <- eigen(t(G) %*% ourSigI %*% P %*% ourSigI %*% G, symmetric=TRUE, only.values=TRUE)
  ourEs <-   tmpEs$values[tmpEs$values>1e-6*tmpEs$values[1]]

pfamskat<-pchisqsum(Q,rep(1,length(evals)),evals,lower=F,method="sad")
pfamskat<-min(pfamskat,1)
pfamskat<-max(pfamskat,0)
res <- data.frame(pvalue=pfamskat,method="famSKAT")

pfamskat<-pchisqsum(Q,rep(1,length(ourEs)),ourEs,lower=F,method="sad")
pfamskat<-min(pfamskat,1)
pfamskat<-max(pfamskat,0)

tmp <- data.frame(pvalue=pfamskat,method="keatsX")
res <- rbind(res,tmp)

write.table(res,file="/projects/geneva/geneva_sata/caitlin/keats_x/tmp_testRes.txt",quote=FALSE,row.names=FALSE)

  print(j)
}

write.table(res,file="/projects/geneva/geneva_sata/caitlin/keats_x/testEvals.txt",quote=FALSE,row.names=FALSE)

q("no")
