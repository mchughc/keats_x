# script to run MONSTER on data, compare with keats_optimal output, fit with just autosomal kinship


library(MASS)
library(GWASTools)
library(SNPRelate)
library(OLGAanalysis)
library(survey)
library(kinship)
library(pbivnorm)
library(SKAT)
setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
source("../mlm_x/allele_drop_functions.R")
source("hapsim_functions.R")
source("keats.R")
source("keats_optimal.R")
source("/projects/geneva/geneva_sata/caitlin/keats_x/estVarComp_combined_v4.R")
wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)


nvar <- 20
maf <- runif(nvar,min=0.005,max=0.05)

ld <- getobj("LDMat_ceu_chr9_100snps.RData")
ld <- ld[1:nvar,1:nvar]
ld[lower.tri(ld,diag=FALSE)] <- t(ld)[lower.tri(ld,diag=FALSE)]
ld <- makepd(ld)

#effVar <- c(0.5,0,0.5)
#effVar <- c(0.2,0.2,0.4)
effVar <- c(0.1,0,0.9)
n <- 2000
nvar <- as.integer(nvar)

SEX <- c("F","M","F","F","M","M","M","F")

kinAuto <- getobj("../mlm_x/1000Peds_8ped_fem_autoKinship.RData")
sex <- rep(SEX,(n/8))
sex.int <- sex
sex.int[sex.int=="M"] <- 1
sex.int[sex.int=="F"] <- 2
kinAuto <- matrix(kinAuto,nrow=8000,ncol=8000)
kinAuto <- kinAuto[1:n,1:n]

colnames(kinAuto) <- 1:n


# simulate genotypes now
genoMat <- matrix(NA,nrow=nvar,ncol=n) # snp x sample; should be sample x snp
colsToGeno <- split(1:n,rep(1:n,each=8,length=n))
geno <- do.call(cbind,lapply(colsToGeno,function(x){genoMat[,x]=Family_alleles_Nmarker_8ped_fem(maf,nvar,sex,ld)}))
MAF<-rowMeans(geno, na.rm = TRUE)/2 # this needs to change for X chr SNPs

s2A <- 0.5 
sigmaAMat <- mvrnorm((n/8),mu=rep(0,8),Sigma=s2A*kinAuto[1:8,1:8])

noise <- rnorm(n,mean=0,sd=1)

vars <- 1:nvar
notCausal <- sample(vars,nvar*as.numeric(effVar[3]))
negVar <- sample(vars[-notCausal],nvar*as.numeric(effVar[1]))
posVar <- sample(vars[-c(notCausal,negVar)],nvar*as.numeric(effVar[2]))

c <- as.numeric(0.2)
beta <- rep(0,nvar)
beta[posVar] <- c*abs(log10(MAF[posVar]))
beta[negVar] <- -c*abs(log10(MAF[negVar]))
beta[MAF==0] <- 0
bgeno <- beta*geno
y <- colSums(bgeno)+noise
pheno <- data.frame(scanID=1:n,y=y,sex=sex)

scan1 <- ScanAnnotationDataFrame(pheno)

genoBurd <- geno*wuweights(MAF)^2
burdRegs<- matrix(apply(genoBurd,2,sum),nrow=1,ncol=n)

genoMt <- MatrixGenotypeReader(genotype=burdRegs,snpID=as.integer(1:nrow(burdRegs)),
                               chromosome=as.integer(rep(1,nrow(burdRegs))),
                               position=as.integer(1:nrow(burdRegs)), scanID=scan1$scanID)

genoData <- GenotypeData(genoMt,scanAnnot=scan1)

#kinAuto <- 2*kinAuto
#diag(kinAuto) <- 0.5 # needs to be 1+inbd=1 or inbd=0; MONSTER uses inbd
kAuto <- list(kinshipAutos=kinAuto)

vcAuto <- estVarComp(scan1,covMatList=kAuto,"y")
cholSigAuto <- vcAuto[["cholSigmaInv"]]

# fit with 2000 individs, 8ped fem structure
dat <- read.table("monster/input_8ped_fem.txt",header=FALSE)
ped <- data.frame(matrix(NA,nrow=2000,ncol=5))
for(i in 1:250){
  ped[(i*8-7):(i*8),] <- dat
  ped[(i*8-7):(i*8),1] <- i
}

ped$phenotype <- y
ped[,2] <- 1:nrow(ped)
write.table(ped,"monster/pheno.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

# geno data file, snp x sample
mgeno <- cbind(1:nrow(geno),geno)
colnames(mgeno) <- c(0,1:nrow(ped))
write.table(mgeno,"monster/geno.txt",col.names=TRUE,quote=FALSE,sep="\t")

# need to make snp list file
snpList <- c("snpset1",1,1:20)
snpList <- rbind(snpList,c("snpset1",0,wuweights(MAF)^2))
write.table(snpList,"monster/SNP.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

# need to make the kinship matrix; try 2*kinship -- doesn't matter if kin or 2*kin
sm <- 0.5*kinAuto[1:8,1:8]
mkin <- read.table("monster/outfile_8ped_fem.txt")
mkin$V4 <- as.vector(sm[lower.tri(sm,diag=TRUE)])
# need to copy this 250x and change family ids, individ ids
mkin[mkin$V2==mkin$V3,"V4"] <- 0
fullKin <- data.frame(matrix(NA,nrow=nrow(mkin)*250,ncol=ncol(mkin)))
for(i in 1:250){
  fullKin[(i*36-35):(i*36),] <- mkin
  fullKin[(i*36-35):(i*36),1] <- i
  # need to map individ ids to id+(8*(i-1)) 
  v2 <- fullKin$X2[(i*36-35):(i*36)]
  fullKin$X2[(i*36-35):(i*36)] <- v2+(8*(i-1))
  v3 <- fullKin$X3[(i*36-35):(i*36)]
  fullKin$X3[(i*36-35):(i*36)] <- v3+(8*(i-1))
}

write.table(mkin,"monster/kin.txt",row.names=FALSE,quote=FALSE,sep="\t")

# now can call monster
# cd /projects/geneva/geneva_sata/caitlin/keats_x/monster/MONSTER
# ./MONSTER -c -p ../pheno.txt -g ../geno.txt -s ../SNP.txt -k ../kin.txt -r    

res <- read.table("monster/MONSTER/MONSTER.out",header=TRUE)
# ok, so have these pvalues to compare

# need to run KEATS and burden with just the auto kinship

keats.fam <- keats(MAF,geno,kAuto,scan1,estVar=vcAuto)

# doesn't matter what is on diags for keats.bt
keats.bt <- assocTestMixedModel(genoData,cholSigmaInv=cholSigAuto,outcome="y",test="Score")

source("famSKAT.R")
fullkins <- makekinship(ped[,1],ped[,2],ped[,3],ped[,4])
out <- famSKAT(phenotype=pheno$y, genotypes=t(geno), id=pheno$scanID,fullkins=fullkins)

source("keats_optimal.R")
keatso <- keats_optimal(MAF,geno,kAuto,scan1,estVar=vcAuto)

obj <- SKAT_Null_Model(y~1,out_type="C")
skato <- SKAT(t(geno),obj,r.corr=seq(from=0,to=1,by=0.1),method="liu")

skato
keats.fam[1:3]
keatso
keats.bt
out[1:3]
res

## first check:
# keatso rho=0 should be equal to keats.fam pvalue - YES
# keatso rho=0 should be near to out pvalue - YES
# keats.fam pvalue should be near to out pvalue - YES

# keatso rho=1 should be equal to keats.bt - YES

# in general, MONSTER is many orders of magnitude different from keats and famSKAT
# res pvalue should be near to keatso pvalue - NOT AT ALL
# res famSKAT pvalue should be near to out pvalue - NOT AT ALL
# res famSKAT pvalue should be near to keats.fam pvalue - NOT AT ALL


# monster is called with kinship coef, 0 on diag
# famSKAT is called with kinship coef, 0.5 on diag
# KEATS is called with 2*kinship coef, 0.5 on diag


rm(list=ls())

#####


library(MASS)
library(GWASTools)
library(SNPRelate)
library(OLGAanalysis)
library(survey)
library(kinship)
library(pbivnorm)
library(SKAT)
setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
source("../mlm_x/allele_drop_functions.R")
source("hapsim_functions.R")
source("keats.R")
source("keats_optimal.R")
source("/projects/geneva/geneva_sata/caitlin/keats_x/estVarComp_combined_v4.R")
wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)


nvar <- 20
maf <- runif(nvar,min=0.005,max=0.05)

ld <- getobj("LDMat_ceu_chr9_100snps.RData")
ld <- ld[1:nvar,1:nvar]
ld[lower.tri(ld,diag=FALSE)] <- t(ld)[lower.tri(ld,diag=FALSE)]
ld <- makepd(ld)

#effVar <- c(0.5,0,0.5)
effVar <- c(0.1,0.1,0.8)
n <- 2000
nvar <- as.integer(nvar)

SEX <- c("F","M","F","F","M","M","M","F")
sex <- rep(SEX,(n/8))

# simulate genotypes now
genoMat <- matrix(NA,nrow=nvar,ncol=n) # snp x sample; should be sample x snp
colsToGeno <- split(1:n,rep(1:n,each=8,length=n))
geno <- do.call(cbind,lapply(colsToGeno,function(x){genoMat[,x]=Family_alleles_Nmarker_unrel(maf,nvar,sex,ld)}))
MAF<-rowMeans(geno, na.rm = TRUE)/2 # this needs to change for X chr SNPs

noise <- rnorm(n,mean=0,sd=1)

vars <- 1:nvar
notCausal <- sample(vars,nvar*as.numeric(effVar[3]))
negVar <- sample(vars[-notCausal],nvar*as.numeric(effVar[1]))
posVar <- sample(vars[-c(notCausal,negVar)],nvar*as.numeric(effVar[2]))

c <- as.numeric(0.2)
beta <- rep(0,nvar)
beta[posVar] <- c*abs(log10(MAF[posVar]))
beta[negVar] <- -c*abs(log10(MAF[negVar]))
beta[MAF==0] <- 0
bgeno <- beta*geno
y <- colSums(bgeno)#+noise
pheno <- data.frame(scanID=1:n,y=y,sex=sex)

scan1 <- ScanAnnotationDataFrame(pheno)

kinAuto <- getobj("../mlm_x/1000Peds_8000unrel_autoKinship.RData")
sex <- rep(SEX,(n/8))
sex.int <- sex
sex.int[sex.int=="M"] <- 1
sex.int[sex.int=="F"] <- 2
kinAuto <- matrix(kinAuto,nrow=8000,ncol=8000)
kinAuto <- kinAuto[1:n,1:n]

colnames(kinAuto) <- 1:n
kinAuto <- 2*kinAuto
diag(kinAuto) <- 0.5 # needs to be 1+inbd=1 or inbd=0; MONSTER uses inbd
kAuto <- list(kinshipAutos=kinAuto)

# make up my own estVar results, only want to use ident as sigma matrix
rm(maf)
maf <- MAF
rm(MAF)

#####
## ok, run some skato functions manually and compare w mine
source("SKAT/R/SKAT_Optimal.R")
obj.res <- SKAT_Null_Model(y~1,out_type="C")
Z <- t(geno)
kernel="linear.weighted"
method="liu"
weights=wuweights(maf)
r.corr=seq(from=0,to=1,by=0.1)
is_check_genotype=TRUE; is_dosage = FALSE; missing_cutoff=0.15; estimate_MAF=1
SetID = NULL; out.z=NULL

res=obj.res$res
X1 <- obj.res$X1
s2=obj.res$s2
res.out=obj.res$res.out
n.Resampling=obj.res$n.Resampling

# now run SKAT_Optimal line by line
n<-dim(Z)[1]
p.m<-dim(Z)[2]	
r.all=r.corr
n.r<-length(r.all)
Z = t(t(Z) * (weights))
## Z is t(geno)*weights

out.Q<-SKAT_Optimal_Get_Q(Z, res, r.all, n.Resampling, res.out)
temp<-t(res) %*% Z
r.corr <- r.all[1]
Q1<-(1-r.corr) * rowSums(temp^2)
Q2<-r.corr * p.m^2 * rowMeans(temp)^2
Q1+Q2

Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) / s2

### compare w my q stats
weights<-wuweights(maf)
G <- t(geno*weights) # is sample x snp now, after t()
n <- nrow(scan1)
nvar <- nrow(geno)
rhoVals <- seq(from=0,to=1,by=0.1)

mod <- lm(y~1,data=pheno)
resids <- mod$residuals

ptA <- G
ptB <- t(G)
H <- X1 %*% solve(t(X1) %*%  X1) %*% t(X1) 

qvals <- data.frame(matrix(NA,nrow=length(rhoVals),ncol=4))
colnames(qvals)<- c("stat","pvalue","rho","qmin")
qvals$rho <- rhoVals
lambdas <- vector("list",length(rhoVals))
for(i in seq_along(rhoVals)){
  rho <- rhoVals[i]
  
  R <- diag(x=1,nrow=length(maf),ncol=length(maf))
  R <- (1-rho)*R + rho*(t(t(rep(1,length(maf))))%*%t(rep(1,length(maf))))
  QMine <- t(resids) %*% ptA %*% R %*% ptB %*% resids
  
  tmpEs <- eigen(ptB %*% (diag(n)-H) %*% G %*% R, symmetric=TRUE, only.values=TRUE)
  ourEs <-   tmpEs$values[tmpEs$values>1e-6*tmpEs$values[1]]
  
  lambdas[[i]] <- ourEs
  pfamskat<-pchisqsum(QMine,rep(1,length(ourEs)),ourEs,lower=FALSE,method="sad")
  pfamskat<-min(pfamskat,1)
  pfamskat<-max(pfamskat,0)
  
  qvals$pvalue[i] <- pfamskat
  qvals$stat[i] <- QMine
} # now we have stat and pvalue for each rho value

out.Q; qvals
# GREAT! they're the same
# skato randomly divides by 2. should i do that?

# use my qvals to move forward
Q.all <- matrix(qvals$stat,nrow=1)

# check evals
n.q<-dim(Q.all)[1]
p.m <- dim(Z1)[2]
Z1 <- Z - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z)
# the skato function calls with Z1/sqrt(2) -- is this normalized??
source("SKAT/R/Function.R")
lambda.all<-list()
for(i in 1:n.r){
  r.corr<-r.all[i]
  R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
  L<-chol(R.M,pivot=TRUE)
  Z2<- Z1 %*% t(L)
  K1<-t(Z2) %*% Z2
  
  lambda.all[[i]]<-Get_Lambda(K1)
}

###
# check with my lambda calculation
cbind(lambdas[[1]],lambda.all[[1]])
# EXACT! YESSS

# use our lambdas for forging ahead
lambda.all <- lambdas
Each_Info<-SKAT_Optimal_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
Each_Info; qvals # nice! so pvalues are pretty close too

# need to calculate qmin values now; see how they compare
T <- min(qvals$pvalue)

for(i in seq_along(rhoVals)){
  rho <- rhoVals[i]
  evals <- lambdas[[i]]
  qvals$qmin[i] <- our_qchisqsum(T,evals)
}

cbind(qvals$qmin,as.vector(Each_Info[[3]])) # also pretty close!

## now just need to find the overall pvalue for the test stat T
pmin.q<-Each_Info$pmin.q
pmin<-Each_Info$pmin
param.m <- SKAT_Optimal_Param(Z1,r.all)
SKAT_Optimal_PValue_Liu(pmin.q[1,],param.m,r.all, pmin) # 0.004550264

# get mine
ourZ <- t(diag(n)-H)  %*% G 
zbar <- rowMeans(ourZ)
M <- zbar %*% solve(crossprod(zbar)) %*% t(zbar)
ident <- diag(n)
eigenMat <- crossprod(ourZ, (ident-M)) %*% ourZ
tmplambdas <- eigen(eigenMat, symmetric=TRUE, only.values=TRUE)
lambdas <- tmplambdas$values[tmplambdas$values > 1e-6*tmplambdas$values[1]]

mua <- sum(lambdas)
sig2a <- 2*sum(lambdas^2)

trMatrix <- crossprod(ourZ,M) %*% ourZ %*% eigenMat
sig2xi <- 4*sum(diag(trMatrix))

kera <- sum(lambdas^4)/(sum(lambdas^2)^2) * 12
ldf <- 12/kera

# now need to find the tau(rho) values
tau <- nvar^2 * rhoVals * crossprod(zbar) + ((1-rhoVals)/crossprod(zbar))*sum(apply(ourZ,2,function(x){(t(zbar)%*%t(t(x)))^2}))

# now need to find the min{(qmin(rho)-rho*chisq_1)/(1-rho)}
# have everything to do the integrate function now
otherParams <- c(mu=mua, degf=ldf, varia=sig2a+sig2xi) 

re <- integrate(integrateFxn, lower=0, upper=40, subdivisions=2000,
                qmin=qvals$qmin,otherParams,tau=tau,rhoVals=rhoVals,abs.tol=10^-25)
pvalue <- 1-re[[1]]

if(T * length(rhoVals) < pvalue){
  pvalue <- T*length(rhoVals)
}
pvalue # 0.004858426
# THE SAME!!!! YESSSSSSS

rm(list=ls())


#####

