##################################################
# SKAT for Quantitative Traits in Family Samples #
# Han Chen (hanchen@bu.edu)                      #
# Version 1.8 Released 04-05-2013                #
##################################################

# NOTES:
# This function (famSKAT) is used to perform SKAT analysis (Wu et al., 2011) for quantitative traits in family data
# It requires the OLD "kinship" R package in the archive:
# http://cran.r-project.org/src/contrib/Archive/kinship/kinship_1.1.0-23.tar.gz
# After downloading the package, you can install it using
# install.packages(pkgs="kinship_1.1.0-23.tar.gz",lib="your_directory/")
# Then you can use it by specifying
# library(kinship,lib.loc="your_directory/")

# This function (famSKAT) requires the R package "survey" to compute the p-value using Kuonen's saddlepoint method (Kuonen, 1999)
# If you would like to compute the p-value using Davies' method (Davies, 1980), you would need the R package "CompQuadForm" as well

library(kinship) # Required. If you also have newer versions of kinship installed, make sure you are using the correct package by specifying lib.loc
library(survey) # Required if you use Kuonen's method to compute the p-value (recommended)
library(CompQuadForm) # This is optional. You can turn it off if you do not use Davies' method to compute the p-value

wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
## If you want to use wuweights, please make sure that your major allele is coded as 0

##### ATTENTION:
##### This program will not run without a full pedigree dataset "fullped", which contains all individuals in the pedigree
##### Even if you only include a subset of your sample in the analysis, you should still use the full pedigree to compute the kinship matrix
##### Otherwise higher degrees of relatedness (for example cousins) will be incorrectly ignored
##### The kinship matrix fullkins should include ALL individuals in your full pedigree
# fullped<-"YOUR FULL PEDIGREE dataset with famid, id, fa, mo"
# fullkins<-makekinship(fullped$famid, fullped$id, fullped$fa, fullped$mo)

# FUNCTION ARGUMENTS:
# phenotype: A vector of quantitative trait in the analysis. The order should match the vector id.
# genotypes: A matrix of genotypes. The order of rows should match the vector id.
# id: A vector of id. Please make sure that all these id's are included in your full pedigree dataset.
# fullkins: Kinship matrix calculated from the full pedigree. Described above.
# covariates: A matrix of covariates. The order of rows should match the vector id. Default NULL.
# h2: A positive number between 0 and 1. A priori knowledge of heritability. If missing, two variance component parameters (random effects of relatedness, and the residual variance) will be estimated from the linear mixed effects model. Default NULL.
# sqrtweights: A vector with the length equal to the number of variants in the test, or a function of the MAF. Default wuweights function as described in Wu's SKAT paper
# binomialimpute: Logical. If FALSE (default) all missing genotypes will be replaced by 0. If TRUE, missing genotypes will be imputed randomly as a binomial distribution with sample MAF, as the original SKAT function does
# method: "Kuonen" (default) or "Davies". Method of computing the p-value.
# acc: Accuracy of numerical integration. Only used in Davies' method. Default NULL.

famSKAT<-function(phenotype, genotypes, id, fullkins, covariates=NULL, h2=NULL, sqrtweights=wuweights, binomialimpute=FALSE, method="Kuonen", acc=NULL) {
	# Regular checks
	if(class(phenotype)!="numeric" && class(phenotype)!="integer") stop("phenotype should be a numeric vector!")
	n<-length(phenotype)
	if(is.data.frame(genotypes)) genotypes<-as.matrix(genotypes)
	if(!is.matrix(genotypes)) stop("genotypes should be a matrix!")
	if(nrow(genotypes)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
	if(!is.vector(id)) stop("id should be a vector!")
	if(any(duplicated(id))) stop("Duplicated id exists. Check your data...")
	if(length(id)!=n) stop("Number of individuals inconsistent between phenotype and id. Check your data...")
	if(!is.bdsmatrix(fullkins)) stop("fullkins should be a kinship matrix from makekinship function!")
	if(!is.null(covariates)) {
		covariates<-as.matrix(covariates)
		if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")
	}
	if(!is.null(h2)) {
		if(h2<0 || h2>1) stop("Heritability should be between 0 and 1")
	}
	if(method=="Davies") {
		if(is.null(acc)) acc<-1e-6
	} else if(method!="Kuonen") {
		stop("Method should be either Kuonen or Davies!")
	}

	# Remove missing phenotype and covariates
	missidx<-is.na(phenotype)
	if(!is.null(covariates)) {
		missidx<-missidx | apply(is.na(covariates), 1, any)
	}
	phenotype<-phenotype[!missidx]
	genotypes<-as.matrix(genotypes[!missidx,])
	id<-id[!missidx]
	if(!is.null(covariates)) covariates<-covariates[!missidx,]
	n<-length(phenotype)

	# Check weights and fullkins
	Z<-genotypes
	nZ<-ncol(Z)
	MAF<-colMeans(Z, na.rm = TRUE)/2
	if(is.function(sqrtweights)) {
		weights<-sqrtweights(MAF)
	} else {
		if(!is.vector(sqrtweights)) stop("Check the class of your variable sqrtweights: should be function or vector!")
		if(length(sqrtweights)!=nZ) stop("Number of variants inconsistent between genotypes and sqrtweights. Check your data...")
		weights<-sqrtweights
	}
	tmpidx<-!is.na(match(dimnames(fullkins)[[1]], id))
	tmpkins<-fullkins[tmpidx, tmpidx]
	tmpidx<-match(id, dimnames(tmpkins)[[1]])
	if(any(is.na(tmpidx))) stop("Some id not exist in fullkins. Check your full pedigree data...")
	kins<-as.matrix(tmpkins)[tmpidx, tmpidx]

	# impute missing genotypes
	for(pp in 1:nZ) {
		IDX<-which(is.na(Z[,pp]))
		if(length(IDX)>0) {
			if(binomialimpute) { # random imputation based on binomial distribution
				Z[IDX,pp]<-rbinom(length(IDX), 2, MAF[pp])
			} else { # default: set missing to 0
				Z[IDX,pp]<-0
			}
		}
	}
	# now Z should not have any missing values
	G<-t(t(Z)*weights)
	if(is.null(h2)) { # default: let lmekin estimate both variance components for the null model
		# null model
		if(is.null(covariates)) {
			X<-matrix(rep(1,n),n,1)
			tmpdata<-data.frame(phenotype,id)
			nullmod<-lmekin(phenotype~1,tmpdata,random=~1|id,varlist=list(tmpkins))
		} else {
			X<-cbind(rep(1,n),as.matrix(covariates))
			tmpdata<-data.frame(phenotype,id,covariates)
			nc<-ncol(covariates)
			if(nc==1) {
				exprs<-"phenotype ~ covariates"
			} else {
				exprs<-paste("phenotype ~", paste(names(tmpdata)[-c(1,2)],collapse=" + "))
			}
			nullmod<-lmekin(as.formula(exprs),tmpdata,random=~1|id,varlist=list(tmpkins))
		}
		res<-nullmod$res
		SIGMA<-nullmod$theta[1]*2*kins+nullmod$theta[2]*diag(n)
	} else { # provide an estimate of heritability externally, only estimate the total variance
		Om<-h2*2*kins+(1-h2)*diag(n)
		Om_i<-solve(Om)
		if(is.null(covariates)) {
			X<-matrix(rep(1,n),n,1)
		} else {
			X<-cbind(rep(1,n),as.matrix(covariates))
		}
		res<-phenotype-as.numeric(X%*%solve(t(X)%*%Om_i%*%X)%*%t(X)%*%Om_i%*%phenotype)
		s2<-as.numeric(t(res)%*%Om_i%*%res)/(n-ncol(X))
		SIGMA<-Om*s2
	}
	SIGMAi<-solve(SIGMA)
	Q<-t(res) %*% SIGMAi %*% G %*% t(G) %*% SIGMAi %*% res
	P<-SIGMA - X %*% solve(t(X) %*% SIGMAi %*% X) %*% t(X)
	eig<-eigen(t(G) %*% SIGMAi %*% P %*% SIGMAi %*% G, symmetric=T, only.values=T)
	evals<-eig$values[eig$values>1e-6*eig$values[1]]
	if(method=="Davies") {
		tmpout<-davies(Q, evals, acc=acc)
		pfamskat<-tmpout$Qq
		message<-tmpout$ifault
	} else {
		pfamskat<-pchisqsum(Q,rep(1,length(evals)),evals,lower=F,method="sad")
		message<-0
	}
	pfamskat<-min(pfamskat,1)
	pfamskat<-max(pfamskat,0)
	return(list(pvalue=pfamskat, N=n, Nmarkers=nZ, Nvariants=sum(MAF>0 & MAF<1), message=message))
}
### OUTPUT:
### pvalue: famSKAT p-value
### N: sample size in the analysis (missing phenotype and covariates removed)
### Nmarkers: number of genetic markers in the test
### Nvariants: number of genetic variants in the test (monomorphic markers removed, should be smaller than or equal to Nmarkers)
### message: error message
###          0: no error
###          1: requested accuracy could not be obtained
###          2: round-off error possibly significant
###          3: invalid parameters
###          4: unable to locate integration parameters

# References:
# Chen H, Meigs JB, Dupuis J. 2013. Sequence kernel association test for quantitative traits in family samples. Genet Epidemiol 37: 196-204.
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Kuonen D. 1999. Saddlepoint approximations for distributions of quadratic forms in normal variables. Biometrika 86:929-935.
# Schifano ED, Epstein MP, Bielak LF, Jhun MA, Kardia SL, Peyser PA, Lin X. 2012. SNP set association analysis for familial data. Genet Epidemiol 36: 797-810.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

### Example
# set.seed(12345)
# fullped<-data.frame(famid=rep(1:250,each=4),id=10001:11000,fa=rep(0,1000),mo=rep(0,1000))
# fullped$fa[(1:250)*4-1]<-fullped$fa[(1:250)*4]<-(1:250)*4+9997
# fullped$mo[(1:250)*4-1]<-fullped$mo[(1:250)*4]<-(1:250)*4+9998
# pheno<-data.frame(id=10001:11000,y=rnorm(1000),sex=rep(1:2,500),bmi=rnorm(1000,25,2))
# geno<-data.frame(id=10001:11000,snp1=sample(0:2,1000,replace=T,prob=c(0.81,0.18,0.01)),snp2=sample(0:2,1000,replace=T,prob=c(0.81,0.18,0.01)),snp3=sample(0:2,1000,replace=T,prob=c(0.81,0.18,0.01)))
# out<-famSKAT(phenotype=pheno$y, genotypes=as.matrix(geno[,-1]), id=pheno$id, fullkins=makekinship(fullped$famid, fullped$id, fullped$fa, fullped$mo), covariates=as.matrix(pheno[,c("sex","bmi")]))
# print(out)





