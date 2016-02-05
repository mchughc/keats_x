
# R script to run KEATS-X on simulated datasets

## Contents:
# 1. Calculate chr 9 LD structure for haps in MXL, CEU, YRI, ASW
# 2. Calculate X chr LD structure for haps in MXL, CEU, YRI, ASW
# 3. Test code to simulate haps with CEU LD structure
# 4. Run nullSims
# 5. Make AR-1 LD matrix; run null sims
# 6. Run powerSims
# 7. More sims
# 8. Histogram of pvalues from 8ped null sims
# 9. Visualize LD structure
# 10. More tyIerr sims
# 11. Combine all null simulations for tyI error table
# 12. Process power results 25/25/50
# 13. Process power results 40/40/20, 50/0/50
# 14. Process power results 40/40/20; diff sigmaA, sigmaX values
# 15. Process power results, keatso x chr, 8ped
# 16. Power results, keatso, chrX for unrelated samples
# 17. Power results, keatso, chrX and auto for different mixes of +/-
# 18. Power results, keatso, chrX and auto for different mixes of +/-
# 19. Compare keatso to keats.bt and keatso rho=0 (keats.skat)
# 20. Power results, keatso, chrX and auto for different mixes of +/- for 8pedFem
# 21. Power results, keatso, chrX and auto for different mixes of +/- for unrelated pedigree
# 22. Make type I error tables for all methods considered
# 23. Make type I error tables for various sig2x, sig2a values
# 24. Get SNPs in known RBC hits from HCHS/SOL data
# 25. Parse type I error results, 50 variants
# 26. Run KEATSO on each RBC gene region in HCHS/SOL
# 27. Power results, keatso, chrX and auto for different mixes of +/- for 8pedFem
# 28. Hist of p-values corresponding to #23.
# 29. Hist of p-values corresponding to #25.
# 30. Type I error results, 8pedFem, auto & x, comparing keatso to keatsbt and keats
# 31. Group SNPs into genes, genome-wide
# 32. Run KEATSO on HCHS/SOL genome-wide
# 33. Power results, keatso, chrX and auto for different mixes of +/- for 8ped
# 34. Hist of optimal rho values, null sims 8pedFem
# 35. KEATSO on HCHS/SOL genome-wide QQ plots
# 36. KEATSO on HCHS/SOL genome-wide single test p-value comparisons
# 37. Make a new HCHS/SOL gene list, NOT excluding SNPs with >5% MAF
# 38. HCHS/SOL genome-wide using MONSTER, and SKAT on unrelateds
# 39. KEATSO on HCHS/SOL genome-wide, NOT filtering out SNPs with >5% MAF
# 40. KEATSO pvalues for genome-wide sig hits, compared to KEATS pvalues
# 41. KEATSO on rare SNPs, platelet HCHS/SOL 
# 42. KEATSO on rare + common SNPs, platelet HCHS/SOL 
# 43. Look at how the variants are mapped to genes
# 44. Follow up on hits from KEATSO rare analysis
# 45. Make table of HBB variants to send to RBC email list
# 46. More power sims to add to dissertation
# 47. Make heatmaps of LD between variants in chr 16, x chr RBC gene hits
# 48. Find GATS and AHI1 rare gene analysis p-values



#####
# 1. Calculate chr 9 LD structure for haps in MXL, CEU, YRI, ASW

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")

# sample annot for impute2 phased files, 1000G phase 1 v3
samp <- read.table("/projects/geneva/gcc-fs2/1000Genomes/20130502/IMPUTE2/1000GP_Phase3/1000GP_Phase3.sample",
                   header=TRUE,as.is=TRUE,sep=" ")
dim(samp); head(samp) # 2504 4
table(samp$POP)

mxl <- samp[samp$POP=="MXL",]
dim(mxl) # 64 4
mxl$sex[mxl$SEX=="male"] <- "M"
mxl$sex[mxl$SEX=="female"] <- "F"

# get phased info for chr 9 for these samples
dat <- read.table(gzfile("/projects/geneva/gcc-fs2/1000Genomes/20130502/IMPUTE2/1000GP_Phase3/1000GP_Phase3_chr19.hap.gz"),
                  header=F,as.is=T,sep=" ",nrow=1000)
dim(dat) # 1000 5008

# find the cols of dat that correspond to the 64 MXL samples
datPop <- rep(samp$POP,each=2)
haps9 <- dat[,datPop=="MXL"]
dim(haps9) # 1000 128

afreq <- apply(haps9,1,function(x){sum(x)/128})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 188 128, so snp x sample (x2)
afreq <- afreq[afreq!=0]

haps9 <- haps9[-c(3,19),]
afreq <- apply(haps9,1,function(x){sum(x)/128})

# want to calculate pairwise LD for these SNPs
library(pbivnorm); library(mvtnorm)
source("hapsim_functions.R")

C <- cor(t(haps9[1:100,]))
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_mxl_chr9_100snps.RData")

# now get for ceu samples
ceu <- samp[samp$POP=="CEU",]
dim(ceu) # 99 4
ceu$sex[ceu$SEX=="male"] <- "M"
ceu$sex[ceu$SEX=="female"] <- "F"

datPop <- rep(samp$POP,each=2)
haps9 <- dat[,datPop=="CEU"]
dim(haps9) # 1000 198

afreq <- apply(haps9,1,function(x){sum(x)/198})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 186 198, so snp x sample (x2)
afreq <- afreq[afreq!=0]

haps9 <- haps9[-c(27),]
afreq <- apply(haps9,1,function(x){sum(x)/198})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 185 198, so snp x sample (x2)
afreq <- afreq[afreq!=0]

C <- cor(t(haps9[1:100,]))
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_ceu_chr9_100snps.RData")

# now for yri samples
yri <- samp[samp$POP=="YRI",]
dim(yri) # 108 4
yri$sex[yri$SEX=="male"] <- "M"
yri$sex[yri$SEX=="female"] <- "F"

datPop <- rep(samp$POP,each=2)
haps9 <- dat[,datPop=="YRI"]
dim(haps9) # 1000 216

afreq <- apply(haps9,1,function(x){sum(x)/216})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 241 216, so snp x sample (x2)
afreq <- afreq[afreq!=0]

haps9 <- haps9[-c(10),]
afreq <- apply(haps9,1,function(x){sum(x)/216})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 240 216, so snp x sample (x2)
afreq <- afreq[afreq!=0]

C <- cor(t(haps9[1:100,]))
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_yri_chr9_100snps.RData")

# finally for asw samples
asw <- samp[samp$POP=="ASW",]
dim(asw) # 61 4
asw$sex[asw$SEX=="male"] <- "M"
asw$sex[asw$SEX=="female"] <- "F"

datPop <- rep(samp$POP,each=2)
haps9 <- dat[,datPop=="ASW"]
dim(haps9) # 1000 122

afreq <- apply(haps9,1,function(x){sum(x)/122})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 220 122, so snp x sample (x2)
afreq <- afreq[afreq!=0]

haps9 <- haps9[-c(7),]
afreq <- apply(haps9,1,function(x){sum(x)/122})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 219 122, so snp x sample (x2)
afreq <- afreq[afreq!=0]

C <- cor(t(haps9[1:100,]))
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_asw_chr9_100snps.RData")

rm(list=ls())


#####
# 2. Calculate X chr LD structure for haps in MXL, CEU, YRI, ASW

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")

# sample annot for impute2 phased files, 1000G phase 1 v3
samp <- read.table("/projects/geneva/gcc-fs2/1000Genomes/20101123/IMPUTE/ALL_1000G_phase1integrated_v3_impute_mar2012/ALL_1000G_phase1integrated_v3.sample",
                   header=TRUE,as.is=TRUE,sep=" ")
dim(samp); head(samp) # 2504 4
table(samp$population)

mxl <- samp[samp$population=="MXL",]
dim(mxl) # 66 4

# get phased info for X chr for these samples
dat <- read.table(gzfile("/projects/geneva/gcc-fs2/1000Genomes/20101123/IMPUTE/ALL_1000G_phase1integrated_v3_impute_mar2012/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz"),
                  header=F,as.is=T,sep=" ",nrow=1000)
dim(dat) # 1000 2184
dat[dat=="-"] <- NA
dat <- as.matrix(dat)
dat <- matrix(as.integer(dat),nrow=1000,ncol=2184)

# find the cols of dat that correspond to the 64 MXL samples
datPop <- rep(samp$population,each=2)
haps9 <- dat[,datPop=="MXL"]
dim(haps9) # 1000 132

afreq <- apply(haps9,1,function(x){sum(x,na.rm=T)/101})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 456 132, so snp x sample (x2)
afreq <- afreq[afreq!=0]

# want to calculate pairwise LD for these SNPs
library(pbivnorm)
source("hapsim_functions.R")

C <- cor(t(haps9[1:100,]),use="complete.obs")
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_mxl_chrX_100snps.RData")

##
# now get for ceu samples
ceu <- samp[samp$population=="CEU",]
dim(ceu) # 85 4

datPop <- rep(samp$population,each=2)
haps9 <- dat[,datPop=="CEU"]
dim(haps9) # 1000 170

afreq <- apply(haps9,1,function(x){sum(x,na.rm=T)/125}) # denom here is # haplotypes -- diff than 170*2 since males only have one
haps9 <- haps9[afreq!=0,]
dim(haps9) # 442 170, so snp x sample (x2)
afreq <- afreq[afreq!=0]

C <- cor(t(haps9[1:100,]),use="complete.obs")
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

afreqS <- afreq[1:100]
save(afreqS,file="hapFreq_ceu_100snps.RData")
save(ld9,file="LDMat_ceu_chrX_100snps.RData")

##
# now for yri samples
yri <- samp[samp$population=="YRI",]
dim(yri) # 88 4

datPop <- rep(samp$population,each=2)
haps9 <- dat[,datPop=="YRI"]
dim(haps9) # 1000 176

afreq <- apply(haps9,1,function(x){sum(x,na.rm=T)/133})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 647 176, so snp x sample (x2)
afreq <- afreq[afreq!=0]

C <- cor(t(haps9[1:100,]),use="complete.obs")
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_yri_chrX_100snps.RData")

##
# finally for asw samples
asw <- samp[samp$population=="ASW",]
dim(asw) # 61 4

datPop <- rep(samp$population,each=2)
haps9 <- dat[,datPop=="ASW"]
dim(haps9) # 1000 122

afreq <- apply(haps9,1,function(x){sum(x,na.rm=T)/98})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 632 122, so snp x sample (x2)
afreq <- afreq[afreq!=0]

C <- cor(t(haps9[1:100,]),use="complete.obs")
nloci=100
P <- afreq[1:100]
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

rownames(ld9) <- rownames(afreq)[1:100]
colnames(ld9) <- rownames(afreq)[1:100]

save(ld9,file="LDMat_asw_chrX_100snps.RData")

rm(list=ls())


#####
# 3. Test code to simulate haps with CEU LD structure

library(MASS)
library(GWASTools)
library(SNPRelate)
library(pbivnorm)
setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
source("../mlm_x/allele_drop_functions.R")
source("hapsim_functions.R")

n <- 8000

SEX <- c("F","M","M","M","F","M","M","M")
sex <- rep(SEX,(n/8))

kinX <- getobj("../mlm_x/1000Peds_8ped_xKinship.RData")
kinAuto <- getobj("../mlm_x/1000Peds_8ped_autoKinship.RData")

kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
kinX <- matrix(kinX,nrow=n,ncol=n)

ldX <- getobj("LDMat_ceu_chrX_100snps.RData")
ldX[lower.tri(ldX,diag=FALSE)] <- t(ldX)[lower.tri(ldX,diag=FALSE)]

ldAuto <- getobj("LDMat_ceu_chr9_100snps.RData")
ldAuto[lower.tri(ldAuto,diag=FALSE)] <- t(ldAuto)[lower.tri(ldAuto,diag=FALSE)]

ldX <- makepd(ldX)
ldAuto <- makepd(ldAuto)

afreq <- getobj("hapFreq_ceu_100snps.RData")

genoMat <- founder_alleles_Nmarker_haplotypes(10000,afreq,100,ldX)
dim(genoMat) # 100 snps x 10000 samples

afreqObs <- apply(genoMat,1,function(x){sum(x,na.rm=T)/(10000)})

C <- cor(t(genoMat))
nloci=100
P <- afreqObs
Q <- qnorm(P)
null.mat <- matrix(0, nrow=nloci,ncol=nloci)
vmat <- covmat(as.integer(nloci),C,P,Q)
ld9 <- matrix(vmat, nrow=nloci, ncol=nloci)

cbind(afreq,afreqObs)
ld9[1:10,1:10]
ldX[1:10,1:10]
# YES! they look great!

rm(list=ls())


#####
# 4. Run nullSims

# submitted:
# cd /projects/geneva/geneva_sata/caitlin/keats_x/
# qsub -q olga.q -p -50 -t 1-2500 -tc 100 -N nullX batch_nullSim_8ped.sh
# qsub -q thornton.q -N null9 -t 1-2500 -tc 50 batch_nullSim_8ped.sh

#cat batch_nullSim_8ped.sh
# #!/bin/bash
# R --vanilla "--args nullSims_8ped_chrXceu/sim $SGE_TASK_ID ped8 FALSE LDMat_ceu_chrX_100snps.RData" </projects/geneva/geneva_sata/caitlin/keats_x/nullSim.R >/projects/geneva/geneva_sata/caitlin/keats_x/nullSim.Rout

# which calls nullSims.R
# can specify 8ped or 8pedFem, also if auto or x chr SNPs, and which LD matrix to simulate from

# stores results in nullSims_8ped_chr9ceu/
# and nullSims_8ped_chrXceu/

# cat *_res* >> allRes.txt
# cat *_maf* >> allMAF.txt

library(readr)
library(dplyr)
library(tidyr)
## need to do nullSims_8ped_fem next
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chrXceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method")
dat <- dat[!duplicated(dat$auto),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  gather(model,pval,-method) %>%
  group_by(model,method)
summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
#  model  method     tyI_err
#1  auto famSKAT 0.073170732
#2  auto   keats 0.076907691
#3     x famSKAT 0.004398241
#4     x   keats 0.004400440
#5  both famSKAT 0.011995202
#6  both   keats 0.013601360
# why is the auto model the only one calibrated??

summarize(tyI_group,tyI_err=sum(pval<0.005)/n())

summarize(tyI_group,n()) # keats 9999, famSKAT 2501

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chr9ceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method")
dat <- dat[!duplicated(dat$auto),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  gather(model,pval,-method) %>%
  group_by(model,method)
summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
#  model method    tyI_err
#1  auto  keats 0.05000500
#2     x  keats 0.00450045
#3  both  keats 0.01090109
# still super undercalibrated, so too conservative -- why??

summarize(tyI_group,tyI_err=sum(pval<0.005)/n())
#  model method    tyI_err
#1  auto  keats 0.00440044
#2     x  keats 0.00010001
#3  both  keats 0.00030003

summarize(tyI_group,n()) # 9999 iterations

rm(list=ls())


#####
# 5. Make AR-1 LD matrix; run null sims

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
ld <- diag(100)
ld <- 0.5^abs(row(ld)-col(ld))
save(ld,file="LDMat_05_ar1_100snps.RData")

rm(list=ls())

##
# qsub -q olga.q -N nullAR1 -t 1-10000 -p -50 batch_nullSim_8ped.sh
# which writes to  nullSims_8ped_chr9_ld05_ar1/sim
# 20 variants, 8ped, auto=TRUE, LDMat_05_ar1_100snps.RData

# qsub -q thornton.q -N nullARX -t 1-10000 -tc 50 batch_nullSim_8ped.sh
# which writes to nullSims_8ped_chrX_ld05_ar1/sim
# 20 variants, 8ped, auto=FALSE, LDMat_05_ar1_100snps.RData

# cat sim*_res.txt >> allRes.txt
##

library(readr)
library(dplyr)
library(tidyr)
## need to do nullSims_8ped_fem next
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method")
dat <- dat[!duplicated(dat$auto),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  filter(auto<=1) %>% # this filters out the Q stat rows
  gather(model,pval,-method) %>%
  group_by(model,method)

summarize(tyI_group,n()) # keats 9998

summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
#  model method    tyI_err
#1  auto  keats 0.04920984
#2     x  keats 0.07131426
#3  both  keats 0.05011002

summarize(tyI_group,tyI_err=sum(pval<0.005)/n())
#  model method    tyI_err
#1  auto  keats 0.00430086
#2     x  keats 0.00840168
#3  both  keats 0.00410082

rm(list=ls())


#####
# 6. Run powerSims

# cd /projects/geneva/geneva_sata/caitlin/keats_x
# qsub -q olga.q -p -50 -N pow9 -t 1-10000 batch_powerSim.sh
# which writes to powerSims_8ped_chr9_ld05_ar1/sim442
# 20 variants, 8ped, auto=TRUE, effVar=c(0.4,0.4,0.2), c=0.3

# qsub -q olga.q -p -50 -N powX -t 1-10000 batch_powerSim.sh
# which writes to powerSims_8ped_chrX_ld05_ar1/sim442
# 20 variants, 8ped, auto=FALSE, effVar=c(0.4,0.4,0.2), c=0.3

# cat sim*_res.txt >> allRes442_c03.txt
##

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

## need to do nullSims_8ped_fem next
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/allRes442_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method","snp")
dat <- dat[!duplicated(dat$auto),]

pow <- dat %>% filter(!is.na(snp)) 
metric <- rep(c("Q","pval"),nrow(pow)/2)
pow <- mutate(pow,"metric"=metric)

pow <- pow %>%
  filter(metric=="pval") %>% # this filters out the Q stat rows
  gather(model,pval,-c(method,snp,metric)) %>%
  group_by(model)

summarize(pow,n()) # 9999 of each type; as expected

# merge in null iterations
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                  delim=" ", col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method")
dat <- dat[!duplicated(dat$auto),]
dat <- filter(dat,!is.na(method))

metric <- rep(c("Q","pval"),nrow(dat)/2)
dat <- mutate(dat,"metric"=metric)

nullSims <- dat %>%
  filter(metric=="pval") %>%
  gather(model,pval,-c(method,metric)) %>%
  group_by(model)
summarize(nullSims,n()) # 9998 of each, great!

nullSims <- mutate(nullSims,pval=as.numeric(pval))
pow <- mutate(pow,pval=as.numeric(pval))

alpha <- seq(from=1e-6,to=0.05,by=5e-06)
fp <- data.frame(matrix(NA,nrow=3,ncol=length(alpha)))
tp <- fp
colnames(fp) <- colnames(tp) <- 1:ncol(fp)

for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pval<alpha[i])/9998)
  fp[1,i] <- as.numeric(s[1,2])
  fp[2,i] <- as.numeric(s[2,2])
  fp[3,i] <- as.numeric(s[3,2])

  s <- summarize(pow,sum(pval<alpha[i])/9999)
  tp[1,i] <- as.numeric(s[1,2])
  tp[2,i] <- as.numeric(s[2,2])
  tp[3,i] <- as.numeric(s[3,2])
}

toPl <- tbl_df(fp)
toPl <- toPl %>%
  mutate(model=c("auto","x","both")) %>%
  gather(alpha,fp_rate,-model) 
  
toPl2 <- tbl_df(tp)
toPl2 <- toPl2 %>%
  mutate(model=c("auto","x","both")) %>%
  gather(alpha,tp_rate,-model) 

allPl <- tbl_df(cbind(toPl,toPl2))

pdf("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/power_442_c02_10Kiters.pdf")
ggplot(allPl,aes(x=fp_rate,y=tp_rate,color=model)) + geom_line(size=1.3) + theme_bw() +
  xlab("False Positive Rate") + ylab("True Positive Rate") 
dev.off()

rm(list=ls())

##

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

## need to do nullSims_8ped_fem next
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes442_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method","snp")
dat <- dat[!duplicated(dat$auto),]

tyI_group <- dat %>% filter(!is.na(snp)) 
metric <- rep(c("Q","pval"),nrow(tyI_group)/2)
tyI_group <- mutate(tyI_group,"metric"=metric)

tyI_group <- tyI_group %>%
  filter(metric=="pval") %>% # this filters out the Q stat rows
  gather(model,pval,-c(method,snp,metric)) %>%
  group_by(model)

summarize(tyI_group,n()) # 10000 of each type; as expected

# merge in null iterations
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ", col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method")
dat <- dat[!duplicated(dat$auto),]
dat <- filter(dat,!is.na(method))

metric <- rep(c("Q","pval"),nrow(dat)/2)
dat <- mutate(dat,"metric"=metric)

nullSims <- dat %>%
  filter(metric=="pval") %>%
  gather(model,pval,-c(method,metric)) %>%
  group_by(model)
summarize(nullSims,n()) # 9527 of each, great!

alpha <- seq(from=1e-05,to=0.5,by=0.0001)

fp <- data.frame(matrix(NA,nrow=3,ncol=1+length(alpha)))
tp <- fp
colnames(tp) <- colnames(fp) <- 1:ncol(tp)

nullSims <- mutate(nullSims,pval=as.numeric(pval))
tyI_group <- mutate(tyI_group,pval=as.numeric(pval))

for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pval<alpha[i])/9527)
  fp[1,i] <- as.numeric(s[1,2])
  fp[2,i] <- as.numeric(s[2,2])
  fp[3,i] <- as.numeric(s[3,2])
  
  s <- summarize(tyI_group,sum(pval<alpha[i])/10000)
  tp[1,i] <- as.numeric(s[1,2])
  tp[2,i] <- as.numeric(s[2,2])
  tp[3,i] <- as.numeric(s[3,2])
}

toPl <- tbl_df(fp)
toPl <- toPl %>%
  mutate(model=c("auto","x","both")) %>%
  gather(alpha,fp_rate,-model) 

toPl2 <- tbl_df(tp)
toPl2 <- toPl2 %>%
  mutate(model=c("auto","x","both")) %>%
  gather(alpha,tp_rate,-model) 

allPl <- tbl_df(cbind(toPl,toPl2))

pdf("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/power_442_c02_10Kiters.pdf")
ggplot(allPl,aes(x=fp_rate,y=tp_rate,color=model)) + geom_line(size=1.3) + theme_bw() +
  xlab("False Positive Rate") + ylab("True Positive Rate") 
dev.off()

# paste sim442_c021*_beta.txt >> allBeta_442_c021.txt
# paste sim442_c022*_beta.txt >> allBeta_442_c022.txt
# paste sim442_c023*_beta.txt >> allBeta_442_c023.txt
# paste sim442_c024*_beta.txt >> allBeta_442_c024.txt
# paste sim442_c025*_beta.txt >> allBeta_442_c025.txt
# paste sim442_c026*_beta.txt >> allBeta_442_c026.txt
# paste sim442_c027*_beta.txt >> allBeta_442_c027.txt
# paste sim442_c028*_beta.txt >> allBeta_442_c028.txt
# paste sim442_c029*_beta.txt >> allBeta_442_c029.txt
# 
# paste allBeta_442_c02*.txt >> betas_442_c02.txt

betas <- read.table("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/betas_442_c02.txt",
                    sep="\t",header=TRUE,as.is=TRUE)
spc <- apply(betas,1,function(x){as.numeric(gsub(".* ","",x))}) # yess
dim(spc) # 10000 20
# now want to take the rowMeans which will be the mean beta value per iteration
betaMn <- apply(spc,1,mean)

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes442_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("auto","x","both","method","snp")
dat <- dat[!duplicated(dat$auto),]

dat <- dat %>% filter(!is.na(snp))
dat <- mutate(dat,metric = rep(c("Q","pval"),nrow(dat)/2))
dat <- dat %>% filter(metric=="pval")
dat <- dat %>% mutate(beta=betaMn)

pow <- dat %>%
  gather(model,pval,-c(method,snp,metric,beta)) 

pdf("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/powerBeta_442_c02_10Kiters.pdf")
ggplot(pow,aes(x=beta,y=-log10(as.numeric(pval)),color=model)) + geom_point() + theme_bw() +
  xlab("Beta") + ylab("Pvalue") 
dev.off()

rm(list=ls())


#####
# 7. More sims

##
# cd /projects/geneva/geneva_sata/caitlin/keats_x/
# qsub -q thornton.q -N nullAR1 -t 1-10000 -tc 50 batch_nullSim_8ped.sh
# which calls nullSim_v2.R & writes to nullSims_8pedfem_chr9_ld05_ar1/sim
# 20 variants, pedFem8, auto=TRUE, LDMat_05_ar1_100snps.RData

# qsub -q thornton.q -N nullX -t 1-10000 -tc 50 batch_nullSim_8ped.sh
# which calls nullSim.R & writes to nullSims_8pedfem_chrX_ld05_ar1/sim
# 20 variants, pedFem8, auto=FALSE, LDMat_05_ar1_100snps.RData

# cat sim*_res.txt >> allRes.txt
##


# cd /projects/geneva/geneva_sata/caitlin/keats_x
# qsub -q olga.q -p -50 -N pow9 -t 1-10000 batch_powerSim.sh
# which calls powerSim.R & writes to powerSims_8ped_chr9_ld05_ar1/sim442_c04
# 20 variants, 8ped, auto=TRUE, effVar=c(0.4,0.4,0.2), c=0.4

# qsub -q olga.q -p -50 -N powX -t 1-10000 batch_powerSim.sh
# which writes to powerSims_8ped_chrX_ld05_ar1/sim442
# 20 variants, 8ped, auto=FALSE, effVar=c(0.4,0.4,0.2), c=0.4

# cat sim*_res.txt >> allRes442_c04.txt
##

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
dat <- read_delim("nullSims_8pedfem_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","method")
dat <- dat[!duplicated(dat$famSKAT),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  filter(famSKAT<=1) %>% # this filters out the Q stat rows
  gather(model,pval,-method) %>%
  mutate(pval=as.numeric(pval)) %>%
  group_by(model,method)

summarize(tyI_group,n()) # keats 9999

summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
#  model method    tyI_err
#1  auto  keats 0.05010501
#2     x  keats 0.06440644
#3  both  keats 0.05030503

summarize(tyI_group,tyI_err=sum(pval<0.005)/n())
#  model method    tyI_err
#1  auto  keats 0.00420042
#2     x  keats 0.00600060
#3  both  keats 0.00410041

chr9 <- tyI_group

pdf("nullSims_8pedfem_chr9_ld05_ar1/pvals_hist.pdf",width=12)
ggplot(tyI_group,aes(x=pval)) + geom_histogram(binwidth=1/10) + facet_wrap(~model)
dev.off()

##
dat <- read_delim("nullSims_8pedfem_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","method")
dat <- dat[!duplicated(dat$famSKAT),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  filter(famSKAT<=1) %>% # this filters out the Q stat rows
  gather(model,pval,-method) %>%
  mutate(pval=as.numeric(pval)) %>%
  group_by(model,method)

summarize(tyI_group,n()) # keats 9999

summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
#  model method    tyI_err
#1  auto  keats 0.06780678
#2     x  keats 0.04860486
#3  both  keats 0.04760476

summarize(tyI_group,tyI_err=sum(pval<0.005)/n())
#  model method    tyI_err
#1  auto  keats 0.00860086
#2     x  keats 0.00520052
#3  both  keats 0.00540054

pdf("nullSims_8pedfem_chrX_ld05_ar1/pvals_hist.pdf",width=12)
ggplot(tyI_group,aes(x=pval)) + geom_histogram(binwidth=1/10) + facet_wrap(~model)
dev.off()

chrX <- tyI_group

### 
# make pvalue histograms of just famSKAT & keats, and combine xchr & chr 9 results into a 2x2 frame
chrX <- mutate(chrX,chr="X Chromosome")
chr9 <- mutate(chr9,chr="Chromosome 9")

tyI_group <- rbind(chrX,chr9)

pdf("nullPvals_8pedfem_chr9X_ld05_ar1.pdf",width=12,height=12)
ggplot(tyI_group[tyI_group$model!="XKC_only",],aes(x=pval)) + geom_histogram(binwidth=1/10) + 
  facet_grid(chr~model) + xlab("p-value") + theme_bw() + ylab("Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(size=18), strip.text.x = element_text(size=18),
        strip.text.y = element_text(size=18))
dev.off()

rm(list=ls())


#####
# 8. Histogram of pvalues from 8ped null sims

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","method")
dat <- dat[!duplicated(dat$famSKAT),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  filter(famSKAT<=1) %>% # this filters out the Q stat rows
  gather(model,pval,-method) %>%
  mutate(pval=as.numeric(pval)) %>%
  group_by(model,method)

summarize(tyI_group,n()) # keats 9998

summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
summarize(tyI_group,tyI_err=sum(pval<0.005)/n())

chr9 <- tyI_group

pdf("nullSims_8ped_chr9_ld05_ar1/pvals_hist.pdf",width=12)
ggplot(tyI_group,aes(x=pval)) + geom_histogram(binwidth=1/10) + facet_wrap(~model)
dev.off()

##
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","method")
dat <- dat[!duplicated(dat$famSKAT),]

tyI_group <- dat %>% filter(!is.na(method)) %>%
  filter(famSKAT<=1) %>% # this filters out the Q stat rows
  gather(model,pval,-method) %>%
  mutate(pval=as.numeric(pval)) %>%
  group_by(model,method)

summarize(tyI_group,n()) # keats 9527

summarize(tyI_group,tyI_err=sum(pval<0.05)/n())
summarize(tyI_group,tyI_err=sum(pval<0.005)/n())  

pdf("nullSims_8ped_chrX_ld05_ar1/pvals_hist.pdf",width=12)
ggplot(tyI_group,aes(x=pval)) + geom_histogram(binwidth=1/10) + facet_wrap(~model)
dev.off()

chrX <- tyI_group

### 
# make pvalue histograms of just famSKAT & keats, and combine xchr & chr 9 results into a 2x2 frame
chrX <- mutate(chrX,chr="X Chromosome")
chr9 <- mutate(chr9,chr="Chromosome 9")

tyI_group <- rbind(chrX,chr9)

pdf("nullPvals_8ped_chr9X_ld05_ar1.pdf",width=12,height=12)
ggplot(tyI_group[tyI_group$model!="XKC_only",],aes(x=pval)) + geom_histogram(binwidth=1/10) + 
  facet_grid(chr~model) + xlab("p-value") + theme_bw() + ylab("Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(size=18), strip.text.x = element_text(size=18),
        strip.text.y = element_text(size=18))
dev.off()

rm(list=ls())


#####
# 9. Visualize LD structure

library(ggplot2)
library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/keats_x")

chr9_ceu <- getobj("LDMat_ceu_chr9_100snps.RData")
chr9_mxl <- getobj("LDMat_mxl_chr9_100snps.RData")

ar1 <- getobj("LDMat_05_ar1_100snps.RData")

# x chr is non-par FYI
x_ceu <- getobj("LDMat_ceu_chrX_100snps.RData")
x_mxl <- getobj("LDMat_mxl_chrX_100snps.RData")


cont = as.table(chr9_ceu)
df = as.data.frame(cont)
df = df[complete.cases(df),]

pdf("ld_chr9_ceu_heatmap.pdf")
ggplot(df,aes(Var1,Var2)) + geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")
dev.off()


cont = as.table(x_ceu)
df = as.data.frame(cont)
df = df[complete.cases(df),]

pdf("ld_chrX_ceu_heatmap.pdf")
ggplot(df,aes(Var1,Var2)) + geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")
dev.off()
##

cont = as.table(chr9_mxl)
df = as.data.frame(cont)
df = df[complete.cases(df),]

pdf("ld_chr9_mxl_heatmap.pdf")
ggplot(df,aes(Var1,Var2)) + geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")
dev.off()


cont = as.table(x_mxl)
df = as.data.frame(cont)
df = df[complete.cases(df),]

pdf("ld_chrX_mxl_heatmap.pdf")
ggplot(df,aes(Var1,Var2)) + geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")
dev.off()


cont = as.table(ar1)
df = as.data.frame(cont)
df = df[complete.cases(df),]

pdf("ld_ar1_05_heatmap.pdf")
ggplot(df,aes(Var1,Var2)) + geom_tile(aes(fill = Freq), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")
dev.off()



###

## look at all together, xaxis = dist between SNPs, y = ld; one line per chr type, hm population
# need to merge in position info
# get phased info for chr 9 for these samples
dat <- read.table(gzfile("/projects/geneva/gcc-fs2/1000Genomes/20130502/IMPUTE2/1000GP_Phase3/1000GP_Phase3_chr19.hap.gz"),
                  header=F,as.is=T,sep=" ",nrow=1000)
dim(dat) # 1000 5008

snp <- read.table(gzfile("/projects/geneva/gcc-fs2/1000Genomes/20130502/IMPUTE2/1000GP_Phase3/1000GP_Phase3_chr19.legend.gz"),
                  header=T,as.is=T,sep=" ",nrow=1000)
dim(snp) # 1000 11

position <- snp$position

samp <- read.table("/projects/geneva/gcc-fs2/1000Genomes/20130502/IMPUTE2/1000GP_Phase3/1000GP_Phase3.sample",
                   header=TRUE,as.is=TRUE,sep=" ")
dim(samp); head(samp) # 2504 4
table(samp$POP)

mxl <- samp[samp$POP=="MXL",]
dim(mxl) # 64 4
mxl$sex[mxl$SEX=="male"] <- "M"
mxl$sex[mxl$SEX=="female"] <- "F"

# find the cols of dat that correspond to the 64 MXL samples
datPop <- rep(samp$POP,each=2)
haps9 <- dat[,datPop=="MXL"]
dim(haps9) # 1000 128
haps9$position <- position

afreq <- apply(haps9[,1:128],1,function(x){sum(x)/128})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 188 128, so snp x sample (x2) + position
afreq <- afreq[afreq!=0]

haps9 <- haps9[-c(3,19),]
# haps9$position is position for MXL chr 9 position
colnames(chr9_mxl) <- rownames(chr9_mxl) <- haps9$position[1:100]

##
ceu <- samp[samp$POP=="CEU",]
dim(ceu) # 99 4
ceu$sex[ceu$SEX=="male"] <- "M"
ceu$sex[ceu$SEX=="female"] <- "F"

datPop <- rep(samp$POP,each=2)
haps9 <- dat[,datPop=="CEU"]
dim(haps9) # 1000 198
haps9 <- as.data.frame(haps9)
haps9$position <- position

afreq <- apply(haps9[,1:198],1,function(x){sum(x)/198})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 186 198, so snp x sample (x2)
afreq <- afreq[afreq!=0]

haps9 <- haps9[-c(27),]
afreq <- apply(haps9[,1:197],1,function(x){sum(x)/198})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 185 198, so snp x sample (x2)
# haps9 is position for CEU chr 9
colnames(chr9_ceu) <- rownames(chr9_ceu) <- haps9$position[1:100]


####
# now for x chr

# get phased info for X chr for these samples
dat <- read.table(gzfile("/projects/geneva/gcc-fs2/1000Genomes/20101123/IMPUTE/ALL_1000G_phase1integrated_v3_impute_mar2012/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz"),
                  header=F,as.is=T,sep=" ",nrow=1000)
dim(dat) # 1000 2184
dat[dat=="-"] <- NA
dat <- as.matrix(dat)
dat <- matrix(as.integer(dat),nrow=1000,ncol=2184)

# get x chr legends file
snp <- read.table(gzfile("/projects/geneva/gcc-fs2/1000Genomes/20101123/IMPUTE/ALL_1000G_phase1integrated_v3_impute_mar2012/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz"),
                  header=T,as.is=T,sep=" ",nrow=1000)
dim(snp) # 1000 12

position <- snp$position

samp <- read.table("/projects/geneva/gcc-fs2/1000Genomes/20101123/IMPUTE/ALL_1000G_phase1integrated_v3_impute_mar2012/ALL_1000G_phase1integrated_v3.sample",
                   header=TRUE,as.is=TRUE,sep=" ")
dim(samp); head(samp) # 2504 4
table(samp$population)

mxl <- samp[samp$population=="MXL",]
dim(mxl) # 66 4

# find the cols of dat that correspond to the 64 MXL samples
datPop <- rep(samp$population,each=2)
haps9 <- dat[,datPop=="MXL"]
dim(haps9) # 1000 132
haps9 <- as.data.frame(haps9)
haps9$position <- position

afreq <- apply(haps9[,1:132],1,function(x){sum(x,na.rm=T)/101})
haps9 <- haps9[afreq!=0,]
dim(haps9) # 456 132, so snp x sample (x2)

rownames(x_mxl) <- colnames(x_mxl) <- haps9$position[1:100]

##
# now get for ceu samples
ceu <- samp[samp$population=="CEU",]
dim(ceu) # 85 4

datPop <- rep(samp$population,each=2)
haps9 <- dat[,datPop=="CEU"]
dim(haps9) # 1000 170
haps9 <- as.data.frame(haps9)
haps9$position <- position

afreq <- apply(haps9[,1:170],1,function(x){sum(x,na.rm=T)/125}) # denom here is # haplotypes -- diff than 170*2 since males only have one
haps9 <- haps9[afreq!=0,]
dim(haps9) # 442 170, so snp x sample (x2)

rownames(x_ceu) <- colnames(x_ceu) <- haps9$position[1:100]

## now plot distance vs ld value
chr9_ceu = as.table(chr9_ceu)
chr9_ceu = as.data.frame(chr9_ceu)
chr9_ceu = chr9_ceu[complete.cases(chr9_ceu),]
chr9_ceu$dist <- as.numeric(levels(chr9_ceu$Var1)[as.integer(chr9_ceu$Var1)])-
  as.numeric(levels(chr9_ceu$Var2)[as.integer(chr9_ceu$Var2)])

tmp = as.table(chr9_mxl)
tmp = as.data.frame(tmp)
chr9_mxl = tmp[complete.cases(tmp),]
chr9_mxl$dist <- as.numeric(levels(chr9_mxl$Var1)[as.integer(chr9_mxl$Var1)])-
  as.numeric(levels(chr9_mxl$Var2)[as.integer(chr9_mxl$Var2)])

tmp = as.table(x_mxl)
tmp = as.data.frame(tmp)
x_mxl = tmp[complete.cases(tmp),]
x_mxl$dist <- as.numeric(levels(x_mxl$Var1)[as.integer(x_mxl$Var1)])-
  as.numeric(levels(x_mxl$Var2)[as.integer(x_mxl$Var2)])

tmp = as.table(x_ceu)
tmp = as.data.frame(tmp)
x_ceu = tmp[complete.cases(tmp),]
x_ceu$dist <- as.numeric(levels(x_ceu$Var1)[as.integer(x_ceu$Var1)])-
  as.numeric(levels(x_ceu$Var2)[as.integer(x_ceu$Var2)])

x_ceu$pop <- "CEU"
chr9_ceu$pop <- "CEU"
chr9_mxl$pop <- "MXL"
x_mxl$pop <- "MXL"
chr9_ceu$chr <- 9
chr9_mxl$chr <- 9
x_mxl$chr <- "X"
x_ceu$chr <- "X"

toPl <- rbind(x_ceu,x_mxl,chr9_ceu,chr9_mxl)
toPl <- toPl[toPl$dist!=0,]
ggplot(toPl,aes(x=-dist/1e6,y=Freq,color=pop))+geom_point()+facet_wrap(~chr,scales="free")+
  theme_bw()
dev.off()

rm(list=ls())


#####
# 10. More tyIerr sims

# called 8ped, chr9 MXL: nullSim_v3.R which creates 1M output when called w -t 1-10000
# called 8ped, chrX MXL: "            "                 "

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")

# chr 9 first
dat <- read_delim("nullSims_8ped_chr9mxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat
str(dat) # great, numeric where we want it

tyI_group <- dat %>% 
  filter(stat=="pval") %>% # this filters out the Q stat rows
  gather(model,pval,-stat) %>%
  mutate(pval=as.numeric(pval)) %>%
  group_by(model)

summarize(tyI_group,n()) # 1M each!

summarize(tyI_group,tyI_err=sum(pval<0.01)/n())
summarize(tyI_group,tyI_err=sum(pval<0.001)/n())
summarize(tyI_group,tyI_err=sum(pval<0.0001)/n())

chr9 <- tyI_group

pdf("nullSims_8ped_chr9mxl/pvals_hist.pdf",width=12)
ggplot(tyI_group,aes(x=pval)) + geom_histogram(binwidth=1/10) + facet_wrap(~model)
dev.off()

# X chr now
dat <- read_delim("nullSims_8ped_chrXmxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat
str(dat) # great, numeric where we want it

tyI_group <- dat %>% 
  filter(stat=="pval") %>% # this filters out the Q stat rows
  gather(model,pval,-stat) %>%
  mutate(pval=as.numeric(pval)) %>%
  group_by(model)

summarize(tyI_group,n()) # 1M each!

summarize(tyI_group,tyI_err=sum(pval<0.01)/n())
summarize(tyI_group,tyI_err=sum(pval<0.001)/n())
summarize(tyI_group,tyI_err=sum(pval<0.0001)/n())

pdf("nullSims_8ped_chrXmxl/pvals_hist.pdf",width=12)
ggplot(tyI_group,aes(x=pval)) + geom_histogram(binwidth=1/10) + facet_wrap(~model)
dev.off()

chrX <- tyI_group

### 
# make pvalue histograms of just famSKAT & keats, and combine xchr & chr 9 results into a 2x2 frame
chrX <- mutate(chrX,chr="X Chromosome")
chr9 <- mutate(chr9,chr="Chromosome 9")

tyI_group <- rbind(chrX,chr9)

pdf("nullPvals_8ped_chr9Xmxl.pdf",width=12,height=12)
ggplot(tyI_group[tyI_group$model!="XKC_only",],aes(x=pval)) + geom_histogram(binwidth=1/10) + 
  facet_grid(chr~model) + xlab("p-value") + theme_bw() + ylab("Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(size=18), strip.text.x = element_text(size=18),
        strip.text.y = element_text(size=18))
dev.off()

rm(list=ls())


#####
# 11. Combine all null simulations for tyI error table

setwd("/projects/geneva/geneva_sata/caitlin/keats_x")
library(ggplot2); library(tidyr)
library(dplyr); library(readr)
library(GWASTools)

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                               delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- filter(dat,!is.na(stat))

chr9_8ped_ar1 <- dat[!duplicated(dat$famSKAT),]
chr9_8ped_ar1 <- chr9_8ped_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8ped_chr9ceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
chr9_8ped_ceu <- dat[!duplicated(dat$famSKAT),]
chr9_8ped_ceu <- chr9_8ped_ceu %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8ped_chr9mxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
chr9_8ped_mxl <- dat[!duplicated(dat$famSKAT),]
chr9_8ped_mxl <- chr9_8ped_mxl %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="MXL")

dat <- read_delim("nullSims_8pedfem_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- filter(dat,!is.na(stat))

chr9_8pedfem_ar1 <- dat[!duplicated(dat$famSKAT),]
chr9_8pedfem_ar1 <- chr9_8pedfem_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8pedfem_chr9ceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
chr9_8pedfem_ceu <- dat[!duplicated(dat$famSKAT),]
chr9_8pedfem_ceu <- chr9_8pedfem_ceu %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8pedfem_chr9mxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
chr9_8pedfem_mxl <- dat[!duplicated(dat$famSKAT),]
chr9_8pedfem_mxl <- chr9_8pedfem_mxl %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="MXL")

dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- filter(dat,!is.na(stat))

x_8ped_ar1 <- dat[!duplicated(dat$famSKAT),]
x_8ped_ar1 <- x_8ped_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8ped_chrXceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
x_8ped_ceu <- dat[!duplicated(dat$famSKAT),]
x_8ped_ceu <- x_8ped_ceu %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8ped_chrXmxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
x_8ped_mxl <- dat[!duplicated(dat$famSKAT),]
x_8ped_mxl <- x_8ped_mxl %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="MXL")

dat <- read_delim("nullSims_8pedfem_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- filter(dat,!is.na(stat))
x_8pedfem_ar1 <- dat[!duplicated(dat$famSKAT),]
x_8pedfem_ar1 <- x_8pedfem_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8pedfem_chrXceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
x_8pedfem_ceu <- dat[!duplicated(dat$famSKAT),]
x_8pedfem_ceu <- x_8pedfem_ceu %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8pedfem_chrXmxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- dat %>% 
  mutate(burden_SKAT=NA) %>%
  mutate(burden_XKC_only=NA) %>%
  mutate(burden_KEATS=NA)
dat <- dat[,c(1,2,3,5,6,7,4)]
dat <- filter(dat,stat=="pval")
x_8pedfem_mxl <- dat[!duplicated(dat$famSKAT),]
x_8pedfem_mxl <- x_8pedfem_mxl %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="MXL")

dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- filter(dat,!is.na(stat))
x_unrel_ar1 <- dat[!duplicated(dat$famSKAT),]
x_unrel_ar1 <- x_unrel_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_unrel_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- filter(dat,!is.na(stat))
chr9_unrel_ar1 <- dat[!duplicated(dat$famSKAT),]
chr9_unrel_ar1 <- chr9_unrel_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5")

# merge all these together
tyI <- rbind(chr9_8ped_ar1,chr9_8ped_ceu,chr9_8ped_mxl,
                   chr9_8pedfem_ar1,chr9_8pedfem_ceu,chr9_8pedfem_mxl,
                   x_8ped_ar1,x_8ped_ceu,x_8ped_mxl,
                   x_8pedfem_ar1,x_8pedfem_ceu,x_8pedfem_mxl,
                   x_unrel_ar1,chr9_unrel_ar1)

tyI_group <- tyI %>% 
  gather(model,pval,-c(stat,chr,ped,ld)) %>%
  mutate(pval=as.numeric(pval)) %>%
  filter(model!="XKC_only") %>%
  group_by(chr,ped,ld,model)

n <- summarize(tyI_group,n()) # 1M each, only 10k of the AR1 iters
res <- summarize(tyI_group, 
          "alpha_0.01"=sum(pval<0.01)/n(),
          "alpha_0.001"=sum(pval<0.001)/n(),
          "alpha_0.0001"=sum(pval<0.0001)/n())

# want the proportions to be "wide" by model
alpha01 <- tyI_group %>% summarize("alpha_0.01"=sum(pval<0.01)/n()) %>%
  spread(model,alpha_0.01)
alpha001 <- tyI_group %>% summarize("alpha_0.01"=sum(pval<0.001)/n()) %>%
  spread(model,alpha_0.01)
alpha0001 <- tyI_group %>% summarize("alpha_0.01"=sum(pval<0.0001)/n()) %>%
  spread(model,alpha_0.01)

# save these as xtable
library(xtable)
xtable(res,digits=5)
print(xtable(cbind(alpha01[,-7],alpha001[,c(4:6,8)],alpha0001[,c(4:6,8)]),digits=c(1,1,1,1,3,3,3,3,4,4,4,4,5,5,5,5)),
      include.rownames=FALSE)

rm(list=ls())


#####
# 12. Process power results 

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chr9_ld05_ar1/allRes_252550_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_8pedfem_ar1_252550 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="25/25/50")

# merge in 8ped iterations
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/allRes_252550_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_8ped_ar1_252550 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="25/25/50")

# merge in unrel iterations
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chr9_ld05_ar1/allRes_252550_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_unrel_ar1_252550 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="25/25/50")

# merge in x chr
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chrX_ld05_ar1/allRes_252550_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_unrel_ar1_252550 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="25/25/50")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes_252550_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_8ped_ar1_252550 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="25/25/50")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chrX_ld05_ar1/allRes_252550_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_8pedfem_ar1_252550 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="25/25/50")

# merge in other prop iterations
# do 40/40/20 since 80/0/20 is BORING

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes_404020_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_8ped_ar1_404020 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="40/40/20")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chrX_ld05_ar1/allRes_404020_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)

dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_8pedfem_ar1_404020 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="40/40/20")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chrX_ld05_ar1/allRes_404020_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_unrel_ar1_404020 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="40/40/20")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chr9_ld05_ar1/allRes_404020_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_8pedfem_ar1_404020 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="40/40/20")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/allRes_404020_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_8ped_ar1_404020 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="40/40/20")


dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chr9_ld05_ar1/allRes_404020_c02.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_unrel_ar1_404020 <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="40/40/20")

pow <- rbind(chr9_unrel_ar1_252550,chr9_8ped_ar1_252550,chr9_8pedfem_ar1_252550,
             x_unrel_ar1_252550,x_8ped_ar1_252550,x_8pedfem_ar1_252550,
             x_8ped_ar1_404020,x_8pedfem_ar1_404020, x_unrel_ar1_404020,
             chr9_8pedfem_ar1_404020,chr9_8ped_ar1_404020,chr9_unrel_ar1_404020)

pow <- pow %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)

n=summarize(pow,n()) # 10K of each type
res <- summarize(pow,"alpha_0.01"=sum(pval<0.01)/n(),
             "alpha_0.001"=sum(pval<0.001)/n())

# make a graph of these results
pdf("power_252550_alpha01_10Kiters.pdf",width=11)
ggplot(data.frame(res),aes(x=ped,y=alpha_0.01)) + geom_point(aes(shape=model))+ facet_wrap(~chr) +
  theme_bw() + ggtitle("h2=0.5, +/-/0=25/25/50, 20 Variants, alpha=0.01")
dev.off()

pdf("power_252550_alpha01_10Kiters_trunc.pdf",width=11)
ggplot(data.frame(res[res$alpha_0.01>0.5,]),aes(x=ped,y=alpha_0.01)) + geom_point(aes(shape=model))+ facet_wrap(~chr) +
  theme_bw() + ggtitle("h2=0.5, +/-/0=25/25/50, 20 Variants, alpha=0.01")
dev.off()

pdf("power_252550_alpha001_10Kiters.pdf",width=11)
ggplot(data.frame(res),aes(x=ped,y=alpha_0.001)) + geom_point(aes(shape=model))+ facet_wrap(~chr) +
  theme_bw() + ggtitle("h2=0.5, +/-/0=25/25/50, 20 Variants, alpha=0.001")
dev.off()

pdf("power_252550_alpha001_10Kiters_trunc.pdf",width=11)
ggplot(data.frame(res[res$alpha_0.001>0.5,]),aes(x=ped,y=alpha_0.001)) + geom_point(aes(shape=model))+ facet_wrap(~chr) +
  theme_bw() + ggtitle("h2=0.5, +/-/0=25/25/50, 20 Variants, alpha=0.001")
dev.off()

## we're showing lower power since tyIerror is off for famSKAT
# need to do the power graphs we do usually
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_unrel_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_unrel_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_8ped_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8pedfem_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

chr9_8pedfem_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8pedfem_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_8pedfem_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_8ped_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_unrel_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]

x_unrel_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

nullSims <- rbind(x_unrel_null,x_8ped_null,x_8pedfem_null,
                  chr9_unrel_null,chr9_8pedfem_null,chr9_8ped_null)

nullSims <- nullSims %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)

summarize(nullSims,n())

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=72,ncol=length(alpha))
falsePos <- matrix(NA,nrow=36,ncol=length(alpha))
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pval<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,5]
  
  s <- summarize(pow,sum(pval<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,5]
}

# make a plot of these two values
s <- summarize(pow,n())
truePos <- cbind(data.frame(s)[,1:4],truePos)
s <- summarize(nullSims,n())
falsePos <- cbind(data.frame(s)[,1:4],falsePos)

# first plot all autosomal results, all peds, prop=25/25/50
smfalse <- falsePos[falsePos$chr=="autosomal",]
smfalse$type = "null"
smtrue <- truePos[truePos$chr=="autosomal"&truePos$prop=="25/25/50",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="25/25/50"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_autoSNP_true252550_trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=25/25/50")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_autoSNP_true252550.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=25/25/50")
dev.off()

# x chr results
smfalse <- falsePos[falsePos$chr=="X",]
smfalse$type = "null"
smtrue <- truePos[truePos$chr=="X"&truePos$prop=="25/25/50",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="25/25/50"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_xSNP_true252550_trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=25/25/50")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_xSNP_true252550.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=25/25/50")
dev.off()

rm(list=ls())


#####
# 13. Process power results 40/40/20, 50/0/50

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_unrel_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
chr9_unrel_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
chr9_8ped_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8pedfem_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
chr9_8pedfem_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8pedfem_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
x_8pedfem_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
x_8ped_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_unrel_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
x_unrel_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

nullSims <- rbind(x_unrel_null,x_8ped_null,x_8pedfem_null,
                  chr9_unrel_null,chr9_8pedfem_null,chr9_8ped_null)
nullSims <- nullSims %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)
summarize(nullSims,n())

# read in power sims for 40/40/20 now

readPow <- function(fn){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  dat <- mutate(dat,X1=NULL)
  colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
  dat <- dat[!duplicated(dat$famSKAT),]
  dat <- dat[dat$famSKAT!="x.res",]
  tmp <- regexpr("_chr",fn)
  chr <- substr(fn,start=tmp[1]+4,stop=tmp[1]+4)
  tmp <- regexpr("powerSims_",fn)
  tmp2 <- regexpr("_chr",fn)
  ped <- substr(fn,start=tmp[1]+10,stop=tmp2[1]-4)
  tmp <- regexpr("chr",fn)
  tmp2 <- regexpr("/allRes",fn)
  ld <- substr(fn,start=tmp[1]+5,stop=tmp2[1]-1)
  tmp <- regexpr("_c0",fn)
  prop <- substr(fn,start=tmp2[1]+8,stop=tmp[1]-1)
  x <- dat %>% 
    mutate(famSKAT=as.numeric(famSKAT)) %>%
    mutate(XKC_only=as.numeric(XKC_only)) %>%
    mutate(KEATS_X=as.numeric(KEATS_X)) %>%
    mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
    mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
    mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
    filter(stat=="pval") %>%
    mutate(chr=chr) %>%
    mutate(ped=ped) %>%
    mutate(ld=ld) %>%
    mutate(prop=prop)
  
  return(x)
}

x_8ped_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes_404020_c02.txt")
x_8pedfem_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chrX_ld05_ar1/allRes_404020_c02.txt")
x_unrel_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chrX_ld05_ar1/allRes_404020_c02.txt")
chr9_8pedfem_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chr9_ld05_ar1/allRes_404020_c02.txt")
chr9_8ped_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/allRes_404020_c02.txt")
chr9_unrel_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chr9_ld05_ar1/allRes_404020_c02.txt")

pow <- rbind(x_8ped_ar1_404020,x_8pedfem_ar1_404020, x_unrel_ar1_404020,
             chr9_8pedfem_ar1_404020,chr9_8ped_ar1_404020,chr9_unrel_ar1_404020)

pow$ped[pow$ped=="8ped"] <- "8pedfem"
pow$ped[pow$ped=="8"] <- "8ped"
pow$ped[pow$ped=="un"] <- "unrel"
pow$chr[pow$chr=="9"] <- "autosomal"
pow <- pow %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)


alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=36,ncol=length(alpha))
falsePos <- matrix(NA,nrow=36,ncol=length(alpha))
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pval<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,5]
  
  s <- summarize(pow,sum(pval<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,5]
}

# make a plot of these two values
s <- summarize(pow,n())
truePos <- cbind(data.frame(s)[,1:4],truePos)
s <- summarize(nullSims,n())
falsePos <- cbind(data.frame(s)[,1:4],falsePos)

# first plot all autosomal results, all peds, prop=404020
smfalse <- falsePos[falsePos$chr=="autosomal",]
smfalse$type = "null"

smtrue <- truePos[truePos$chr=="autosomal"&truePos$prop=="404020",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="404020"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_autoSNP_true404020trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=40/40/20")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_autoSNP_true404020.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=40/40/20")
dev.off()

# x chr results
smfalse <- falsePos[falsePos$chr=="X",]
smfalse$type = "null"
smtrue <- truePos[truePos$chr=="X"&truePos$prop=="404020",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="404020"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_xSNP_true404020_trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=40/40/20")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_xSNP_true404020.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=40/40/20")
dev.off()

## read in 50/0/50 results

x_8ped_ar1_50050 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes_50050_c02.txt")
x_8pedfem_ar1_50050 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chrX_ld05_ar1/allRes_50050_c02.txt")
x_unrel_ar1_50050 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chrX_ld05_ar1/allRes_50050_c02.txt")
chr9_8pedfem_ar1_50050 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chr9_ld05_ar1/allRes_50050_c02.txt")
chr9_8ped_ar1_50050 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/allRes_50050_c02.txt")
chr9_unrel_ar1_50050 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chr9_ld05_ar1/allRes_50050_c02.txt")

pow <- rbind(x_8ped_ar1_50050,x_8pedfem_ar1_50050, x_unrel_ar1_50050,
             chr9_8pedfem_ar1_50050,chr9_8ped_ar1_50050,chr9_unrel_ar1_50050)

pow$ped[pow$ped=="8ped"] <- "8pedfem"
pow$ped[pow$ped=="8"] <- "8ped"
pow$ped[pow$ped=="un"] <- "unrel"
pow$chr[pow$chr=="9"] <- "autosomal"
pow <- pow %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)


alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=36,ncol=length(alpha))
falsePos <- matrix(NA,nrow=36,ncol=length(alpha))
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pval<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,5]
  
  s <- summarize(pow,sum(pval<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,5]
}

# make a plot of these two values
s <- summarize(pow,n())
truePos <- cbind(data.frame(s)[,1:4],truePos)
s <- summarize(nullSims,n())
falsePos <- cbind(data.frame(s)[,1:4],falsePos)

# first plot all autosomal results, all peds, prop=50050
smfalse <- falsePos[falsePos$chr=="autosomal",]
smfalse$type = "null"

smtrue <- truePos[truePos$chr=="autosomal"&truePos$prop=="50050",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="50050"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_autoSNP_true50050trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=50/0/50")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_autoSNP_true50050.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=50/0/50")
dev.off()

# x chr results
smfalse <- falsePos[falsePos$chr=="X",]
smfalse$type = "null"
smtrue <- truePos[truePos$chr=="X"&truePos$prop=="50050",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="50050"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_xSNP_true50050_trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=50/0/50")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_xSNP_true50050.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=50/0/50")
dev.off()

rm(list=ls())


#####
# 14. Process power results 40/40/20; diff sigmaA, sigmaX values

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_unrel_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
chr9_unrel_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
chr9_8ped_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8pedfem_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
chr9_8pedfem_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8pedfem_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
x_8pedfem_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
x_8ped_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

dat <- read_delim("/projects/geneva/geneva_sata/caitlin/keats_x/nullSims_unrel_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
dat <- dat[!duplicated(dat$famSKAT),]
dat <- dat[dat$famSKAT!="x.res",]
x_unrel_null <- dat %>% 
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
  mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
  mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
  filter(stat=="pval") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel") %>%
  mutate(ld="AR1_0.5") %>%
  mutate(prop="null")

nullSims <- rbind(x_unrel_null,x_8ped_null,x_8pedfem_null,
                  chr9_unrel_null,chr9_8pedfem_null,chr9_8ped_null)
nullSims <- nullSims %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)
summarize(nullSims,n())

# read in power sims for 40/40/20, h06 

readPow <- function(fn){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  dat <- mutate(dat,X1=NULL)
  colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","burden_SKAT","burden_XKC_only","burden_KEATS","stat")
  dat <- dat[!duplicated(dat$famSKAT),]
  dat <- dat[dat$famSKAT!="x.res",]
  tmp <- regexpr("_chr",fn)
  chr <- substr(fn,start=tmp[1]+4,stop=tmp[1]+4)
  tmp <- regexpr("powerSims_",fn)
  tmp2 <- regexpr("_chr",fn)
  ped <- substr(fn,start=tmp[1]+10,stop=tmp2[1]-4)
  tmp <- regexpr("chr",fn)
  tmp2 <- regexpr("/allRes",fn)
  ld <- substr(fn,start=tmp[1]+5,stop=tmp2[1]-1)
  tmp <- regexpr("_c0",fn)
  prop <- substr(fn,start=tmp2[1]+8,stop=tmp[1]-1)
  x <- dat %>% 
    mutate(famSKAT=as.numeric(famSKAT)) %>%
    mutate(XKC_only=as.numeric(XKC_only)) %>%
    mutate(KEATS_X=as.numeric(KEATS_X)) %>%
    mutate(burden_SKAT=as.numeric(burden_SKAT)) %>%
    mutate(burden_XKC_only=as.numeric(burden_XKC_only)) %>%
    mutate(burden_KEATS=as.numeric(burden_KEATS)) %>%
    filter(stat=="pval") %>%
    mutate(chr=chr) %>%
    mutate(ped=ped) %>%
    mutate(ld=ld) %>%
    mutate(prop=prop)
  
  return(x)
}

x_8ped_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chrX_ld05_ar1/allRes_404020_c02_h06.txt")
x_8pedfem_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chrX_ld05_ar1/allRes_404020_c02_h06.txt")
x_unrel_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chrX_ld05_ar1/allRes_404020_c02_h06.txt")
chr9_8pedfem_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8pedfem_chr9_ld05_ar1/allRes_404020_c02_h06.txt")
chr9_8ped_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_8ped_chr9_ld05_ar1/allRes_404020_c02_h06.txt")
chr9_unrel_ar1_404020 <- readPow("/projects/geneva/geneva_sata/caitlin/keats_x/powerSims_unrel_chr9_ld05_ar1/allRes_404020_c02_h06.txt")

pow <- rbind(x_8ped_ar1_404020,x_8pedfem_ar1_404020, x_unrel_ar1_404020,
             chr9_8pedfem_ar1_404020,chr9_8ped_ar1_404020,chr9_unrel_ar1_404020)

pow$ped[pow$ped=="8ped"] <- "8pedfem"
pow$ped[pow$ped=="8"] <- "8ped"
pow$ped[pow$ped=="un"] <- "unrel"
pow$chr[pow$chr=="9"] <- "autosomal"
pow <- pow %>%
  gather(model,pval,-c(stat,chr,ped,ld,prop)) %>%
  group_by(model,chr,ped,prop)


alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=36,ncol=length(alpha))
falsePos <- matrix(NA,nrow=36,ncol=length(alpha))
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pval<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,5]
  
  s <- summarize(pow,sum(pval<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,5]
}

# make a plot of these two values
s <- summarize(pow,n())
truePos <- cbind(data.frame(s)[,1:4],truePos)
s <- summarize(nullSims,n())
falsePos <- cbind(data.frame(s)[,1:4],falsePos)

# first plot all autosomal results, all peds, prop=404020
smfalse <- falsePos[falsePos$chr=="autosomal",]
smfalse$type = "null"

smtrue <- truePos[truePos$chr=="autosomal"&truePos$prop=="404020",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="404020"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_autoSNP_true404020_h06_trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=40/40/20")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_autoSNP_true404020_h06.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 Autosomal Variants, +/-/0=40/40/20")
dev.off()

# x chr results
smfalse <- falsePos[falsePos$chr=="X",]
smfalse$type = "null"
smtrue <- truePos[truePos$chr=="X"&truePos$prop=="404020",]
smtrue$type = "true"
allres <- rbind(smfalse,smtrue)

allres <- allres %>%
  gather(iter,rate,-c(model,chr,ped,prop))
allres$prop[allres$prop=="404020"]<- "true"
allres$rate <- as.numeric(allres$rate)

toPl <- allres %>%
  spread(prop,rate) 

pdf("power_xSNP_true404020_h06_trunc.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=40/40/20")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

pdf("power_xSNP_true404020_h06.pdf",width=20)
ggplot(toPl,aes(x=null,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + 
  facet_wrap(~ped) + ggtitle("20 X Chromosome Variants, +/-/0=40/40/20")
dev.off()

rm(list=ls())


#####
# 15. Process power results, keatso x chr, 8ped

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
  mutate(qstat=as.numeric(qstat)) %>%
  mutate(pvalue=as.numeric(pvalue)) %>%
  mutate(rho=as.numeric(rho)) %>%
  mutate(iter=rep(1:n,each=12*4)) %>%
  mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
  filter(finalpval==TRUE)
  return(datNew)
}

chrX_8ped_ar1_208 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim208_c02_allOptRes.txt")
chrX_8ped_ar1_406 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim406_c02_allOptRes.txt")
chrX_8ped_ar1_604 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim604_c02_allOptRes.txt")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim226_c02_allOptRes.txt",n=9892)

chrX_8ped_ar1_208 <- chrX_8ped_ar1_208 %>%
  mutate(prop="20/0/80")
chrX_8ped_ar1_406 <- chrX_8ped_ar1_406 %>%
  mutate(prop="40/0/60")
chrX_8ped_ar1_604 <- chrX_8ped_ar1_604 %>%
  mutate(prop="60/0/40")
chrX_8ped_ar1_226 <- chrX_8ped_ar1_226 %>%
  mutate(prop="20/20/60")

allPow <- rbind(chrX_8ped_ar1_208,chrX_8ped_ar1_406,
                chrX_8ped_ar1_604,chrX_8ped_ar1_226)

# want to calculate the power at alpha=1e-4
pow <- allPow %>%
  group_by(model,prop)
summarize(pow,n()) # 10K for each model at each prop

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  group_by(model)

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=16,ncol=length(alpha))
falsePos <- matrix(NA,nrow=4,ncol=length(alpha))
colnames(truePos) <- paste0("alpha.",alpha)
colnames(falsePos) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pvalue<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,2]
  rownames(falsePos) <- data.frame(s)[,1]
  
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,3]
  rownames(truePos) <- paste(data.frame(s)[,1],data.frame(s)[,2],sep=".")
}

# make a 2x2 plot, one for each proportion
# plot each model in the panel

falsePos <- data.frame(falsePos)
falsePos$type <- "false"
rnms <- paste(c("auto","both","unrel","x"),"20/0/80",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"

smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
falsePos$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

allres <- rbind(smtrue,falsePos)

allres <- allres %>%
  gather(alpha,rate,-c(model,type))

toPl <- allres %>%
  spread(type,rate) 

# plot these results now
pdf("power_xSNP_trunc.pdf")
par(mfrow=c(2,2))
ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=20/0/80")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"40/0/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=40/0/60")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"60/0/40",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=60/0/40")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"20/20/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=20/20/60")+ scale_x_continuous(limits=c(0,0.1))

dev.off()


#####
# now do the same thing for autosomal SNPs

chrX_8ped_ar1_208 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim208_c02_allOptRes.txt")
chrX_8ped_ar1_406 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim406_c02_allOptRes.txt")
chrX_8ped_ar1_604 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim604_c02_allOptRes.txt")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim226_c02_allOptRes.txt")

chrX_8ped_ar1_208 <- chrX_8ped_ar1_208 %>%
  mutate(prop="20/0/80")
chrX_8ped_ar1_406 <- chrX_8ped_ar1_406 %>%
  mutate(prop="40/0/60")
chrX_8ped_ar1_604 <- chrX_8ped_ar1_604 %>%
  mutate(prop="60/0/40")
chrX_8ped_ar1_226 <- chrX_8ped_ar1_226 %>%
  mutate(prop="20/20/60")

allPow <- rbind(chrX_8ped_ar1_208,chrX_8ped_ar1_406,
                chrX_8ped_ar1_604,chrX_8ped_ar1_226)

# want to calculate the power at alpha=1e-4
pow <- allPow %>%
  group_by(model,prop)
summarize(pow,n()) # 10K for each model at each prop

# now read in null simulations
dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  group_by(model)

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=16,ncol=length(alpha))
falsePos <- matrix(NA,nrow=4,ncol=length(alpha))
colnames(truePos) <- paste0("alpha.",alpha)
colnames(falsePos) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pvalue<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,2]
  rownames(falsePos) <- data.frame(s)[,1]
  
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,3]
  rownames(truePos) <- paste(data.frame(s)[,1],data.frame(s)[,2],sep=".")
}

# make a 2x2 plot, one for each proportion
# plot each model in the panel

falsePos <- data.frame(falsePos)
falsePos$type <- "false"
rnms <- paste(c("auto","both","unrel","x"),"20/0/80",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"

smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
falsePos$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

allres <- rbind(smtrue,falsePos)

allres <- allres %>%
  gather(alpha,rate,-c(model,type))

toPl <- allres %>%
  spread(type,rate) 

# plot these results now
pdf("power_autoSNP_trunc.pdf")
par(mfrow=c(2,2))
ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Chr Variants, +/-/0=20/0/80")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"40/0/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Chr Variants, +/-/0=40/0/60")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"60/0/40",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Chr Variants, +/-/0=60/0/40")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"20/20/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Chr Variants, +/-/0=20/20/60")+ scale_x_continuous(limits=c(0,0.1))

dev.off()


#####
# do with the extreme sigma values

chrX_8ped_ar1_208 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim208_x3a03_allOptRes.txt")
chrX_8ped_ar1_406 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim406_x3a03_allOptRes.txt")
chrX_8ped_ar1_604 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim604_x3a03_allOptRes.txt")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim226_x3a03_allOptRes.txt")

chrX_8ped_ar1_208 <- chrX_8ped_ar1_208 %>%
  mutate(prop="20/0/80")
chrX_8ped_ar1_406 <- chrX_8ped_ar1_406 %>%
  mutate(prop="40/0/60")
chrX_8ped_ar1_604 <- chrX_8ped_ar1_604 %>%
  mutate(prop="60/0/40")
chrX_8ped_ar1_226 <- chrX_8ped_ar1_226 %>%
  mutate(prop="20/20/60")

allPow <- rbind(chrX_8ped_ar1_208,chrX_8ped_ar1_406,
                chrX_8ped_ar1_604,chrX_8ped_ar1_226)

# want to calculate the power at alpha=1e-4
pow <- allPow %>%
  group_by(model,prop)
summarize(pow,n()) # 10K for each model at each prop

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/sim_x3a03_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  group_by(model)

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=16,ncol=length(alpha))
falsePos <- matrix(NA,nrow=4,ncol=length(alpha))
colnames(truePos) <- paste0("alpha.",alpha)
colnames(falsePos) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pvalue<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,2]
  rownames(falsePos) <- data.frame(s)[,1]
  
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,3]
  rownames(truePos) <- paste(data.frame(s)[,1],data.frame(s)[,2],sep=".")
}

falsePos <- data.frame(falsePos)
falsePos$type <- "false"
rnms <- paste(c("auto","both","unrel","x"),"20/0/80",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"

smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
falsePos$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

allres <- rbind(smtrue,falsePos)

allres <- allres %>%
  gather(alpha,rate,-c(model,type))

toPl <- allres %>%
  spread(type,rate) 

# plot these results now
pdf("power_xSNP_x3a03_trunc.pdf")
par(mfrow=c(2,2))
ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.5,1)) +
  ggtitle("20 X Chr Variants, +/-/0=20/0/80")+ scale_x_continuous(limits=c(0,0.5))

rnms <- paste(c("auto","both","unrel","x"),"40/0/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.5,1)) +
  ggtitle("20 X Chr Variants, +/-/0=40/0/60")+ scale_x_continuous(limits=c(0,0.5))

rnms <- paste(c("auto","both","unrel","x"),"60/0/40",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.5,1)) +
  ggtitle("20 X Chr Variants, +/-/0=60/0/40")+ scale_x_continuous(limits=c(0,0.5))

rnms <- paste(c("auto","both","unrel","x"),"20/20/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.5,1)) +
  ggtitle("20 X Chr Variants, +/-/0=20/20/60")+ scale_x_continuous(limits=c(0,0.5))

dev.off()


#####
# 16. Power results, keatso, chr9 and auto for unrelated samples

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE)
  return(datNew)
}

chrX_8ped_ar1_226 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim226_h00_allOptRes.txt")
chrX_8ped_ar1_226 <- chrX_8ped_ar1_226 %>%
  mutate(prop="20/20/60")

allPow <- chrX_8ped_ar1_226

# want to calculate the power at alpha=1e-4
pow <- allPow %>%
  group_by(model,prop)
summarize(pow,n()) # 10K for each model at each prop

# now read in null simulations
dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/h00_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  group_by(model)

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=4,ncol=length(alpha))
falsePos <- matrix(NA,nrow=4,ncol=length(alpha))
colnames(truePos) <- paste0("alpha.",alpha)
colnames(falsePos) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pvalue<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,2]
  rownames(falsePos) <- data.frame(s)[,1]
  
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,3]
  rownames(truePos) <- paste(data.frame(s)[,1],data.frame(s)[,2],sep=".")
}

# make a 2x2 plot, one for each proportion
# plot each model in the panel

falsePos <- data.frame(falsePos)
falsePos$type <- "false"
falsePos$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

smtrue <- data.frame(truePos)
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))

toPl <- allres %>%
  spread(type,rate) 

# plot these results now
pdf("power_xSNP_unrel_h00_trunc.pdf")
ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=20/20/60")+ scale_x_continuous(limits=c(0,0.1))
dev.off()

rm(list=ls())


#####
# 17. Power results, keatso, chrX and auto for different mixes of +/-

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim15256_c02_allOptRes.txt")
chrX_8ped_ar1_028 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim028_c02_allOptRes.txt")
chrX_8ped_ar1_064 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim064_c02_allOptRes.txt")

chrX_8ped_ar1_15256 <- chrX_8ped_ar1_15256 %>%
  mutate(prop="15/25/60")
chrX_8ped_ar1_028 <- chrX_8ped_ar1_028 %>%
  mutate(prop="0/20/80")
chrX_8ped_ar1_064 <- chrX_8ped_ar1_064 %>%
  mutate(prop="0/60/40")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_028,
                chrX_8ped_ar1_064)

# want to calculate the power at alpha=1e-4
pow <- allPow %>%
  group_by(model,prop)
summarize(pow,n()) # 10K for each model at each prop

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  group_by(model)

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=12,ncol=length(alpha))
falsePos <- matrix(NA,nrow=4,ncol=length(alpha))
colnames(truePos) <- paste0("alpha.",alpha)
colnames(falsePos) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pvalue<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,2]
  rownames(falsePos) <- data.frame(s)[,1]
  
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,3]
  rownames(truePos) <- paste(data.frame(s)[,1],data.frame(s)[,2],sep=".")
}

falsePos <- data.frame(falsePos)
falsePos$type <- "false"
rnms <- paste(c("auto","both","unrel","x"),"15/25/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"

smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
falsePos$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

allres <- rbind(smtrue,falsePos)

allres <- allres %>%
  gather(alpha,rate,-c(model,type))

toPl <- allres %>%
  spread(type,rate) 

# plot these results now
pdf("power_xSNP_moreProps_trunc.pdf")
par(mfrow=c(2,2))
ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=15/25/60")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"0/20/80",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=0/20/80")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"0/60/40",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 X Chr Variants, +/-/0=0/60/40")+ scale_x_continuous(limits=c(0,0.1))

dev.off()

## now for autosomal

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim15256_c02_allOptRes.txt")
chrX_8ped_ar1_028 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim028_c02_allOptRes.txt")
chrX_8ped_ar1_064 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim064_c02_allOptRes.txt")

chrX_8ped_ar1_15256 <- chrX_8ped_ar1_15256 %>%
  mutate(prop="15/25/60")
chrX_8ped_ar1_028 <- chrX_8ped_ar1_028 %>%
  mutate(prop="0/20/80")
chrX_8ped_ar1_064 <- chrX_8ped_ar1_064 %>%
  mutate(prop="0/60/40")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_028,
                chrX_8ped_ar1_064)

# want to calculate the power at alpha=1e-4
pow <- allPow %>%
  group_by(model,prop)
summarize(pow,n()) # 10K for each model at each prop

# now read in null simulations
dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  group_by(model)

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
truePos <- matrix(NA,nrow=12,ncol=length(alpha))
falsePos <- matrix(NA,nrow=4,ncol=length(alpha))
colnames(truePos) <- paste0("alpha.",alpha)
colnames(falsePos) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(nullSims,sum(pvalue<alpha[i])/n())
  falsePos[,i] <- data.frame(s)[,2]
  rownames(falsePos) <- data.frame(s)[,1]
  
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  truePos[,i] <- data.frame(s)[,3]
  rownames(truePos) <- paste(data.frame(s)[,1],data.frame(s)[,2],sep=".")
}

falsePos <- data.frame(falsePos)
falsePos$type <- "false"
rnms <- paste(c("auto","both","unrel","x"),"15/25/60",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"

smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
falsePos$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")

allres <- rbind(smtrue,falsePos)

allres <- allres %>%
  gather(alpha,rate,-c(model,type))

toPl <- allres %>%
  spread(type,rate) 

# plot these results now
pdf("power_autoSNP_moreProps_trunc.pdf")
ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Variants, +/-/0=15/25/60")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"0/20/80",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Variants, +/-/0=0/20/80")+ scale_x_continuous(limits=c(0,0.1))

rnms <- paste(c("auto","both","unrel","x"),"0/60/40",sep=".")
smtrue <- data.frame(truePos[is.element(rownames(truePos),rnms),])
smtrue$type <- "true"
smtrue$model <- c("MONSTER","KEATS-O","SKATO","KEATS-O_xOnly")
allres <- rbind(smtrue,falsePos)
allres <- allres %>%
  gather(alpha,rate,-c(model,type))
toPl <- allres %>%
  spread(type,rate) 

ggplot(toPl,aes(x=false,y=true,color=model)) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  ggtitle("20 Autosomal Variants, +/-/0=0/60/40")+ scale_x_continuous(limits=c(0,0.1))

dev.off()

rm(list=ls())


#####
# 18. Power results, keatso, chrX and auto for different mixes of +/-

# a 2x3 plot, auto vs x on the rows, cols corresponding to props of:
# 5/35/60, 15/25/60, 20/20/60

# need to put null and true sim results together,
# indicators for null/true, x/auto and proportion, as well as model

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim226_c02_allOptRes.txt",n=9892,prop="20/20/60",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
summarize(pow,n()) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKAT"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_8ped_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
 scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

# make histograms of the rho values at the same configurations
powTr <- pow[pow$type==TRUE&pow$model=="both",]
powTr$prop <- ordered(powTr$prop,levels=c("5/35/60","15/25/60","20/20/60"))
pdf("rho_ped8_keatso_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("")
dev.off()

powTr <- pow[pow$type==TRUE&pow$model=="auto",]
powTr$prop <- ordered(powTr$prop,levels=c("5/35/60","15/25/60","20/20/60"))
pdf("rho_ped8_monster_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("")
dev.off()

rm(list=ls())


#####
# 19. Compare keatso to keats.bt and keatso rho=0 (keats.skat)

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr,bt=NULL){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  if(is.null(bt)){
    colnames(dat) <- c("qstat","pvalue","rho","model")
    dat <- dat[dat$model!="model",]
    datNew <- dat %>%
      mutate(qstat=as.numeric(qstat)) %>%
      mutate(pvalue=as.numeric(pvalue)) %>%
      mutate(rho=as.numeric(rho)) %>%
      #mutate(iter=rep(1:n,each=12*4)) %>%
      mutate(finalpval=rep(c(TRUE,rep(FALSE,10),TRUE),4*n)) %>%
      filter(finalpval==TRUE) %>%
      filter(model=="both") 
  }else{
    colnames(dat) <- c("auto","x","both","unrel","stat")
    dat <- dat[dat$auto!="auto.burden",]
    datNew <- dat %>%
      filter(stat=="pval") %>%
      mutate(both=as.numeric(both)) %>%
      mutate(auto=NULL) %>%
      mutate(x=NULL) %>%
      mutate(unrel=NULL)
  }
  datNew <- datNew %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim226_c02_allOptRes.txt",n=9892,prop="20/20/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_406 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_604 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                             type=TRUE,chr="X")
chrX_8ped_ar1_208 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="X")

# same for autosomal
ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="Autosomal")
ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="Autosomal")
ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")
ped_ar1_406 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                             type=TRUE,chr="Autosomal")
ped_ar1_604 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                             type=TRUE,chr="Autosomal")
ped_ar1_208 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chrX_8ped_ar1_406,chrX_8ped_ar1_604,chrX_8ped_ar1_208,
                ped_ar1_15256,ped_ar1_5356,ped_ar1_226,
                ped_ar1_406,ped_ar1_604,ped_ar1_208)
allPow$model <- rep(c("KEATS","KEATS-O"),nrow(allPow)/2)
allPow$qstat <- NULL
allPow$rho <- NULL; allPow$finalpval <- NULL

# need to get keats.bt results now
chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim15256_c02_allBtRes.txt",prop="15/25/60",
                               type=TRUE,chr="X",bt=TRUE)
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim5356_c02_allBtRes.txt",prop="5/35/60",
                              type=TRUE,chr="X",bt=TRUE)
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim226_c02_allBtRes.txt",n=9892,prop="20/20/60",
                             type=TRUE,chr="X",bt=TRUE)
chrX_8ped_ar1_406 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim406_c02_allBtRes.txt",prop="40/0/60",
                             type=TRUE,chr="X",bt=TRUE)
chrX_8ped_ar1_604 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim604_c02_allBtRes.txt",prop="60/0/40",
                             type=TRUE,chr="X",bt=TRUE)
chrX_8ped_ar1_208 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim208_c02_allBtRes.txt",prop="20/0/80",
                             type=TRUE,chr="X",bt=TRUE)

# same for autosomal
ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim15256_c02_allBtRes.txt",prop="15/25/60",
                         type=TRUE,chr="Autosomal",bt=TRUE)
ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim5356_c02_allBtRes.txt",prop="5/35/60",
                        type=TRUE,chr="Autosomal",bt=TRUE)
ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim226_c02_allBtRes.txt",prop="20/20/60",
                       type=TRUE,chr="Autosomal",bt=TRUE)
ped_ar1_406 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim406_c02_allBtRes.txt",prop="40/0/60",
                       type=TRUE,chr="Autosomal",bt=TRUE)
ped_ar1_604 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim604_c02_allBtRes.txt",prop="60/0/40",
                       type=TRUE,chr="Autosomal",bt=TRUE)
ped_ar1_208 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim208_c02_allBtRes.txt",prop="20/0/80",
                       type=TRUE,chr="Autosomal",bt=TRUE)

allBt <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
               chrX_8ped_ar1_406,chrX_8ped_ar1_604,chrX_8ped_ar1_208,
               ped_ar1_15256,ped_ar1_5356,ped_ar1_226,
               ped_ar1_406,ped_ar1_604,ped_ar1_208)
allBt$stat <- "KEATS-bt"
colnames(allBt)[1:2] <- c("pvalue","model")

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model=="both") %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(pvalue=as.numeric(pvalue)) %>%
  mutate(chr="X") 
nullSims$qstat <- NULL; nullSims$rho <- NULL
nullSims$model <- "KEATS-O"

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"
nullSimsc <- nullSimsa
nullSimsc$prop <- "40/0/60"
nullSimsd <- nullSimsa
nullSimsd$prop <- "60/0/40"
nullSimse <- nullSimsa
nullSimse$prop <- "20/0/80"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model=="both") %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") %>%
  mutate(pvalue=as.numeric(pvalue))
nullSims2$qstat <- NULL; nullSims2$rho <- NULL
nullSims2$model <- "KEATS-O"

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"
nullSims5 <- nullSims2
nullSims5$prop <- "40/0/60"
nullSims6 <- nullSims2
nullSims6$prop <- "60/0/40"
nullSims7 <- nullSims2
nullSims7$prop <- "20/0/80"

allRes <- rbind(allPow,allBt,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4,nullSimsc,nullSimsd,nullSimse,
                nullSims5,nullSims6,nullSims7)

## now read in null.bt
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allBtRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("auto","x","both","unrel","stat")
dat <- dat[dat$auto!="auto.burden",]
dat$auto <- NULL; dat$x <- NULL; dat$unrel <- NULL
nullSims <- dat %>%
  filter(stat=="pval") %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") %>%
  mutate(both=as.numeric(both))
nullSims$stat <- "KEATS-bt"
colnames(nullSims)[1:2] <- c("pvalue","model")

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"
nullSimsc <- nullSimsa
nullSimsc$prop <- "40/0/60"
nullSimsd <- nullSimsa
nullSimsd$prop <- "60/0/40"
nullSimse <- nullSimsa
nullSimse$prop <- "20/0/80"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allBtRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("auto","x","both","unrel","stat")
dat <- dat[dat$auto!="auto.burden",]
dat$auto <- NULL; dat$x <- NULL; dat$unrel <- NULL
nullSims2 <- dat %>%
  filter(stat=="pval") %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") %>%
  mutate(both=as.numeric(both))
nullSims2$stat <- "KEATS-bt"
colnames(nullSims2)[1:2] <- c("pvalue","model")

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"
nullSims5 <- nullSims2
nullSims5$prop <- "40/0/60"
nullSims6 <- nullSims2
nullSims6$prop <- "60/0/40"
nullSims7 <- nullSims2
nullSims7$prop <- "20/0/80"

allRes <- rbind(allRes,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4,nullSimsc,nullSimsd,nullSimse,
                nullSims5,nullSims6,nullSims7)

# now read in null skat results
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allSkatRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model=="both") %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(pvalue=as.numeric(pvalue)) %>%
  mutate(chr="X") 
nullSims$qstat <- NULL; nullSims$rho <- NULL
nullSims$model <- "KEATS"

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"
nullSimsc <- nullSimsa
nullSimsc$prop <- "40/0/60"
nullSimsd <- nullSimsa
nullSimsd$prop <- "60/0/40"
nullSimse <- nullSimsa
nullSimse$prop <- "20/0/80"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allSkatRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model=="both") %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") %>%
  mutate(pvalue=as.numeric(pvalue))
nullSims2$qstat <- NULL; nullSims2$rho <- NULL
nullSims2$model <- "KEATS"

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"
nullSims5 <- nullSims2
nullSims5$prop <- "40/0/60"
nullSims6 <- nullSims2
nullSims6$prop <- "60/0/40"
nullSims7 <- nullSims2
nullSims7$prop <- "20/0/80"

allRes <- rbind(allRes,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4,nullSimsc,nullSimsd,nullSimse,
                nullSims5,nullSims6,nullSims7)

write.table(allRes,file="allRes_8ped_nullPow_ar1_s2as2x05.txt",row.names=FALSE,quote=FALSE)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop, 100K for the null sims

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=72,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60","20/0/80","40/0/60","60/0/40"))

finalR2 <- finalR[is.element(finalR$prop,c("20/0/80","40/0/60","60/0/40")),]
# plot these results now
pdf("power_8ped_keatsbtskat_trunc.pdf",width=14)
ggplot(finalR2,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

rm(list=ls())


#####
# 20. Power results, keatso, chrX and auto for different mixes of +/- for 8pedFem

# a 2x3 plot, auto vs x on the rows, cols corresponding to props of:
# 5/35/60, 15/25/60, 20/20/60

# need to put null and true sim results together,
# indicators for null/true, x/auto and proportion, as well as model

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_8pedFem_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

# make histograms of the rho values at the same configurations
powTr <- pow[pow$type==TRUE&pow$model=="both",]
powTr$prop <- ordered(powTr$prop,levels=c("5/35/60","15/25/60","20/20/60"))
pdf("rho_unrel_keatso_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("")
dev.off()

powTr <- pow[pow$type==TRUE&pow$model=="auto",]
powTr$prop <- ordered(powTr$prop,levels=c("5/35/60","15/25/60","20/20/60"))
pdf("rho_unrel_monster_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("")
dev.off()

rm(list=ls())


#####
# 21. Power results, keatso, chrX and auto for different mixes of +/- for unrelated pedigree

# a 2x3 plot, auto vs x on the rows, cols corresponding to props of:
# 5/35/60, 15/25/60, 20/20/60

# need to put null and true sim results together,
# indicators for null/true, x/auto and proportion, as well as model

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_unrel_chr9_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_unrel_chr9_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_unrel_chr9_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_unrel_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
#  s <- summarize(pow,n())
#  plo[,1] <- data.frame(s)[,5]
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,(i)] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
#plo$sampSize <- data.frame(summarize(pow,n()))[,5]
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
allres$alpha <- rep(alpha,each=48)
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_unrel_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()


# make cis for power plots
# need to calculate the se for each alpha level
# then plot the power est +/- 1.96*se
s <- summarize(pow,n())
s$model[s$model=="auto"] <- "MONSTER"
s$model[s$model=="both"] <- "KEATS-O"
s$model[s$model=="x"] <- "KEATS-O X only"
s$model[s$model=="unrel"] <- "SKATO"
s$type[s$type==TRUE] <- "power"
s$type[s$type==FALSE] <- "null"
# now, get se for each s value, each alpha level
# n is either 10K or 100K, so get 2 se values, one for each of those
se.10k <- 1.96*sqrt(alpha*(1-alpha)/10000)
se.100k <- 1.96*sqrt(alpha*(1-alpha)/100000)

finalR$ciL.power <- finalR$power-se.10k
finalR$ciU.power <- finalR$power+se.10k
finalR$ciL.null <- finalR$null-se.100k
finalR$ciU.null <- finalR$null+se.100k

## plot with cis now
pdf("power_unrel_trunc_cis.pdf",width=14)
ggplot(finalR,aes(x=null,color=model)) + facet_grid(chr~prop) + 
  geom_line(aes(y=ciL.power),size=1) + geom_line(aes(y=ciU.power),size=1)+
  theme_bw() + scale_y_continuous(limits=c(0.9,1)) +
  scale_x_continuous(limits=c(0,0.05)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

rm(list=ls())


#####
# 22. Make type I error tables for all methods considered

library(dplyr); library(readr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

# compare keats-o with monster and skato, also keats-o x only
# make table with auto vs x, within each of the 3 pedigree structures
# consider alpha=0.001 and 1e-4

# now read in null simulations
dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(ped="unrel")

dat <- read_delim("nullSims_unrel_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(ped="unrel")

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims3 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(ped="8ped")

dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims4 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped")

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(ped="8pedFem")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(ped="8pedFem")

allTy <- rbind(nullSims,nullSims2,nullSims3,nullSims4,nullSims5,nullSims6)

pow <- allTy %>%
  group_by(model,ped,chr)
data.frame(summarize(pow,n())) # 200K for each configuration

alpha <- c(0.01,0.001)
plo <- matrix(NA,nrow=24,ncol=length(alpha))

for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,4]
}
plo <- cbind(s[,1:3],plo)
n <- data.frame(summarize(pow,n()))[,4]
se.1 <- 1.96*sqrt(alpha[1]*(1-alpha[1])/n)
se.2 <- 1.96*sqrt(alpha[2]*(1-alpha[2])/n)

cbind(plo,inci=plo[,4]+se.1>alpha[1]&plo[,4]-se.1<alpha[1])


## look at the extreme sigma null values

dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/sim_x3a03_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped")

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/x3a03_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(ped="8ped")

allTy <- rbind(nullSims5,nullSims6)

pow <- allTy %>%
  group_by(model,ped,chr)
data.frame(summarize(pow,n())) # 200K for each configuration

alpha <- c(0.01,0.001)
plo <- matrix(NA,nrow=8,ncol=length(alpha))

for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,4]
}
plo <- cbind(s[,1:3],plo)
n <- data.frame(summarize(pow,n()))[,4]
se.1 <- 1.96*sqrt(alpha[1]*(1-alpha[1])/n)
se.2 <- 1.96*sqrt(alpha[2]*(1-alpha[2])/n)

cbind(plo,inci=plo[,4]+se.1>alpha[1]&plo[,4]-se.1<alpha[1])

## look at the s2x=0, s2a=1 results

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h01_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(ped="8ped")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h01_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(ped="8ped")

allTy <- rbind(nullSims5,nullSims6)

pow <- allTy %>%
  group_by(model,ped,chr)
data.frame(summarize(pow,n())) # 200K for each configuration

alpha <- c(0.01,0.001)
plo <- matrix(NA,nrow=8,ncol=length(alpha))

for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,4]
}
plo <- cbind(s[,1:3],plo)
n <- data.frame(summarize(pow,n()))[,4]
se.1 <- 1.96*sqrt(alpha[1]*(1-alpha[1])/n)
se.2 <- 1.96*sqrt(alpha[2]*(1-alpha[2])/n)

cbind(plo,inci=plo[,4]+se.1>alpha[1]&plo[,4]-se.1<alpha[1])
# ok, these look good! as expected

rm(list=ls())


#####
# 23. Make type I error tables for various sig2x, sig2a values

library(dplyr); library(readr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

# compare keats-o with monster and skato, also keats-o x only
# make table with auto vs x - only considering the 8pedFem pedigree structure now
# consider alpha=0.01 and 0.001

# first, look at s2a=1, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h01_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="01") 

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h01_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="01")

# now, s2a=0, s2x=1

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h10_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims3 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="10")

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h10_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims4 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="10")

# now, s2a=1, s2x=1

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h11_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="11")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h11_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="11")

# now, s2a=0, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h00_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims7 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="00")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h00_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims8 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="00")

# now, extreme sigma values

dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/sim_x3a03_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims9 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="303")

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/x3a03_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims10 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="303")

allTy <- rbind(nullSims,nullSims2,nullSims3,nullSims4,nullSims5,nullSims6,nullSims7,nullSims8,
               nullSims9,nullSims10)
allTy$model[allTy$model=="auto"] <- "MONSTER"
allTy$model[allTy$model=="both"] <- "KEATS-O"
allTy$model[allTy$model=="x"] <- "KEATS-O X only"
allTy$model[allTy$model=="unrel"] <- "SKATO"

pow <- allTy %>%
  group_by(chr,s2xa,model)
data.frame(summarize(pow,n())) # 100K for each configuration

alpha <- c(0.05,0.01,0.001)
plo <- matrix(NA,nrow=40,ncol=length(alpha))

for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,4]
}
plo <- cbind(s[,1:3],plo)
n <- 100000
se.0 <- 1.96*sqrt(alpha[1]*(1-alpha[1])/n)
se.1 <- 1.96*sqrt(alpha[2]*(1-alpha[2])/n)
se.2 <- 1.96*sqrt(alpha[3]*(1-alpha[3])/n)

cbind(plo,inci.0=plo[,4]+se.0>alpha[1]&plo[,4]-se.0<alpha[1],
      inci.1=plo[,5]+se.1>alpha[2]&plo[,5]-se.1<alpha[2],
      inci.2=plo[,6]+se.2>alpha[3]&plo[,6]-se.2<alpha[3])

colnames(plo)[4:6] <- paste("alpha",alpha,sep="_")

# alpha=0.05
totex <- plo %>% 
  mutate(alpha_0.01=NULL) %>%
  mutate(alpha_0.001=NULL) %>%
  spread(model,alpha_0.05)

# make this into latex table
library(xtable)
x <- xtable(totex,digits=3)
print(x,include.rownames=FALSE)

# alpha=0.01
totex <- plo %>% 
  mutate(alpha_0.05=NULL) %>%
  mutate(alpha_0.001=NULL) %>%
  spread(model,alpha_0.01)
x <- xtable(totex,digits=3)
print(x,include.rownames=FALSE)

# alpha=0.001
totex <- plo %>% 
  mutate(alpha_0.05=NULL) %>%
  mutate(alpha_0.01=NULL) %>%
  spread(model,alpha_0.001)
x <- xtable(totex,digits=4)
print(x,include.rownames=FALSE)

# make graphs of these results, too
tyIerrSm <- plo[,c(1:3,6)]
tyIerrSm <- tyIerrSm[!is.element(tyIerrSm$s2xa,"303"),]
tyIerrSm <- tyIerrSm[tyIerrSm$chr=="X",]
se <- 1.96*sqrt(0.001*(1-0.001)/100000)
tyIerrSm$lower <- tyIerrSm$alpha_0.001-se
tyIerrSm$upper <- tyIerrSm$alpha_0.001+se
tyIerrSm$sigma_x <- substr(tyIerrSm$s2xa,start=1,stop=1)  
tyIerrSm$sigma_a <- substr(tyIerrSm$s2xa,start=2,stop=nchar(tyIerrSm$s2xa))

plot_labeller <- function(variable,value){
  labels <- c('0'='0','1'='1')
  if (variable=='sigma_x') {
    v1 <- bquote(sigma[X]^2~"="~.(labels[value[1]]))
    v2 <- bquote(sigma[X]^2~"="~.(labels[value[2]]))
    return(c(v1,v2))
  } else {
    v1 <- bquote(sigma[A]^2~"="~.(labels[value[1]]))
    v2 <- bquote(sigma[A]^2~"="~.(labels[value[2]]))
    return(c(v1,v2))
  }
}
pdf("typeIErr_8pedFem_001.pdf")
ggplot(tyIerrSm,aes(x=model,y=alpha_0.001)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=0.001,color="gray"))+
  facet_grid(sigma_x~sigma_a,labeller=plot_labeller) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Model, Testing X Chr SNP, 20 Variants") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.001")))
dev.off()

tyIerrSm <- plo[,c(1:3,6)]
tyIerrSm <- tyIerrSm[!is.element(tyIerrSm$s2xa,"303"),]
tyIerrSm <- tyIerrSm[tyIerrSm$chr=="Autosomal",]
se <- 1.96*sqrt(0.001*(1-0.001)/100000)
tyIerrSm$lower <- tyIerrSm$alpha_0.001-se
tyIerrSm$upper <- tyIerrSm$alpha_0.001+se
tyIerrSm$sigma_x <- substr(tyIerrSm$s2xa,start=1,stop=1)  
tyIerrSm$sigma_a <- substr(tyIerrSm$s2xa,start=2,stop=nchar(tyIerrSm$s2xa))

pdf("typeIErr_8pedFem_001_autoSNP.pdf")
ggplot(tyIerrSm,aes(x=model,y=alpha_0.001)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=0.001,color="gray"))+
  facet_grid(sigma_x~sigma_a,labeller=plot_labeller) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Model, Testing Autosomal SNP, 20 Variants") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.001")))
dev.off()

rm(list=ls())


#####
# 24. Get SNPs in known RBC hits from HCHS/SOL data

library(dplyr); library(tidyr)
library(ggplot2); library(readr)

library(GenomicFeatures); library(QCannot)
library(biomaRt)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")
hits <- read.table("olga_rbc_assoc/RBC_trait_gwasCatalog_results_22sept15.tsv",sep="\t",header=T,as.is=T)
dim(hits) # 31 33
table(hits$Disease.Trait) # maybe subset to red blood cell count + red blood cell traits, only?
table(hits$Initial.Sample.Size) # european, but also af am and asian
summary(hits$p.Value) # all less than 8e-06
table(hits$Mapped_gene,hits$Reported.Gene.s.) 
table(hits$Reported.Gene.s.)
table(hits$Chr_id)
#1  2  4  5  6  7  9 12 16 20 23 
#3  1  3  1  9  3  2  2  3  1  3

# try subsetting
hitsSm <- hits[is.element(hits$Disease.Trait,c("Red blood cell count","Red blood cell traits")),]
dim(hitsSm) # 19 33
table(hitsSm$Initial.Sample.Size) # european, but also af am and asian
summary(hitsSm$p.Value) # all less than 4e-08
table(hitsSm$Mapped_gene,hitsSm$Reported.Gene.s.) 
table(hitsSm$Reported.Gene.s.)
table(hitsSm$Chr_id)
#1  2  4  5  6  7  9 12 16 23 
#2  1  2  1  5  1  2  2  1  2 

# this is a good place to start
# join the mapped_gene and reported.gene.s columns
hitsSm$Mapped_gene[hitsSm$Reported.Gene.s.=="CCND2"] <- "CCND2"
hitsSm$Mapped_gene[hitsSm$Reported.Gene.s.=="HBS1L, MYB"] <- "MYB"
hitsSm$Mapped_gene[hitsSm$Reported.Gene.s.=="ABO"] <- "ABO"
hitsSm$Mapped_gene[hitsSm$Reported.Gene.s.=="PDGFRA, HK1"] <- "PDGFRA"
hitsSm$Mapped_gene[hitsSm$Reported.Gene.s.=="KIT"] <- "KIT"
hitsSm$Mapped_gene[hitsSm$Reported.Gene.s.=="KITLG"] <- "KITLG"


# extract all polymorphic sites genotyped in HCHS/SOL that are within 100kb of these genes
geneList <- hitsSm[,c("Region", "Chr_id","Chr_pos","Mapped_gene","SNPs")]
# collapse the duplicate entries
dups <- duplicated(paste0(geneList$Chr_id,geneList$Chr_pos))
geneList <- geneList[!dups,]
geneList <- geneList[order(geneList$Chr_id,geneList$Chr_pos),]
# the 2 chr 4 hits are one SNP away in the same gene
# same with 2 of the chr 6 hits
# as are the 2 chr 9 hits
dups <- duplicated(geneList$Region)
geneList <- geneList[!dups,]
# great! looks good

# first, check snpannot for these rsIDs
snpAnnot <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData")

sum(is.element(geneList$SNPs,snpAnnot$rsID)) # 8
geneList$pos_37 <- geneList$Chr_pos
geneList$pos_37[geneList$SNPs=="rs7529925"] <- 199007208
# rs3811444 used to be rs17672233; is that one in SOL?
sum(is.element("rs17672233",snpAnnot$rsID)) # nope
geneList$pos_37[geneList$SNPs=="rs3811444"] <- 248039451
# rs7775698 used to be rs59556344; is that one in SOL?
sum(is.element("rs59556344",snpAnnot$rsID)) # nope
geneList$pos_37[geneList$SNPs=="rs7775698"] <- 135418635
# rs11611647 used to be rs57707142; is that in SOL?
sum(is.element("rs57707142",snpAnnot$rsID)) # nope
geneList$pos_37[geneList$SNPs=="rs11611647"] <- 4333919
# rs13335629 used to be rs56444421; is that in SOL?
sum(is.element("rs56444421",snpAnnot$rsID))
geneList$pos_37[geneList$SNPs=="rs13335629"] <- 310380
# rs1050828 used to be rs1894404, rs2230034, rs3191188; are any of these in SOL?
sum(is.element(c("rs1894404","rs2230034","rs3191188"),snpAnnot$rsID))
geneList$pos_37[geneList$SNPs=="rs1050828"] <- 154536002

# doesn't have first 2 on chr1, the last one on chr6, nor the first one on chr12
# check for SNPs nearby these 4 hits
sum(snpAnnot$chromosome==1&snpAnnot$position<(geneList$Chr_pos[1]+50)&
      snpAnnot$position>(geneList$Chr_pos[1]-50)) # 1; check it
pData(snpAnnot)[snpAnnot$chromosome==1&snpAnnot$position<(geneList$Chr_pos[1]+50)&
           snpAnnot$position>(geneList$Chr_pos[1]-50),] # looks good but has high MAF

geneList$sol.snps <- geneList$SNPs
geneList$sol.snps[1] <- snpAnnot$rsID[snpAnnot$chromosome==1&snpAnnot$position<(geneList$Chr_pos[1]+50)&
                                          snpAnnot$position>(geneList$Chr_pos[1]-50)]
# examine the other chr 1 hit
sum(snpAnnot$chromosome==1&snpAnnot$position<(geneList$Chr_pos[2]+100)&
      snpAnnot$position>(geneList$Chr_pos[2]-100)) # none; hmm
geneList$sol.snps[2] <- NA

# last chr6 hit
sum(snpAnnot$chromosome==6&snpAnnot$position<(geneList$Chr_pos[8]+100)&
      snpAnnot$position>(geneList$Chr_pos[8]-100)) # none; hmm
geneList$sol.snps[8] <- NA

# first chr 12 hit
sum(snpAnnot$chromosome==12&snpAnnot$position<(geneList$Chr_pos[11]+100)&
      snpAnnot$position>(geneList$Chr_pos[11]-100)) # none; hmm
geneList$sol.snps[11] <- NA

# the chr 1 rsID has gene N/A
# from dbSNP, a quick search said it's in the LINC01221 gene
geneList$Mapped_gene[geneList$Region=="1q32.1"] <- "LINC01221"

# remove the geneList entries for rsIDs that have Mapped_gene "-"
# they are intergenic and don't map to a gene
geneList<- geneList[!grepl("-",geneList$Mapped_gene),]
# the Mapped_gene with "," make 2 lines
geneList <- rbind(geneList,geneList[8,])
geneList$Mapped_gene[8] <- "TFR2"
geneList$Mapped_gene[14] <- "ACTL6B"

# now, need to add in the entrezgene ids so we can map to olga genotyped snps in these genes
# do by hand using geneCards online
geneList$entrezgene[geneList$Mapped_gene=="LINC01221"] <- 104266961
geneList$entrezgene[geneList$Mapped_gene=="TRIM58"] <- 25893 
geneList$entrezgene[geneList$Mapped_gene=="PRKCE"] <- 5581
geneList$entrezgene[geneList$Mapped_gene=="PDGFRA"] <- 5156
geneList$entrezgene[geneList$Mapped_gene=="TERT"] <- 7015
geneList$entrezgene[geneList$Mapped_gene=="CCND3"] <- 896
geneList$entrezgene[geneList$Mapped_gene=="MYB"] <- 4602
geneList$entrezgene[geneList$Mapped_gene=="TFR2"] <- 7036
geneList$entrezgene[geneList$Mapped_gene=="ABO"] <- 28
geneList$entrezgene[geneList$Mapped_gene=="CCND2"] <- 894
geneList$entrezgene[geneList$Mapped_gene=="KITLG"] <- 4254
geneList$entrezgene[geneList$Mapped_gene=="ITFG3"] <- 83986
geneList$entrezgene[geneList$Mapped_gene=="G6PD"] <- 2539
geneList$entrezgene[geneList$Mapped_gene=="ACTL6B"] <- 51412

# add a new SNP, reported in "GWAS of blood cell traits identifies novel associated loci and epistatic interactions in Caucasian and African-American children."
geneList <- rbind(geneList,c("16p13.3",16,124390,"NPRL3","rs7203560",NA,NA,8131))

# ok, see how many SOL SNPs are in these regions +/- 100 kb
txdb <- getTxDb("hg19", "knownGene")

# only keep monomorphic SNPs; also remove SNPs not in chrs 1-23
smSNP <- snpAnnot[snpAnnot$MAF.study>0&snpAnnot$MAF.study<0.05&snpAnnot$composite.filter,] # 2440128 36
dim(smSNP) # 897163

out <- defineExomeVars(smSNP, txdb)
dim(out) # 897163 8; same as snpAnnot
sum(is.na(out$geneID)) # 7190

## now, check out against geneList and see number of variants we have in each gene
geneList$numVars <- apply(geneList,1,function(x){sum(is.element(out$geneID,as.integer(x["entrezgene"])))})
# great, so these are all SNPs with MAF<0.05 that aren't monomorphic and pass QC

# get start/end position for genes
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
info <- getGene(geneList$Mapped_gene,type="hgnc_symbol",mart=mart)
info <- info[-c(4,10),] # dup entries

# use the genome build 37 mart
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
info37 <- getBM(attributes=c('hgnc_symbol','entrezgene','chromosome_name','start_position','end_position'),
                filters='entrezgene',values=geneList$entrezgene,mart=grch37)
info37 <- info37[!is.element(info37$chromosome_name,c("HG1497_PATCH","HG79_PATCH","LRG_148","LRG_343")),]

geneList <- merge(geneList,info37[,c("hgnc_symbol","start_position","end_position")],by.x="Mapped_gene",
                  by.y="hgnc_symbol",all.x=TRUE)
# got NA for LINC01221; look up on UCSC and add in
geneList$start_position[geneList$Mapped_gene=="LINC01221"] <- info$start_position[info$hgnc_symbol=="LINC01221"]
geneList$end_position[geneList$Mapped_gene=="LINC01221"] <- info$end_position[info$hgnc_symbol=="LINC01221"]

geneList$num.vars.100kb <- apply(geneList,1,function(x){sum(smSNP$chromosome==as.integer(x["Chr_id"])&
                                                              smSNP$position>as.integer(x["start_position"])-100&
                                                              smSNP$position<as.integer(x["end_position"])+100)})

# save the two lists -- one of vectors of SNPs in out for each of the genes 
# and one of the vectors of SNPs w/in 100Mb of each of the gene start & end positions
snp_map <- vector("list",nrow(geneList))
snp_100kb <- vector("list",nrow(geneList))
for(i in 1:nrow(geneList)){
  snp_map[[i]] <- out$snpID[is.element(out$geneID,geneList$entrezgene[i])]
  snp_100kb[[i]] <- smSNP$snpID[smSNP$chromosome==geneList$Chr_id[i]&
                                  smSNP$position>geneList$start_position[i]-100&
                                  smSNP$position<geneList$end_position[i]+100]
}

# ok, so check the overlap between these
# loop through by hand; store the ids
allRes <- vector("list",nrow(geneList))
for(i in 1:nrow(geneList)){
  allRes[[i]] <- unique(c(snp_map[[i]],snp_100kb[[i]]))
}
geneList$numVars.toTest <- sapply(allRes,length)
names(allRes) <- geneList$Mapped_gene

# ok, nice! use the larger of the two counts
# save a list of gene names, and variants that map to them
write.table(geneList,"olga_rbc_assoc/geneList_processed.txt",quote=FALSE,row.names=FALSE)
save(allRes,file="olga_rbc_assoc/snpIDs_100kb_RBChits.RData")

rm(list=ls())


#####
# 25. Parse type I error results, 50 variants

library(dplyr); library(readr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

# compare keats-o with monster and skato, also keats-o x only
# make table with auto vs x - only considering the 8pedFem pedigree structure now
# consider alpha=0.01 and 0.001

# first, look at s2a=1, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h01_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="01") 

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h01_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="01")

# now, s2a=0, s2x=1

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h10_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims3 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="10")

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h10_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims4 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="10")

# now, s2a=1, s2x=1

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h11_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="11")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h11_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="11")

# now, s2a=0, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h00_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims7 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="00")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h00_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims8 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="00")

allTy <- rbind(nullSims,nullSims2,nullSims3,nullSims4,nullSims5,nullSims6,nullSims7,nullSims8)
allTy$model[allTy$model=="auto"] <- "MONSTER"
allTy$model[allTy$model=="both"] <- "KEATS-O"
allTy$model[allTy$model=="x"] <- "KEATS-O X only"
allTy$model[allTy$model=="unrel"] <- "SKATO"

pow <- allTy %>%
  group_by(chr,s2xa,model)
data.frame(summarize(pow,n())) # 100K for each configuration

alpha <- c(0.05,0.01,0.001)
plo <- matrix(NA,nrow=32,ncol=length(alpha))

for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,4]
}
plo <- cbind(s[,1:3],plo)
n <- 100000
se.0 <- 1.96*sqrt(alpha[1]*(1-alpha[1])/n)
se.1 <- 1.96*sqrt(alpha[2]*(1-alpha[2])/n)
se.2 <- 1.96*sqrt(alpha[3]*(1-alpha[3])/n)

cbind(plo,inci.0=plo[,4]+se.0>alpha[1]&plo[,4]-se.0<alpha[1],
      inci.1=plo[,5]+se.1>alpha[2]&plo[,5]-se.1<alpha[2],
      inci.2=plo[,6]+se.2>alpha[3]&plo[,6]-se.2<alpha[3])

colnames(plo)[4:6] <- paste("alpha",alpha,sep="_")

# alpha=0.05
totex <- plo %>% 
  mutate(alpha_0.01=NULL) %>%
  mutate(alpha_0.001=NULL) %>%
  spread(model,alpha_0.05)

# make this into latex table
library(xtable)
x <- xtable(totex,digits=3)
print(x,include.rownames=FALSE)

# alpha=0.01
totex <- plo %>% 
  mutate(alpha_0.05=NULL) %>%
  mutate(alpha_0.001=NULL) %>%
  spread(model,alpha_0.01)
x <- xtable(totex,digits=3)
print(x,include.rownames=FALSE)

# alpha=0.001
totex <- plo %>% 
  mutate(alpha_0.05=NULL) %>%
  mutate(alpha_0.01=NULL) %>%
  spread(model,alpha_0.001)
x <- xtable(totex,digits=4)
print(x,include.rownames=FALSE)

# make graphs of these results, too
tyIerrSm <- plo[,c(1:3,6)]
tyIerrSm <- tyIerrSm[tyIerrSm$chr=="X",]
se <- 1.96*sqrt(0.001*(1-0.001)/100000)
tyIerrSm$lower <- tyIerrSm$alpha_0.001-se
tyIerrSm$upper <- tyIerrSm$alpha_0.001+se
tyIerrSm$sigma_x <- substr(tyIerrSm$s2xa,start=1,stop=1)  
tyIerrSm$sigma_a <- substr(tyIerrSm$s2xa,start=2,stop=nchar(tyIerrSm$s2xa))

plot_labeller <- function(variable,value){
  labels <- c('0'='0','1'='1')
  if (variable=='sigma_x') {
    v1 <- bquote(sigma[X]^2~"="~.(labels[value[1]]))
    v2 <- bquote(sigma[X]^2~"="~.(labels[value[2]]))
    return(c(v1,v2))
  } else {
    v1 <- bquote(sigma[A]^2~"="~.(labels[value[1]]))
    v2 <- bquote(sigma[A]^2~"="~.(labels[value[2]]))
    return(c(v1,v2))
  }
}
pdf("typeIErr_8pedFem_001_50vars.pdf")
ggplot(tyIerrSm,aes(x=model,y=alpha_0.001)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=0.001,color="gray"))+
  facet_grid(sigma_x~sigma_a,labeller=plot_labeller) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Model, Testing X Chr SNP, 50 Variants") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.001")))
dev.off()

tyIerrSm <- plo[,c(1:3,6)]
tyIerrSm <- tyIerrSm[tyIerrSm$chr=="Autosomal",]
se <- 1.96*sqrt(0.001*(1-0.001)/100000)
tyIerrSm$lower <- tyIerrSm$alpha_0.001-se
tyIerrSm$upper <- tyIerrSm$alpha_0.001+se
tyIerrSm$sigma_x <- substr(tyIerrSm$s2xa,start=1,stop=1)  
tyIerrSm$sigma_a <- substr(tyIerrSm$s2xa,start=2,stop=nchar(tyIerrSm$s2xa))

pdf("typeIErr_8pedFem_001_autoSNP_50vars.pdf")
ggplot(tyIerrSm,aes(x=model,y=alpha_0.001)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=0.001,color="gray"))+
  facet_grid(sigma_x~sigma_a,labeller=plot_labeller) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Model, Testing Autosomal SNP, 50 Variants") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.001")))
dev.off()

rm(list=ls())


#####
# 26. Run KEATSO on each RBC gene region in HCHS/SOL

# ran /projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal/olga_rbc_assoc/geneList_rbc_assoc.R
# in batch

library(GWASTools)
library(OLGApipeline); library(OLGAanalysis)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal/")
olgaRes <- get(load("olga_rbc_assoc/rbc_snpIDs_100kb_RBChits_keatsoResults.RData"))
geneAnnot <- read.table("olga_rbc_assoc/geneList_processed.txt",header=T,as.is=T)

# first, merge in rsIDs from the snpIDs we have
snpList <- get(load("olga_rbc_assoc/snpIDs_100kb_RBChits.RData"))
snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
for(i in 1:length(snpList)){
  snpList[[i]] <- data.frame(snpID=snpList[[i]])
  snpList[[i]] <- merge(snpList[[i]],pData(snpAnnot)[,c("snpID","rsID")],by="snpID")
}

# merge in snpIDs that we need for olgaData to get the single SNP assoc ids
sort(unique(geneAnnot$Chr_id)) #  1  2  4  5  6  7  9 12 16 23
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=1) 
# need to merge in snpIDs for chr 1 only
i <- which(is.element(geneAnnot$Chr_id,1))
for(j in i){
  snpList[[j]] <- merge(snpList[[j]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")
}

snpAnnot <- getSnpAnnotation(olgaData, chromosome=2) 
# need to merge in snpIDs for chr 1 only
i <- which(is.element(geneAnnot$Chr_id,2))
snpList[[i]] <- merge(snpList[[i]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")

snpAnnot <- getSnpAnnotation(olgaData, chromosome=4) 
i <- which(is.element(geneAnnot$Chr_id,4))
snpList[[i]] <- merge(snpList[[i]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")

snpAnnot <- getSnpAnnotation(olgaData, chromosome=5) 
i <- which(is.element(geneAnnot$Chr_id,5))
snpList[[i]] <- merge(snpList[[i]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")

snpAnnot <- getSnpAnnotation(olgaData, chromosome=6) 
i <- which(is.element(geneAnnot$Chr_id,6))
for(j in i){
  snpList[[j]] <- merge(snpList[[j]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")
}

snpAnnot <- getSnpAnnotation(olgaData, chromosome=7) 
i <- which(is.element(geneAnnot$Chr_id,7))
for(j in i){
  snpList[[j]] <- merge(snpList[[j]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")
}

snpAnnot <- getSnpAnnotation(olgaData, chromosome=9) 
i <- which(is.element(geneAnnot$Chr_id,9))
snpList[[i]] <- merge(snpList[[i]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")

snpAnnot <- getSnpAnnotation(olgaData, chromosome=12) 
i <- which(is.element(geneAnnot$Chr_id,12))
for(j in i){
  snpList[[j]] <- merge(snpList[[j]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")
}

snpAnnot <- getSnpAnnotation(olgaData, chromosome=16) 
i <- which(is.element(geneAnnot$Chr_id,16))
for(j in i){
  snpList[[j]] <- merge(snpList[[j]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")
}

snpAnnot <- getSnpAnnotation(olgaData, chromosome=23) 
i <- which(is.element(geneAnnot$Chr_id,23))
snpList[[i]] <- merge(snpList[[i]],pData(snpAnnot)[,c("snpID","rsID")],by="rsID")

# alright, now snpList has snpIDs for both olgaData and the usual snpAnnot

# in geneAnnot, store rho value and pvalue from keats-o
geneAnnot$rho <- sapply(olgaRes,function(x){x[[1]]$rho[which.min(x[[1]]$pvalue)]})
geneAnnot$keats.pvalue <- sapply(olgaRes,function(x){x[[2]]})

# now, compare the keatso pvalue with the single variant results for each of the sets of snps,
# what is the smallest pvalue out of all of them?
# also merge in maf.study
# argh, need to get the other snpAnnot, actually, then go by rsID to get the snpAnnot Ids from what we need
for(i in 1:length(snpList)){
  chr <- geneAnnot$Chr_id[i]
  res <- getobj(paste0("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/Assoc/assoc_316987_chr",
                       chr,".RData"))
  snpList[[i]] <- merge(snpList[[i]],res[,c("snpID","MAF","n","pval")],by.x="snpID.y",by.y="snpID")
  rm(res)
}
# great! now we have the single snp pvalues too

geneAnnot$single.min.pval <- sapply(snpList,function(x){min(x$pval,na.rm=T)})
# hmm. not as promising as i wanted...

write.table(geneAnnot,"olga_rbc_assoc/geneList_processed.txt",row.names=FALSE,quote=FALSE)
save(snpList,file="olga_rbc_assoc/snpIDs_100kb_RBChits.RData")

# make an xtable of the geneAnnot results
library(xtable)
geneAnnot <- geneAnnot[order(geneAnnot$Chr_id),]
geneAnnot$keats.pvalue <- format(geneAnnot$keats.pvalue,scientific=TRUE,digits=3)
geneAnnot$single.min.pval <- format(geneAnnot$single.min.pval,scientific=TRUE,digits=3)
x <- xtable(geneAnnot[,c("Mapped_gene","Chr_id","Region","numVars.toTest","keats.pvalue","rho",
                         "single.min.pval")])
print(x,include.rownames=FALSE)

rm(list=ls())


#####
# 27. Power results, keatso, chrX and auto for different mixes of +/- for 8pedFem

# a 2x3 plot, auto vs x on the rows, cols corresponding to props of:
# 20/0/80, 40/0/60, 60/0/40

# need to put null and true sim results together,
# indicators for null/true, x/auto and proportion, as well as model

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_8pedFem_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

# make histograms of the rho values at the same configurations
powTr <- pow[pow$type==TRUE&pow$model=="both",]
powTr$prop <- ordered(powTr$prop,levels=c("20/0/80","40/0/60","60/0/40"))
pdf("rho_8pedFem_keatso_btProps_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("")
dev.off()

powTr <- pow[pow$type==TRUE&pow$model=="auto",]
powTr$prop <- ordered(powTr$prop,levels=c("20/0/80","40/0/60","60/0/40"))
pdf("rho_8pedFem_monster_btProps_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("")
dev.off()

rm(list=ls())


#####
# 28. Hist of p-values corresponding to #23.

library(dplyr); library(readr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

# compare keats-o with monster and skato, also keats-o x only
# make table with auto vs x - only considering the 8pedFem pedigree structure now
# consider alpha=0.01 and 0.001

# first, look at s2a=1, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h01_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="01") 

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h01_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="01")

# now, s2a=0, s2x=1

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h10_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims3 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="10")

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h10_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims4 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="10")

# now, s2a=1, s2x=1

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h11_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="11")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h11_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="11")

# now, s2a=0, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h00_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims7 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="00")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h00_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims8 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="00")

allTy <- rbind(nullSims,nullSims2,nullSims3,nullSims4,nullSims5,nullSims6,nullSims7,nullSims8)
allTy$model[allTy$model=="auto"] <- "MONSTER"
allTy$model[allTy$model=="both"] <- "KEATS-O"
allTy$model[allTy$model=="x"] <- "KEATS-O X only"
allTy$model[allTy$model=="unrel"] <- "SKATO"

pow <- allTy %>%
  group_by(chr,s2xa,model)
data.frame(summarize(pow,n())) # 100K for each configuration

# want histograms, one file of autosomal SNPs
# facet_wrap by sig2xa

# one file of x chr SNPs
# facet_wrap by sig2xa

tyIerrSm <- pow %>%
  filter(chr=="X") %>%
  filter(s2xa=="11")

pdf("typeIErr_8pedFem_001_pvalueHist.pdf")
ggplot(tyIerrSm,aes(x=pvalue)) + geom_histogram(binwidth=0.05) + 
  theme_bw() + facet_wrap(~model,scales="free_y") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Frequency") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("P-Value, Testing X Chr SNP, 20 Variants") + 
  ggtitle(expression(paste("P-Values for Null Sims, ", sigma[X]^2, "=1, ",sigma[A]^2,"=1")))
dev.off()


tyIerrSm <- pow %>%
  filter(chr=="Autosomal") %>%
  filter(s2xa=="11")

pdf("typeIErr_8pedFem_001_autoSNP_pvalueHist.pdf")
ggplot(tyIerrSm,aes(x=pvalue)) + geom_histogram(binwidth=0.05) + 
  theme_bw() + facet_wrap(~model,scales="free_y") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Frequency") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("P-Value, Testing Autosomal SNP, 20 Variants") + 
  ggtitle(expression(paste("P-Values for Null Sims, ", sigma[X]^2, "=1, ",sigma[A]^2,"=1")))
dev.off()

rm(list=ls())


#####
# 29. Hist of p-values corresponding to #25.

library(dplyr); library(readr)
library(tidyr); library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/keats_x/keats_optimal")

# compare keats-o with monster and skato, also keats-o x only
# make table with auto vs x - only considering the 8pedFem pedigree structure now
# consider alpha=0.01 and 0.001

# first, look at s2a=1, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h01_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="01") 

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h01_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="01")

# now, s2a=0, s2x=1

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h10_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims3 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="10")

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h10_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims4 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="10")

# now, s2a=1, s2x=1

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h11_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="11")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h11_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="11")

# now, s2a=0, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h00_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims7 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="00")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h00_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims8 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="00")

allTy <- rbind(nullSims,nullSims2,nullSims3,nullSims4,nullSims5,nullSims6,nullSims7,nullSims8)
allTy$model[allTy$model=="auto"] <- "MONSTER"
allTy$model[allTy$model=="both"] <- "KEATS-O"
allTy$model[allTy$model=="x"] <- "KEATS-O X only"
allTy$model[allTy$model=="unrel"] <- "SKATO"

pow <- allTy %>%
  group_by(chr,s2xa,model)
data.frame(summarize(pow,n())) # 100K for each configuration

# want histograms, one file of autosomal SNPs
# facet_wrap by sig2xa

# one file of x chr SNPs
# facet_wrap by sig2xa

tyIerrSm <- pow %>%
  filter(chr=="X") %>%
  filter(s2xa=="11")

pdf("typeIErr_8pedFem_001_50vars_pvalueHist.pdf")
ggplot(tyIerrSm,aes(x=pvalue)) + geom_histogram(binwidth=0.05) + 
  theme_bw() + facet_wrap(~model,scales="free_y") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Frequency") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("P-Value, Testing X Chr SNP, 50 Variants") + 
  ggtitle(expression(paste("P-Values for Null Sims, ", sigma[X]^2, "=1, ",sigma[A]^2,"=1")))
dev.off()

tyIerrSm <- pow %>%
  filter(chr=="Autosomal") %>%
  filter(s2xa=="11")

pdf("typeIErr_8pedFem_001_autoSNP_50vars_pvalueHist.pdf")
ggplot(tyIerrSm,aes(x=pvalue)) + geom_histogram(binwidth=0.05) + 
  theme_bw() + facet_wrap(~model,scales="free_y") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Frequency") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("P-Value, Testing Autosomal SNP, 50 Variants") + 
  ggtitle(expression(paste("P-Values for Null Sims, ", sigma[X]^2, "=1, ",sigma[A]^2,"=1")))
dev.off()

rm(list=ls())


#####
# 30. 

library(dplyr); library(readr)
library(tidyr); library(ggplot2)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

# compare keats-o with monster and skato, also keats-o x only
# make table with auto vs x - only considering the 8pedFem pedigree structure now
# consider alpha=0.01 and 0.001

# first, look at s2a=1, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h01_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="01") 

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h01_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="01")

# now, s2a=0, s2x=1

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h10_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims3 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="10")

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h10_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims4 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="10")

# now, s2a=1, s2x=1

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h11_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims5 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="11")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h11_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims6 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="11")

# now, s2a=0, s2x=0

dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/sim_h00_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims7 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="X") %>%
  mutate(s2xa="00")

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/sim_h00_v50_allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims8 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(chr="Autosomal") %>%
  mutate(s2xa="00")

allTy <- rbind(nullSims,nullSims2,nullSims3,nullSims4,nullSims5,nullSims6,nullSims7,nullSims8)
allTy$model[allTy$model=="auto"] <- "MONSTER"
allTy$model[allTy$model=="both"] <- "KEATS-O"
allTy$model[allTy$model=="x"] <- "KEATS-O X only"
allTy$model[allTy$model=="unrel"] <- "SKATO"

pow <- allTy %>%
  group_by(chr,s2xa,model)
data.frame(summarize(pow,n())) # 100K for each configuration

alpha <- c(0.05,0.01,0.001)
plo <- matrix(NA,nrow=32,ncol=length(alpha))

for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,4]
}
plo <- cbind(s[,1:3],plo)
n <- 100000
se.0 <- 1.96*sqrt(alpha[1]*(1-alpha[1])/n)
se.1 <- 1.96*sqrt(alpha[2]*(1-alpha[2])/n)
se.2 <- 1.96*sqrt(alpha[3]*(1-alpha[3])/n)

cbind(plo,inci.0=plo[,4]+se.0>alpha[1]&plo[,4]-se.0<alpha[1],
      inci.1=plo[,5]+se.1>alpha[2]&plo[,5]-se.1<alpha[2],
      inci.2=plo[,6]+se.2>alpha[3]&plo[,6]-se.2<alpha[3])

colnames(plo)[4:6] <- paste("alpha",alpha,sep="_")

# alpha=0.05
totex <- plo %>% 
  mutate(alpha_0.01=NULL) %>%
  mutate(alpha_0.001=NULL) %>%
  spread(model,alpha_0.05)

# make this into latex table
library(xtable)
x <- xtable(totex,digits=3)
print(x,include.rownames=FALSE)

# alpha=0.01
totex <- plo %>% 
  mutate(alpha_0.05=NULL) %>%
  mutate(alpha_0.001=NULL) %>%
  spread(model,alpha_0.01)
x <- xtable(totex,digits=3)
print(x,include.rownames=FALSE)

# alpha=0.001
totex <- plo %>% 
  mutate(alpha_0.05=NULL) %>%
  mutate(alpha_0.01=NULL) %>%
  spread(model,alpha_0.001)
x <- xtable(totex,digits=4)
print(x,include.rownames=FALSE)

# make graphs of these results, too
tyIerrSm <- plo[,c(1:3,6)]
tyIerrSm <- tyIerrSm[tyIerrSm$chr=="X",]
se <- 1.96*sqrt(0.001*(1-0.001)/100000)
tyIerrSm$lower <- tyIerrSm$alpha_0.001-se
tyIerrSm$upper <- tyIerrSm$alpha_0.001+se
tyIerrSm$sigma_x <- substr(tyIerrSm$s2xa,start=1,stop=1)  
tyIerrSm$sigma_a <- substr(tyIerrSm$s2xa,start=2,stop=nchar(tyIerrSm$s2xa))

plot_labeller <- function(variable,value){
  labels <- c('0'='0','1'='1')
  if (variable=='sigma_x') {
    v1 <- bquote(sigma[X]^2~"="~.(labels[value[1]]))
    v2 <- bquote(sigma[X]^2~"="~.(labels[value[2]]))
    return(c(v1,v2))
  } else {
    v1 <- bquote(sigma[A]^2~"="~.(labels[value[1]]))
    v2 <- bquote(sigma[A]^2~"="~.(labels[value[2]]))
    return(c(v1,v2))
  }
}
pdf("typeIErr_8pedFem_001_50vars.pdf")
ggplot(tyIerrSm,aes(x=model,y=alpha_0.001)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=0.001,color="gray"))+
  facet_grid(sigma_x~sigma_a,labeller=plot_labeller) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Model, Testing X Chr SNP, 50 Variants") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.001")))
dev.off()

tyIerrSm <- plo[,c(1:3,6)]
tyIerrSm <- tyIerrSm[tyIerrSm$chr=="Autosomal",]
se <- 1.96*sqrt(0.001*(1-0.001)/100000)
tyIerrSm$lower <- tyIerrSm$alpha_0.001-se
tyIerrSm$upper <- tyIerrSm$alpha_0.001+se
tyIerrSm$sigma_x <- substr(tyIerrSm$s2xa,start=1,stop=1)  
tyIerrSm$sigma_a <- substr(tyIerrSm$s2xa,start=2,stop=nchar(tyIerrSm$s2xa))

pdf("typeIErr_8pedFem_001_autoSNP_50vars.pdf")
ggplot(tyIerrSm,aes(x=model,y=alpha_0.001)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=0.001,color="gray"))+
  facet_grid(sigma_x~sigma_a,labeller=plot_labeller) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  xlab("Model, Testing Autosomal SNP, 50 Variants") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.001")))
dev.off()

rm(list=ls())


#####
# 31. Group SNPs into genes, genome-wide
# use freeze1 annotation files!

# ran in batch, by chr:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_rbc_assoc
# qsub -q olga.q -t 1-23 -N geneList batch_make_geneList.sh
# which calls make_geneList.R

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools)
library(ggplot2)

# read in all gene lists to make plots of some summaries
geneAnnot <- NULL
geneList <- NULL
for(i in 1:23){
  dat <- read.table(paste0("olga_rbc_assoc/geneList_chr",i,".txt"),header=T,as.is=T)
  geneAnnot <- rbind(geneAnnot,dat)
  
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_chr",i,".RData"))
  geneList <- c(geneList,dat)
}
dim(geneAnnot) # 16500 2
length(geneList) # 16500

t <- sapply(geneList,length)

library(ggplot2)
pdf("olga_rbc_assoc/rbc_assoc_geneList_hist_16500.pdf")
ggplot(data.frame(t),aes(t)) + geom_histogram(binwidth=10) + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Frequency") + xlab("SNPs per Gene") + 
  ggtitle("Number of SNPs per 16,500 Genes")
dev.off()


# make a plot of genes by chromosome
pdf("olga_rbc_assoc/rbc_assoc_geneList_histByChr_16500.pdf")
ggplot(geneAnnot,aes(chromosome)) + geom_histogram(binwidth=1) + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(1,23,1)) +
  ylab("Frequency") + xlab("Genes per Chromosome") + 
  ggtitle("Number of Genes per Chromosome")
dev.off()

rm(list=ls())


#####
# 32. Run KEATSO on HCHS/SOL genome-wide

# use gds files per chromosome that have all genotyped SNPs that pass quality filter, with sporadic missingness imputed
# by shapeit
# /projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/observed/SOL_obs_from_imp_chr-10.gds 

# need to write the script to loop through by chromosome
# then through genes that are on a particular chromosome

# called:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_rbc_assoc/
# qsub -q olga.q N rbcKeatsX batch_geneList_rbc_genomeWide_assoc.sh

# want results for:
# plot pvalue by genomic position - manh plot
# look at rho values for each gene
# look at nvar for each gene, also sample size
# compare with rho=0, rho=1 vs keatso pvalues
# compare with GWAS p-values?

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)
library(MASS)
library(SNPRelate)
library(survey)
library(kinship)
library(pbivnorm)
library(SKAT)
library(dplyr); library(tidyr)
library(ggplot2); library(readr)
library(GenomicFeatures); library(QCannot)
library(biomaRt)

geneAnnot <- NULL
for(i in 1:23){
  tmp <- read.table(paste0("olga_rbc_assoc/geneList_chr",i,".txt"),header=T,as.is=T)
  geneAnnot <- rbind(geneAnnot,tmp)
}
geneAnnot$p.value <- NA
geneAnnot$nvar <- NA
geneAnnot$nsamp <- NA
geneAnnot$optimal.rho <- NA

chrs <- 1:23
for(i in chrs){
  fn <- paste0("olga_rbc_assoc/rbc_snpIDs_genomewide_keatsoResults_chr",i,".RData")
  dat <- getobj(fn)
  # 2nd value of each list item is the gene p-value
  # $nvar is the number of vars used
  # $nsamp is number of samples used
  # do another loop through the list of genes
  for(j in 1:length(dat)){
    ind <- which(geneAnnot$geneID==dat[[j]]$gene)
    geneAnnot$p.value[ind] <- dat[[j]][[2]]
    geneAnnot$nvar[ind] <- dat[[j]]$nvar
    geneAnnot$nsamp[ind] <- dat[[j]]$nsamp
    opt <- dat[[j]][[1]]
    geneAnnot$optimal.rho[ind] <- opt$rho[which.min(opt$pvalue)]
  }
}

dim(geneAnnot) # 16498 6
sum(is.na(geneAnnot$p.value)) # 2


# make a manh plot of the p-values
png("olga_rbc_assoc/genomeWide_keatso_manh.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500)
dev.off()

# plot the nvar by p.value
pdf("olga_rbc_assoc/genomeWide_pvalue_byNvar.pdf")
ggplot(geneAnnot,aes(x=-log10(p.value),y=nvar))+geom_point()+theme_bw()
dev.off() # no trend, many of the sig pvalues are on the low end of the nvar

pdf("olga_rbc_assoc/genomeWide_optimalRho_hist.pdf")
ggplot(geneAnnot,aes(x=optimal.rho))+geom_histogram(binwidth=0.1)+theme_bw()+
  xlab(expression(paste("Optimal ",rho))) + scale_x_continuous(breaks=seq(0,1,0.1)) 
dev.off() 

write.table(geneAnnot,file="olga_rbc_assoc/geneList_genomewide_results.txt",row.names=FALSE,quote=FALSE)

# get the significant gene names
0.05/16500 # 3.03e-06
geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]
#       geneID chromosome      p.value nvar nsamp optimal.rho
# 9064    3043         11 5.527472e-12    4 12488         0.2
# 12274  55692         16 1.045970e-16   39 12489         0.0 
# 12275  83986         16 2.041289e-13   13 12489         0.0
# 12279   9727         16 1.575372e-08   43 12489         0.0
# 16494   2539         23 3.582132e-16    5 12488         0.0
# 16497   2157         23 1.324010e-07   15 12488         0.0

library(xtable)
print(xtable(geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

# add mtxt for the 6 most sig p-values
png("olga_rbc_assoc/genomeWide_keatso_manh_annotated.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500)
N <- nrow(geneAnnot)
chromstart <- which(c(1, diff(geneAnnot$chromosome)) == 1)
chromend <- c(chromstart[-1], N)
geneAnnot$x <- (1:N) + geneAnnot$chromosome * (chromend[1]/6)
hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]
text(x=hits$x,y=-log10(hits$p.value)-0.6,labels=c("HBB","LUC7L","ITFG3","RAB11FIP3","G6PD","F8"))
dev.off()

# 3043: HBB protein-coding gene
#The alpha (HBA) and beta (HBB) loci determine the structure of the 2 types of polypeptide chains in adult hemoglobin, Hb A. 
#The normal adult hemoglobin tetramer consists of two alpha chains and two beta chains. Mutant beta globin causes sickle cell anemia.
#Absence of beta chain causes beta-zero-thalassemia. Reduced amounts of detectable beta globin causes beta-plus-thalassemia. 
#The order of the genes in the beta-globin cluster is 5'-epsilon -- gamma-G -- gamma-A -- delta -- beta--3'. 
#[provided by RefSeq, Jul 2008]

# 55692: LUC7L protein-coding gene
#The LUC7L gene may represent a mammalian heterochromatic gene, encoding a putative RNA-binding protein similar to the 
#yeast Luc7p subunit of the U1 snRNP splicing complex that is normally required for 5-prime splice site selection 
#(Tufarelli et al., 2001 [PubMed 11170747]).[supplied by OMIM, Mar 2008]

# 83986: ITFG3 protein-coding gene -- * included in the previously published associations

# 9727: RAB11FIP3 protein-coding gene
# Proteins of the large Rab GTPase family (see RAB1A; MIM 179508) have regulatory roles in the formation, targeting, 
#and fusion of intracellular transport vesicles. RAB11FIP3 is one of many proteins that interact with and regulate Rab GTPases 
#(Hales et al., 2001 [PubMed 11495908]).[supplied by OMIM, Mar 2008]

# 2539: G6PD -- * included in the previously published associations

# 2157: F8 protein-coding gene
# This gene encodes coagulation factor VIII, which participates in the intrinsic pathway of blood coagulation; 
#factor VIII is a cofactor for factor IXa which, in the presence of Ca+2 and phospholipids, converts factor X to the 
#activated form Xa. This gene produces two alternatively spliced transcripts. Transcript variant 1 encodes a large glycoprotein, 
#isoform a, which circulates in plasma and associates with von Willebrand factor in a noncovalent complex. 
#This protein undergoes multiple cleavage events. Transcript variant 2 encodes a putative small protein, isoform b, 
#which consists primarily of the phospholipid binding domain of factor VIIIc. 
#This binding domain is essential for coagulant activity. Defects in this gene results in hemophilia A, a common recessive 
#X-linked coagulation disorder. [provided by RefSeq, Jul 2008]
# * mentioned in COGENT paper

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
info <- getGene(geneList$Mapped_gene,type="hgnc_symbol",mart=mart)
hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]
info <- getBM(attributes=c('hgnc_symbol','entrezgene','chromosome_name','start_position','end_position'),
      filters='entrezgene',values=hits$geneID,mart=mart)

hits <- merge(hits,info,by.x="geneID",by.y="entrezgene")
# all the chr 16 hits are pretty close to eachother, the x chr hits are also neighboring

library(xtable)
print(xtable(hits[,c(7,8,9,4,3,6)],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

# read in single snp association p-values for these regions
# record the most sig in the region, also the bonf corrected sig p-value, based on #vars in gene region
# want chrs 11, 16, 23
# need SNPs for each gene region, too
geneList <- NULL
for(i in c(11,16,23)){
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_chr",i,".RData"))
  geneList <- c(geneList,dat)
}
geneList <- lapply(geneList,function(x){unique(x)})

assc11 <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr11.RData")
assc16 <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr16.RData")
assc23 <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr23.RData")

hits$sigP <- NA
hits$sigP.bonf <- NA
for(i in 1:nrow(hits)){
  snps <- geneList[as.character(hits$geneID[i])]  
  if(hits$chromosome[i]==11){
    pvals <- assc11$pval[is.element(assc11$snpID,unlist(snps))]
  }
  if(hits$chromosome[i]==16){
    pvals <- assc16$pval[is.element(assc16$snpID,unlist(snps))]
  }
  if(hits$chromosome[i]==23){
    pvals <- assc23$pval[is.element(assc23$snpID,unlist(snps))]
  }
  stopifnot(length(pvals)==length(unlist(snps)))
  hits$sigP[i] <- min(pvals,na.rm=T)
  hits$sigP.bonf[i] <- min(pvals,na.rm=T)*length(pvals)
  rm(snps)
}

#print(xtable(hits[,c(7,8,9,4,3,6,11,12)],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

## look up results for genes that were sig in the rare + common analysis
geneAnnot[is.element(geneAnnot$geneID,c(4354,100132963)),]

rm(list=ls())


#####
# 33. Power results, keatso, chrX and auto for different mixes of +/- for 8ped

# a 2x3 plot, auto vs x on the rows, cols corresponding to props of:
# 5/35/60, 15/25/60, 20/20/60

# need to put null and true sim results together,
# indicators for null/true, x/auto and proportion, as well as model

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="X",n=9892)
chr9_8ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim15256_c02_allOptRes.txt",prop="15/25/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim5356_c02_allOptRes.txt",prop="5/35/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim226_c02_allOptRes.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_8ped_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

# make histograms of the rho values at the same configurations
powTr <- pow[pow$type==TRUE&pow$model=="both",]
powTr$prop <- ordered(powTr$prop,levels=c("5/35/60","15/25/60","20/20/60"))
pdf("rho_8ped_keatso_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  theme_bw() +  ylab("")
dev.off()

powTr <- pow[pow$type==TRUE&pow$model=="auto",]
powTr$prop <- ordered(powTr$prop,levels=c("5/35/60","15/25/60","20/20/60"))
pdf("rho_8ped_monster_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("") + theme_bw()
dev.off()

rm(list=ls())

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_8ped_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

# make histograms of the rho values at the same configurations
powTr <- pow[pow$type==TRUE&pow$model=="both",]
powTr$prop <- ordered(powTr$prop,levels=c("20/0/80","40/0/60","60/0/40"))
pdf("rho_8ped_keatso_btProps_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("") + theme_bw()
dev.off()

powTr <- pow[pow$type==TRUE&pow$model=="auto",]
powTr$prop <- ordered(powTr$prop,levels=c("20/0/80","40/0/60","60/0/40"))
pdf("rho_8ped_monster_btProps_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("") + theme_bw()
dev.off()

rm(list=ls())


#####


library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/users/caitlin/keats_x/keats_optimal")


readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_unrel_chrX_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_unrel_chr9_ld05_ar1/sim406_c02_allOptRes.txt",prop="40/0/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_unrel_chr9_ld05_ar1/sim604_c02_allOptRes.txt",prop="60/0/40",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_unrel_chr9_ld05_ar1/sim208_c02_allOptRes.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_unrel_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_unrel_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

# make histograms of the rho values at the same configurations
powTr <- pow[pow$type==TRUE&pow$model=="both",]
powTr$prop <- ordered(powTr$prop,levels=c("20/0/80","40/0/60","60/0/40"))
pdf("rho_unrel_keatso_btProps_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("") + theme_bw()
dev.off()

powTr <- pow[pow$type==TRUE&pow$model=="auto",]
powTr$prop <- ordered(powTr$prop,levels=c("20/0/80","40/0/60","60/0/40"))
pdf("rho_unrel_monster_btProps_hist.pdf")
ggplot(powTr,aes(x=rho)) + geom_histogram(binwidth=0.1) + facet_grid(chr~prop,scales="free_y") + xlab(expression(paste(rho))) +
  ylab("") + theme_bw()
dev.off()

rm(list=ls())


#####
# 34. Hist of optimal rho values, null sims 8pedFem

library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  colnames(dat) <- c("qstat","pvalue","rho","model")
  dat <- dat[dat$model!="model",]
  
  datNew <- dat %>%
    mutate(qstat=as.numeric(qstat)) %>%
    mutate(pvalue=as.numeric(pvalue)) %>%
    mutate(rho=as.numeric(rho)) %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

# now read in null simulations
dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/allOptRes.txt",delim=" ",skip=1,col_names=FALSE)
colnames(dat) <- c("qstat","pvalue","rho","model")
nullSims2 <- dat %>%
  filter(!is.na(model)) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop; 200K for some null iterations

plo <- pow

plo <- data.frame(plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

# want histogram of plo$rho, stratified by chr
# only for model=KEATS-O
sm <- plo[plo$model=="KEATS-O",]
xres <- sm[sm$chr=="X",]
autores <- sm[sm$chr=="Autosomal",]
sm <- rbind(xres[1:100000,],autores[1:100000,])
pdf("rho_8pedFem_nullSims_keatso.pdf",width=14)
ggplot(sm,aes(x=rho)) +geom_histogram(binwidth=0.1) + theme_bw() + facet_grid(~chr) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text=element_text(size=16)) +
  scale_x_continuous(limits=c(0,1)) + xlab(expression(rho)) +ylab("Frequency")
dev.off()

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_8pedFem_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0.8,1)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

rm(list=ls())


#####
# 35. KEATSO on HCHS/SOL genome-wide QQ plots

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=T,as.is=T)

dim(geneAnnot) # 16498 6
sum(is.na(geneAnnot$p.value))

# make a qq plot of the p-values
png("olga_rbc_assoc/genomeWide_keatso_qq.png",width=600)
qqPlot(geneAnnot$p.value)
dev.off()

# can't calculate lambda_GC, ie. genomic control inflation factor, as we don't have the test stats

rm(list=ls())


#####
# 36. KEATSO on HCHS/SOL genome-wide single test p-value comparisons

# get the gene-based p-values, compare with min p-values per region, corrected for num vars in that gene region
# plot manh of the single test values, also a x-y plot with the corrected single test and the keatso values

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=T,as.is=T)

# read in all single test p-values and add to geneAnnot
geneList <- NULL
for(i in 1:23){
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_chr",i,".RData"))
  dat <- lapply(dat,function(x){unique(x)})
  assc <- getobj(paste0("../../olga_xchr_assoc/Assoc/assoc_316987_chr",i,".RData"))

  pvals <- lapply(dat,function(x){min(assc$pval[is.element(assc$snpID,x)],na.rm=T)})
  pres <- data.frame("singleTest.min.pval"=unlist(pvals))
  pres$gene <- names(pvals)
  
  geneAnnot <- merge(geneAnnot,pres,by.x="geneID",by.y="gene",all.x=TRUE)
}
# need to merge all the 'singleTest.min.pval.x and .y columns into one
for(i in 1:20){
geneAnnot$singleTest.min.pval[!is.na(geneAnnot[,7])] <- geneAnnot[!is.na(geneAnnot[,7]),7]
geneAnnot[,7] <- NULL
}

sum(is.na(geneAnnot$singleTest.min.pval)) # 0!
geneAnnot$singleTest.corr.pval <- geneAnnot$singleTest.min.pval*geneAnnot$nvar
geneAnnot$singleTest.corr.pval[geneAnnot$singleTest.corr.pval>1] <- 1

# reorder geneAnnot
tmp <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=T,as.is=T)
tmp$ord <- 1:nrow(tmp)
geneAnnot <- merge(geneAnnot,tmp[,c("geneID","ord")],by="geneID",all=TRUE)
geneAnnot <- geneAnnot[geneAnnot$ord,]

# save the new geneAnnot
write.table(geneAnnot,file="olga_rbc_assoc/geneList_genomewide_results.txt",row.names=FALSE,quote=FALSE)

# make manh plot of corr single test p-values
png("olga_rbc_assoc/genomeWide_singleTestCorrPval_manh.png",width=600)
manhattanPlot(geneAnnot$singleTest.corr.pval,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500)
dev.off()

# make both in same file
png("olga_rbc_assoc/genomeWide_keatso_singleTestCorr_manh.png",width=600,height=500)
par(mfrow=c(2,1))
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500,main="KEATS-O P-Value")
manhattanPlot(geneAnnot$singleTest.corr.pval,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500,main="Single Test Corrected P-Value")
dev.off()

# make x vs y plot of pvalues; WOW! 
png("olga_rbc_assoc/genomeWide_keatso_singleTestCorr_comparison.png",width=600,height=500)
ggplot(geneAnnot,aes(x=singleTest.corr.pval,y=p.value)) + geom_point() +
  theme_bw() + xlab("Single Test Corrected P-Value") + ylab("KEATS-O P-Value") +
  geom_abline(intercept=0,slope=1)
dev.off()

# color them by optimal rho
png("olga_rbc_assoc/genomeWide_keatso_singleTestCorr_comparison_colored.png",width=600,height=500)
ggplot(geneAnnot,aes(x=singleTest.corr.pval,y=p.value)) + geom_point(aes(color=optimal.rho)) +
  theme_bw() + xlab("Single Test Corrected P-Value") + ylab("KEATS-O P-Value") +
  geom_abline(intercept=0,slope=1)
dev.off()

rm(list=ls())


#####
# 37. Make a new HCHS/SOL gene list, NOT excluding SNPs with >5% MAF

# use freeze1 annotation files!

# ran in batch, by chr:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_rbc_assoc
# qsub -q olga.q -p -50 -t 1-23 -N geneList batch_make_geneList_inclSmallMAF.sh
# which calls make_geneList_inclSmallMAF.R

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools)
library(ggplot2)

# read in all gene lists to make plots of some summaries
geneAnnot <- NULL
geneList <- NULL
for(i in 1:23){
  dat <- read.table(paste0("olga_rbc_assoc/geneList_inclSmallMAF_chr",i,".txt"),header=T,as.is=T)
  geneAnnot <- rbind(geneAnnot,dat)
  
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_inclSmallMAF_chr",i,".RData"))
  geneList <- c(geneList,dat)
}
dim(geneAnnot) # 19629 2
length(geneList) # 19629

t <- sapply(geneList,length)
summary(t)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.00   11.00   24.00   58.76   55.00 4951.00 
which.max(t)# gene 64478; CSMD1 on 8p23.2
sum(t) # 1153339

library(ggplot2)
pdf("olga_rbc_assoc/rbc_assoc_geneList_inclSmallMAF_hist_19629.pdf")
ggplot(data.frame(t),aes(t)) + geom_histogram(binwidth=10) + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Frequency") + xlab("SNPs per Gene") + 
  ggtitle("Number of SNPs per 19,629 Genes")
dev.off()


# make a plot of genes by chromosome
pdf("olga_rbc_assoc/rbc_assoc_geneList_inclSmallMAF_histByChr_19629.pdf")
ggplot(geneAnnot,aes(chromosome)) + geom_histogram(binwidth=1) + theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(1,23,1)) +
  ylab("Frequency") + xlab("Genes per Chromosome") + 
  ggtitle("Number of Genes per Chromosome")
dev.off()

rm(list=ls())


#####
# 38. HCHS/SOL genome-wide using MONSTER, and SKAT on unrelateds

# MONSTER is calling KEATS-O with the auto KC matrix only
# called:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_rbc_assoc/
# qsub -q olga.q -p -50 -t 1-23 -N rbcKEATSx batch_geneList_rbc_genomeWide_autoKConly_assoc.sh

# SKAT is calling KEATS-O with no KC matrix, subset to unrelated samples
# called:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_rbc_assoc/
# qsub -q olga.q -p -50 -t 1-23 -N skatoRBC batch_geneList_rbc_genomeWide_skato_assoc.sh

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)
library(ggplot2)

# process skato results
geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=T,as.is=T)

geneAnnot$skato.nvar <- NA
geneAnnot$skato.nsamp <- NA
geneAnnot$skato.pval <- NA
geneAnnot$skato.optimalRho <- NA

chrs <- 1:23
for(i in chrs){
  fn <- paste0("olga_rbc_assoc/rbc_snpIDs_genomewide_skatoResults_chr",i,".RData")
  dat <- getobj(fn)
  # 2nd value of each list item is the gene p-value
  # $nvar is the number of vars used
  # $nsamp is number of samples used
  # do another loop through the list of genes
  for(j in 1:length(dat)){
    ind <- which(geneAnnot$geneID==dat[[j]]$gene)
    geneAnnot$p.value[ind] <- dat[[j]][[2]]
    geneAnnot$nvar[ind] <- dat[[j]]$nvar
    geneAnnot$nsamp[ind] <- dat[[j]]$nsamp
    opt <- dat[[j]][[1]]
    geneAnnot$optimal.rho[ind] <- opt$rho[which.min(opt$pvalue)]
  }
}

dim(geneAnnot) # 16498 6
sum(is.na(geneAnnot$p.value))




#####
# 39. KEATSO on HCHS/SOL genome-wide, NOT filtering out SNPs with >5% MAF

# use gds files per chromosome that have all genotyped SNPs that pass quality filter, with sporadic missingness imputed
# by shapeit
# /projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/observed/SOL_obs_from_imp_chr-10.gds 

# need to write the script to loop through by chromosome
# then through genes that are on a particular chromosome

# called:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_rbc_assoc/
# qsub -q olga.q -p -50 -t 1-23 -N rbcKEsm batch_geneList_rbc_genomeWide_assoc.sh

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)
library(MASS)
library(SNPRelate)
library(survey)
library(kinship)
library(pbivnorm)
library(SKAT)
library(dplyr); library(tidyr)
library(ggplot2); library(readr)
library(GenomicFeatures); library(QCannot)
library(biomaRt)

geneAnnot <- NULL
for(i in 1:23){
  tmp <- read.table(paste0("olga_rbc_assoc/geneList_inclSmallMAF_chr",i,".txt"),header=T,as.is=T)
  geneAnnot <- rbind(geneAnnot,tmp)
}
geneAnnot$p.value <- NA
geneAnnot$nvar <- NA
geneAnnot$nsamp <- NA
geneAnnot$optimal.rho <- NA
dim(geneAnnot) # 19629 6

chrs <- 1:23
for(i in chrs){
  fn <- paste0("olga_rbc_assoc/rbc_snpIDs_genomewide_inclSmallMAF_keatsoResults",i,".RData")
  dat <- getobj(fn)
  # 2nd value of each list item is the gene p-value
  # $nvar is the number of vars used
  # $nsamp is number of samples used
  # do another loop through the list of genes
  for(j in 1:length(dat)){
    ind <- which(geneAnnot$geneID==dat[[j]]$gene)
    geneAnnot$p.value[ind] <- dat[[j]][[2]]
    geneAnnot$nvar[ind] <- dat[[j]]$nvar
    geneAnnot$nsamp[ind] <- dat[[j]]$nsamp
    opt <- dat[[j]][[1]]
    geneAnnot$optimal.rho[ind] <- opt$rho[which.min(opt$pvalue)]
  }
}

dim(geneAnnot) # 19629 6
sum(is.na(geneAnnot$p.value)) # 0; good

# make a manh plot of the p-values
png("olga_rbc_assoc/genomeWide_inclSmallMAF_keatso_manh.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/19629)
dev.off()

# plot the nvar by p.value
pdf("olga_rbc_assoc/genomeWide_inclSmallMAF_pvalue_byNvar.pdf")
ggplot(geneAnnot,aes(x=-log10(p.value),y=nvar))+geom_point()+theme_bw()
dev.off() # no trend, many of the sig pvalues are on the low end of the nvar

pdf("olga_rbc_assoc/genomeWide_inclSmallMAF_optimalRho_hist.pdf")
ggplot(geneAnnot,aes(x=optimal.rho))+geom_histogram(binwidth=0.1)+theme_bw()+
  xlab(expression(paste("Optimal ",rho))) + scale_x_continuous(breaks=seq(0,1,0.1)) 
dev.off() 

write.table(geneAnnot,file="olga_rbc_assoc/geneList_inclSmallMAF_genomewide_results.txt",row.names=FALSE,quote=FALSE)

# get the significant gene names
0.05/19629 # 3.03e-06
geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>-log10(0.05/19629),]
#          geneID chromosome      p.value nvar nsamp optimal.rho
# 14229     55692         16 1.384604e-16  109 12489           0
# 14230     83986         16 1.640016e-13   26 12489           0
# 14240      9727         16 9.835044e-09   77 12489           0
# 19617      2539         23 3.264979e-16   12 12488           0
# 19620      4354         23 1.715940e-08   10 12488           0
# 19621 100132963         23 2.320606e-10    7 12488           0
# 19622      2157         23 8.232680e-08   43 12488           0

library(xtable)
print(xtable(geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

# add mtxt for the 6 most sig p-values
png("olga_rbc_assoc/genomeWide_inclSmallMAF_keatso_manh_annotated.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/19629)
N <- nrow(geneAnnot)
chromstart <- which(c(1, diff(geneAnnot$chromosome)) == 1)
chromend <- c(chromstart[-1], N)
geneAnnot$x <- (1:N) + geneAnnot$chromosome * (chromend[1]/6)
hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]
text(x=hits$x,y=-log10(hits$p.value)-0.4,labels=c("LUC7L","ITFG3","RAB11FIP3","G6PD","MPP1","SMIM9","F8"))
dev.off()

# read in single snp association p-values for these regions
# record the most sig in the region, also the bonf corrected sig p-value, based on #vars in gene region
# want chrs 16, 23
# need SNPs for each gene region, too
geneList <- NULL
for(i in c(16,23)){
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_inclSmallMAF_chr",i,".RData"))
  geneList <- c(geneList,dat)
}
geneList <- lapply(geneList,function(x){unique(x)})

assc16 <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr16.RData")
assc23 <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr23.RData")

hits$sigP <- NA
hits$sigP.bonf <- NA
for(i in 1:nrow(hits)){
  snps <- geneList[as.character(hits$geneID[i])]  
  if(hits$chromosome[i]==16){
    pvals <- assc16$pval[is.element(assc16$snpID,unlist(snps))]
  }
  if(hits$chromosome[i]==23){
    pvals <- assc23$pval[is.element(assc23$snpID,unlist(snps))]
  }
  stopifnot(length(pvals)==length(unlist(snps)))
  hits$sigP[i] <- min(pvals,na.rm=T)
  hits$sigP.bonf[i] <- min(pvals,na.rm=T)*length(pvals)
  rm(snps)
}

#print(xtable(hits[,c(7,8,9,4,3,6,11,12)],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

# get the pvalues for the rare analysis sig genes for table in dissertation
geneAnnot[is.element(geneAnnot$geneID,c(3043,83986)),]

rm(list=ls())


#####
# 40. KEATSO pvalues for genome-wide sig hits, compared to KEATS pvalues

# get the KEATS-O p-values for gene-based genome-wide results
# compare to KEATS-O p-values for when rho=0 and rho=1

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=T,as.is=T)

# merge in KEATS-bt and KEATS p-values
geneAnnot$keats.bt.pval <- NA
geneAnnot$keats.pval <- NA

chrs <- 1:23
for(i in chrs){
  fn <- paste0("olga_rbc_assoc/rbc_snpIDs_genomewide_keatsoResults_chr",i,".RData")
  dat <- getobj(fn)
  # 2nd value of each list item is the gene p-value
  # $nvar is the number of vars used
  # $nsamp is number of samples used
  # do another loop through the list of genes
  for(j in 1:length(dat)){
    ind <- which(geneAnnot$geneID==dat[[j]]$gene)
    geneAnnot$keats.pval[ind] <- dat[[j]][[1]][1,2]
    geneAnnot$keats.bt.pval[ind] <- dat[[j]][[1]][11,2]
    }
}

# so these are bt p-values using rho=1 in the keats-o analysis

# save the new geneAnnot
write.table(geneAnnot,file="olga_rbc_assoc/geneList_genomewide_results.txt",row.names=FALSE,quote=FALSE)

# look up the pvalues for the six sig genes
library(xtable)
geneAnnot[-log10(geneAnnot$p.value)>6,c("geneID","chromosome","chromosome","nvar","p.value","optimal.rho","keats.bt.pval","keats.pval")]
print(xtable(geneAnnot[-log10(geneAnnot$p.value)>6,c("geneID","chromosome","chromosome","nvar","p.value","optimal.rho","keats.bt.pval","keats.pval")],
             digits=c(0,0,0,0,0,0,0,0,1)),include.rownames=FALSE)

rm(list=ls())


#####
# 41. KEATSO on rare SNPs, platelet HCHS/SOL 

# analysis id 636058
# called:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_platelet_assoc/
# qsub -q olga.q -N platKeatsX -t 1-23 batch_geneList_platelet_genomeWide_assoc.sh

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)
library(MASS)
library(SNPRelate)
library(survey)c
library(kinship)
library(pbivnorm)
library(SKAT)
library(dplyr); library(tidyr)
library(ggplot2); library(readr)
library(GenomicFeatures); library(QCannot)
library(biomaRt)

geneAnnot <- NULL
for(i in 1:23){
  tmp <- read.table(paste0("olga_rbc_assoc/geneList_chr",i,".txt"),header=T,as.is=T)
  geneAnnot <- rbind(geneAnnot,tmp)
}
geneAnnot$p.value <- NA
geneAnnot$nvar <- NA
geneAnnot$nsamp <- NA
geneAnnot$optimal.rho <- NA

chrs <- 1:23
for(i in chrs){
  fn <- paste0("olga_platelet_assoc/platelet_snpIDs_genomewide_keatsoResults",i,".RData")
  dat <- getobj(fn)
  # 2nd value of each list item is the gene p-value
  # $nvar is the number of vars used
  # $nsamp is number of samples used
  # do another loop through the list of genes
  for(j in 1:length(dat)){
    ind <- which(geneAnnot$geneID==dat[[j]]$gene)
    geneAnnot$p.value[ind] <- dat[[j]][[2]]
    geneAnnot$nvar[ind] <- dat[[j]]$nvar
    geneAnnot$nsamp[ind] <- dat[[j]]$nsamp
    opt <- dat[[j]][[1]]
    geneAnnot$optimal.rho[ind] <- opt$rho[which.min(opt$pvalue)]
  }
}

dim(geneAnnot) # 16500 6
sum(is.na(geneAnnot$p.value)) # 114
table(is.na(geneAnnot$p.value),geneAnnot$chromosome) # chr 20

# make a manh plot of the p-values
png("olga_platelet_assoc/genomeWide_keatso_manh.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500)
dev.off()
# hmm only one hit, chr 19

# plot the nvar by p.value
pdf("olga_platelet_assoc/genomeWide_pvalue_byNvar.pdf")
ggplot(geneAnnot,aes(x=-log10(p.value),y=nvar))+geom_point()+theme_bw()
dev.off() # no trend, many of the sig pvalues are on the low end of the nvar

pdf("olga_platelet_assoc/genomeWide_optimalRho_hist.pdf")
ggplot(geneAnnot,aes(x=optimal.rho))+geom_histogram(binwidth=0.1)+theme_bw()+
  xlab(expression(paste("Optimal ",rho))) + scale_x_continuous(breaks=seq(0,1,0.1)) 
dev.off() 

write.table(geneAnnot,file="olga_platelet_assoc/geneList_genomewide_results.txt",row.names=FALSE,quote=FALSE)

# get the significant gene names
0.05/16500 # 3.03e-06
geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>5,]
#      geneID chromosome      p.value nvar nsamp optimal.rho
#14356   7171         19 1.893137e-06   10 12478           0

# is gene TPM4; already known...

library(xtable)
print(xtable(geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

# add mtxt for the sig p-value
png("olga_platelet_assoc/genomeWide_keatso_manh_annotated.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/16500)
N <- nrow(geneAnnot)
chromstart <- which(c(1, diff(geneAnnot$chromosome)) == 1)
chromend <- c(chromstart[-1], N)
geneAnnot$x <- (1:N) + geneAnnot$chromosome * (chromend[1]/6)
hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>5,]
text(x=hits$x,y=-log10(hits$p.value)+0.6,labels=c("TPM4"))
dev.off()

rm(list=ls())


#####
# 42. KEATSO on rare + common SNPs, platelet HCHS/SOL 

# analysis id 636058
# called:
# cd /projects/users/caitlin/keats_x/keats_optimal/olga_platelet_assoc/
# qsub -q olga.q -t 1-23 -N platKeSm batch_geneList_platelet_genomeWide_assoc.sh

setwd("/projects/users/caitlin/keats_x/keats_optimal")
library(GWASTools); library(QCpipeline)
library(OLGApipeline); library(OLGAanalysis)
library(MASS)
library(SNPRelate)
library(survey)
library(kinship)
library(pbivnorm)
library(SKAT)
library(dplyr); library(tidyr)
library(ggplot2); library(readr)
library(GenomicFeatures); library(QCannot)
library(biomaRt)

geneAnnot <- NULL
for(i in 1:23){
  tmp <- read.table(paste0("olga_rbc_assoc/geneList_inclSmallMAF_chr",i,".txt"),header=T,as.is=T)
  geneAnnot <- rbind(geneAnnot,tmp)
}
geneAnnot$p.value <- NA
geneAnnot$nvar <- NA
geneAnnot$nsamp <- NA
geneAnnot$optimal.rho <- NA

chrs <- 1:23
for(i in chrs){
  fn <- paste0("olga_platelet_assoc/platelet_snpIDs_genomewide_inclSmallMAF_keatsoResults",i,".RData")
  dat <- getobj(fn)
  # 2nd value of each list item is the gene p-value
  # $nvar is the number of vars used
  # $nsamp is number of samples used
  # do another loop through the list of genes
  for(j in 1:length(dat)){
    ind <- which(geneAnnot$geneID==dat[[j]]$gene)
    geneAnnot$p.value[ind] <- dat[[j]][[2]]
    geneAnnot$nvar[ind] <- dat[[j]]$nvar
    geneAnnot$nsamp[ind] <- dat[[j]]$nsamp
    opt <- dat[[j]][[1]]
    geneAnnot$optimal.rho[ind] <- opt$rho[which.min(opt$pvalue)]
  }
}

dim(geneAnnot) # 19629 6
sum(is.na(geneAnnot$p.value)) # 0

# make a manh plot of the p-values
png("olga_platelet_assoc/genomeWide_keatso_inclSmallMAF_manh.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/nrow(geneAnnot))
dev.off()
# hmm only one hit, chr 19

# plot the nvar by p.value
pdf("olga_platelet_assoc/genomeWide_inclSmallMAF_pvalue_byNvar.pdf")
ggplot(geneAnnot,aes(x=-log10(p.value),y=nvar))+geom_point()+theme_bw()
dev.off() # no trend, many of the sig pvalues are on the low end of the nvar

pdf("olga_platelet_assoc/genomeWide_inclSmallMAF_optimalRho_hist.pdf")
ggplot(geneAnnot,aes(x=optimal.rho))+geom_histogram(binwidth=0.1)+theme_bw()+
  xlab(expression(paste("Optimal ",rho))) + scale_x_continuous(breaks=seq(0,1,0.1)) 
dev.off() 

write.table(geneAnnot,file="olga_platelet_assoc/geneList_inclSmallMAF_genomewide_results.txt",row.names=FALSE,quote=FALSE)

# get the significant gene names
0.05/nrow(geneAnnot) # 2.547e-06
geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>5,]
#      geneID chromosome      p.value nvar nsamp optimal.rho
#14356   7171         19 1.893137e-06   10 12478           0

# is gene TPM4; already known...

library(xtable)
print(xtable(geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,],digits=c(0,0,0,0,0,0,1)),include.rownames=FALSE)

# add mtxt for the sig p-value
png("olga_platelet_assoc/genomeWide_keatso_inclSmallMAF_manh_annotated.png",width=600)
manhattanPlot(geneAnnot$p.value,geneAnnot$chromosome,trunc.lines=FALSE,ylim=c(0,16),
              signif=0.05/nrow(geneAnnot))
N <- nrow(geneAnnot)
chromstart <- which(c(1, diff(geneAnnot$chromosome)) == 1)
chromend <- c(chromstart[-1], N)
geneAnnot$x <- (1:N) + geneAnnot$chromosome * (chromend[1]/6)
hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>5,]
text(x=hits$x,y=-log10(hits$p.value)+0.6,labels=c("TPM4"))
dev.off()

rm(list=ls())


#####
# 43. Look at how the variants are mapped to genes

# how many of the variants are mapped to 1, 2, ... genes?
library(GenomicFeatures); library(QCannot)
library(biomaRt); library(ggplot2)
library(GWASTools)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

# first look at the rare gene regions
geneList <- NULL
for(i in 1:23){
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_chr",i,".RData"))
  geneList <- c(geneList,dat)
}
geneList <- lapply(geneList,function(x){unique(x)})
length(geneList) # 16500

# get a table of the count of each time the variant shows up in the geneList
t <- table(unlist(geneList))
table(t)
#     1      2      3      4      5      6      7      8      9     10     11 
#411829  18584   1030     98     20      8     12     24     41      8      3 
#12     13     14     15     16     17     18     19     20     21     22 
#17      9     24     57      2      1      5      7      6      2     25 

# so 411829 show up once
# 18584 show up twice
# 1030 show up 3x
# 98 show up 4x
# 25 show up 22x!

library(xtable)
print(xtable(table(t)),row.names=FALSE)

pdf("olga_rbc_assoc/hist_tableOfVars_perGene.pdf")
ggplot(data.frame(t)) + geom_histogram(aes(x=Freq)) + theme_bw()
dev.off()

toPl <- data.frame(t,stringsAsFactors=FALSE)
toPl$Freq <- as.integer(toPl$Freq)
toPl$Var1 <- as.integer(toPl$Var1)

pdf("olga_rbc_assoc/hist_tableOfVars_perGene_trunc.pdf")
ggplot(toPl,aes(x=Freq)) + ggtitle("Frequency of 431,812 SNPs Appearing in 16,500 Genes") +
  geom_histogram(binwidth=1) + theme_bw() + coord_cartesian(ylim=c(0,110)) + #ylim(c(0,100)) + 
  xlim(c(1,22)) +
  stat_bin(geom="text", aes(label=..count.., vjust=-1),binwidth=1) +
  annotate("text", x = 1.9, y = 90, label = "411829",size=3.8,color="white") + 
  annotate("text", x = 2.2, y = 85, label = "18584",size=3.8,color="white") +
  annotate("text", x = 3, y = 80, label = "1030",size=3.8,color="white")
dev.off()

# which are the snps that are mapped to many genes?
highFreq <- names(t)[t==22]
diff(as.integer(highFreq)) # they're all pretty close together

# try to figure out which genes they are in to find which chr SNP annot to look at
genIn <- sapply(geneList,function(x){sum(is.element(x,highFreq))})
genIn[genIn>0] # so here are the 22 genes
# a gene on chr 5, i think

#read in snp annot
snp <- getobj(paste0("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/annotation/SOL_freeze1_snpAnnot_chr-5.RData"))
snpSm <- snp[is.element(snp$snpID,highFreq),]
pData(snpSm)
# all on chr 5, between 140868085 - 140889260

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ann <- getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position",
               "ensembl_gene_id"), "entrezgene", names(genIn)[genIn>0], mart)
ann$hgnc_symbol
#[1] "PCDHGA12" "PCDHGC3"  "PCDHGC5"  "PCDHGC4"  "PCDHGB7"  "PCDHGB6" 
#[7] "PCDHGB5"  "PCDHGB3"  "PCDHGB2"  "PCDHGB1"  "PCDHGA11" "PCDHGA10"
#[13] "PCDHGA9"  "PCDHGA7"  "PCDHGA6"  "PCDHGA5"  "PCDHGA4"  "PCDHGA3" 
#[19] "PCDHGA2"  "PCDHGA1"  "PCDHGB4"  "PCDHGA8" 

##
# now look at rare + common genes
geneList <- NULL
for(i in 1:23){
  dat <- getobj(paste0("olga_rbc_assoc/snpIDs_genes_inclSmallMAF_chr",i,".RData"))
  geneList <- c(geneList,dat)
}
geneList <- lapply(geneList,function(x){unique(x)})
length(geneList) # 19629

# get a table of the count of each time the variant shows up in the geneList
t <- table(unlist(geneList))
table(t)
#      1       2       3       4       5       6       7       8       9      10 
#1000236   44884    2461     295      58      22      21      40      62      12 
#11      12      13      14      15      16      17      18      19      20 
# 6      23      17      46      73       2       2       9      18      11 
#21      22 
# 4      43 

# so 1000236 show up once
# 44884 show up twice
# 2461 show up 3x
# 295 show up 4x
# 58 show up 22x!

library(xtable)
print(xtable(table(t)),row.names=FALSE)

pdf("olga_rbc_assoc/hist_tableOfVars_perGene_inclSmallMAF.pdf")
ggplot(data.frame(t)) + geom_histogram(aes(x=Freq)) + theme_bw()
dev.off()

toPl <- data.frame(t,stringsAsFactors=FALSE)
toPl$Freq <- as.integer(toPl$Freq)
toPl$Var1 <- as.integer(toPl$Var1)

pdf("olga_rbc_assoc/hist_tableOfVars_perGene_inclSmallMAF_trunc.pdf")
ggplot(toPl,aes(x=Freq)) + ggtitle("Frequency of 1,048,345 SNPs Appearing in 19,629 Genes") +
  geom_histogram(binwidth=1) + theme_bw() + coord_cartesian(ylim=c(0,340)) + xlim(c(1,22)) +
  stat_bin(geom="text", aes(label=..count.., vjust=-1),binwidth=1) +
  annotate("text", x = 2.0, y = 330, label = "1000236",size=3.8,color="white") + 
  annotate("text", x = 2.2, y = 320, label = "44884",size=3.8,color="white") +
  annotate("text", x = 3, y = 310, label = "2461",size=3.8,color="white")
dev.off()

# which are the snps that are mapped to many genes?
highFreq <- names(t)[t==22]
diff(as.integer(highFreq)) # they're all pretty close together

# try to figure out which genes they are in to find which chr SNP annot to look at
genIn <- sapply(geneList,function(x){sum(is.element(x,highFreq))})
genIn[genIn>0] # so here are the 22 genes
# a gene on chr 5, i think

#read in snp annot
snp <- getobj(paste0("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/annotation/SOL_freeze1_snpAnnot_chr-5.RData"))
snpSm <- snp[is.element(snp$snpID,highFreq),]
# all on chr 5, between 140867347 - 140890770

ann <- getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position",
               "ensembl_gene_id"), "entrezgene", names(genIn)[genIn>0], mart)
ann$hgnc_symbol
#[1] "PCDHGA12" "PCDHGC3"  "PCDHGC5"  "PCDHGC4"  "PCDHGB7"  "PCDHGB6" 
#[7] "PCDHGB5"  "PCDHGB3"  "PCDHGB2"  "PCDHGB1"  "PCDHGA11" "PCDHGA10"
#[13] "PCDHGA9"  "PCDHGA7"  "PCDHGA6"  "PCDHGA5"  "PCDHGA4"  "PCDHGA3" 
#[19] "PCDHGA2"  "PCDHGA1"  "PCDHGB4"  "PCDHGA8" 
#protocadherin gamma subfamily A, 12; B, 3; C, 4; C, 5; ...

rm(list=ls())


#####
# 44. Follow up on hits from KEATSO rare analysis

# for each of the sig genes, look at the single var test results 
# what are the MAF of each of the variants?
# what is the p-value? 
# what do the cluster plots look like?
# how much of the variants are overlapping in the chr 16 and x chr genes?

library(GWASTools)
library(OLGApipeline); library(OLGAanalysis)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=TRUE,as.is=TRUE)
dim(geneAnnot) # 16502 13
head(geneAnnot)

hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]

# get the mapping of snps to genes
map16 <- get(load("olga_rbc_assoc/snpIDs_genes_chr16.RData"))
map16 <- sapply(map16,function(x){unique(x)})

hitsMap <- map16[is.element(names(map16),hits$geneID)]
snpRes <- data.frame(unlist(hitsMap))
len <- sapply(hitsMap,function(x){length(x)})
snpRes$gene <- rep(names(len),len)
snpRes <- merge(snpRes,hits,by.x="gene",by.y="geneID",all.x=TRUE)
colnames(snpRes)[2] <- "snpID"

## ok, do this for the chr 23 hits too
# get the mapping of snps to genes
map23 <- get(load("olga_rbc_assoc/snpIDs_genes_chr23.RData"))
map23 <- sapply(map23,function(x){unique(x)})

hitsMap <- map23[is.element(names(map23),hits$geneID)]
snpRes2 <- data.frame(unlist(hitsMap))
len <- sapply(hitsMap,function(x){length(x)})
snpRes2$gene <- rep(names(len),len)
snpRes2 <- merge(snpRes2,hits,by.x="gene",by.y="geneID",all.x=TRUE)
colnames(snpRes2)[2] <- "snpID"

snpRes <- rbind(snpRes,snpRes2)

sum(is.na(snpRes$singleTest.min.pval)) # 0
snpRes$singleTest.min.pval.x <- snpRes$singleTest.min.pval.y <- NULL

# make sure we aren't missing any single vars that were sig in the single SNP test in these regions
assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr16.RData")
assc <- assc[!is.na(assc$pval),]
assc <- assc[assc$pval<5e-8|is.element(assc$snpID,snpRes$snpID),]
dim(assc) # 296 17

snpRes <- merge(snpRes,assc[,c("snpID","n","Stat","pval")],by="snpID",all=TRUE)

assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr23.RData")
assc <- assc[!is.na(assc$pval),]
assc <- assc[assc$pval<5e-8|is.element(assc$snpID,snpRes$snpID),]
dim(assc) # 158 17

snpRes <- merge(snpRes,assc[,c("snpID","n","Stat","pval")],by="snpID",all=TRUE)

snpRes$n <- snpRes$n.x
snpRes$n[is.na(snpRes$n)] <- snpRes$n.y[is.na(snpRes$n)]
snpRes$n.x <- snpRes$n.y <- NULL

snpRes$Stat <- snpRes$Stat.x
snpRes$Stat[is.na(snpRes$Stat)] <- snpRes$Stat.y[is.na(snpRes$Stat)]
snpRes$Stat.x <- snpRes$Stat.y <- NULL

snpRes$single.test.pval <- snpRes$pval.x
snpRes$single.test.pval[is.na(snpRes$single.test.pval)] <- snpRes$pval.y[is.na(snpRes$single.test.pval)]
snpRes$pval.x <- snpRes$pval.y <- NULL

dim(snpRes) # 468 15

# get the snp annot
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=16) 

snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","snpID","position","chromosome")],by="snpID",all.x=TRUE)
dim(snpRes) # 468 18

snpRes$chromosome <- snpRes$chromosome.x
snpRes$chromosome[is.na(snpRes$chromosome)] <- snpRes$chromosome.y[is.na(snpRes$chromosome)]
snpRes$chromosome.x <- snpRes$chromosome.y <- NULL

snpAnnot <- getSnpAnnotation(olgaData, chromosome=23) 

snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","snpID","position","chromosome")],by="snpID",all.x=TRUE)
dim(snpRes) # 468 20

snpRes$position <- snpRes$position.x
snpRes$position[is.na(snpRes$position)] <- snpRes$position.y[is.na(snpRes$position)]
snpRes$position.x <- snpRes$position.y <- NULL

snpRes$rsID <- snpRes$rsID.x
snpRes$rsID[is.na(snpRes$rsID)] <- snpRes$rsID.y[is.na(snpRes$rsID)]
snpRes$rsID.x <- snpRes$rsID.y <- NULL

snpRes$chromosome <- snpRes$chromosome.x
snpRes$chromosome[is.na(snpRes$chromosome)] <- snpRes$chromosome.y[is.na(snpRes$chromosome)]
snpRes$chromosome.x <- snpRes$chromosome.y <- NULL


# now merge in the old SNP annot so we can get MAF
snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","exclude","MAF.study","composite.filter","quality.filter")],by="rsID",all.x=TRUE)
dim(snpRes) # 468 21

snpRes$geneName <- NA
snpRes$geneName[snpRes$gene==2157] <- "F8"
snpRes$geneName[snpRes$gene==55692] <- "LUC7L"
snpRes$geneName[snpRes$gene==83986] <- "ITFG3"
snpRes$geneName[snpRes$gene==9727] <- "RAB11FIP3"
snpRes$geneName[snpRes$gene==2539] <- "G6PD"

snpRes$geneName[is.na(snpRes$geneName)] <- "notMapped"

t <- table(snpRes$snpID,snpRes$geneName)
t[rowSums(t[,c(1:2)])>1,] # none
t[rowSums(t[,c(3:5)])>1,] # some SNPs overlap between itfg3 and luc7l
#          F8 G6PD ITFG3 LUC7L RAB11FIP3
# 22666746  0    0     1     1         0
# 22666782  0    0     1     1         0
# 22666842  0    0     1     1         0
# 22666903  0    0     1     1         0
# 22666907  0    0     1     1         0
# 22666972  0    0     1     1         0
# 22667033  0    0     1     1         0
# 22667039  0    0     1     1         0
# 22667083  0    0     1     1         0
# 22667103  0    0     1     1         0
# 22667108  0    0     1     1         0
# 22667132  0    0     1     1         0
# 22667146  0    0     1     1         0
# 

# make table of chr 16 significant from single snp vs LUC7L vs ITFG3 vs RAB11FIP3
sm <- snpRes[snpRes$chromosome==16,]
table(sm$geneName,exclude=NULL)

assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr16.RData")
assc <- assc[!is.na(assc$pval)&assc$pval<5e-8,]
dim(assc) # 225 17

library(VennDiagram)
todr <- list("SingleSNP"=assc$snpID,"LUC7L"=sm$snpID[sm$geneName=="LUC7L"],
             "ITFG3"=sm$snpID[sm$geneName=="ITFG3"],"RAB11FIP3"=sm$snpID[sm$geneName=="RAB11FIP3"])
venn.diagram(todr,filename="olga_rbc_assoc/venn_sigHits_chr16.tiff")

# do the same for the chr 23 hits
sm <- snpRes[snpRes$chromosome==23,]
table(sm$geneName,exclude=NULL)

assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr23.RData")
assc <- assc[!is.na(assc$pval)&assc$pval<5e-8,]
dim(assc) # 141 17

todr <- list("F8"=sm$snpID[sm$geneName=="F8"],"SingleSNP"=assc$snpID,
             "G6PD"=sm$snpID[sm$geneName=="G6PD"])
venn.diagram(todr,filename="olga_rbc_assoc/venn_sigHits_chr23.tiff")

# now, look at pvalues for single snp hits vs gene region test
# make latex tables for these genes
library(xtable)
snpRes$single.test.pval <- format(snpRes$single.test.pval,digits=3)
snpRes$MAF.study <- format(snpRes$MAF.study,digits=3)

print(xtable(snpRes[snpRes$geneName=="G6PD",c("rsID","single.test.pval","MAF.study")],digits=c(1,0,2,2)),
      include.rownames=FALSE)
print(xtable(snpRes[snpRes$geneName=="F8",c("rsID","single.test.pval","MAF.study")],digits=c(1,0,6,6)),
      include.rownames=FALSE)
print(xtable(snpRes[snpRes$geneName=="LUC7L",c("rsID","single.test.pval","MAF.study")],digits=c(1,0,6,6)),
      include.rownames=FALSE)
print(xtable(snpRes[snpRes$geneName=="ITFG3",c("rsID","single.test.pval","MAF.study")],digits=c(1,0,6,6)),
      include.rownames=FALSE)
print(xtable(snpRes[snpRes$geneName=="RAB11FIP3",c("rsID","single.test.pval","MAF.study")],digits=c(1,0,6,6)),
      include.rownames=FALSE)

##
# do this same thing for the HBB chr 11 hit
geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=TRUE,as.is=TRUE)
dim(geneAnnot) # 16502 13
head(geneAnnot)

hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]

# get the mapping of snps to genes
map11 <- get(load("olga_rbc_assoc/snpIDs_genes_chr11.RData"))
map11 <- sapply(map11,function(x){unique(x)})

hitsMap <- map11[is.element(names(map11),hits$geneID)]
snpRes <- data.frame(unlist(hitsMap))
len <- sapply(hitsMap,function(x){length(x)})
snpRes$gene <- rep(names(len),len)
snpRes <- merge(snpRes,hits,by.x="gene",by.y="geneID",all.x=TRUE)
colnames(snpRes)[2] <- "snpID"

sum(is.na(snpRes$singleTest.min.pval)) # 0
snpRes$singleTest.min.pval.x <- snpRes$singleTest.min.pval.y <- NULL

# make sure we aren't missing any single vars that were sig in the single SNP test in these regions
assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr11.RData")
assc <- assc[!is.na(assc$pval),]
assc <- assc[assc$pval<5e-8|is.element(assc$snpID,snpRes$snpID),]
dim(assc) # 13 17

snpRes <- merge(snpRes,assc[,c("snpID","n","Stat","pval")],by="snpID",all=TRUE)

# get the snp annot
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=11) 

snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","snpID","position","chromosome")],by="snpID",all.x=TRUE)
dim(snpRes) # 468 18

snpRes$chromosome <- snpRes$chromosome.x
snpRes$chromosome[is.na(snpRes$chromosome)] <- snpRes$chromosome.y[is.na(snpRes$chromosome)]
snpRes$chromosome.x <- snpRes$chromosome.y <- NULL

# now merge in the old SNP annot so we can get MAF
snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","exclude","MAF.study","composite.filter","quality.filter")],by="rsID",all.x=TRUE)
dim(snpRes) # 13 21

snpRes$geneName <- NA
snpRes$geneName[snpRes$gene==3043] <- "HBB"
snpRes$geneName[is.na(snpRes$geneName)] <- "notMapped"

t <- table(snpRes$snpID,snpRes$geneName)

# now, look at pvalues for single snp hits vs gene region test
# make latex tables for these genes
snpRes$pval <- format(snpRes$pval,digits=3)
snpRes$MAF.study <- format(snpRes$MAF.study,digits=3)
print(xtable(snpRes[snpRes$geneName=="HBB",c("rsID","pval","MAF.study")],digits=c(1,0,2,2)),
      include.rownames=FALSE)

rm(list=ls())


#####
# 45. Make table of HBB variants to send to RBC email list

library(GWASTools)
library(OLGApipeline); library(OLGAanalysis)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=TRUE,as.is=TRUE)
dim(geneAnnot) # 16502 13
head(geneAnnot)

hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]

# want the variants mapped to gene HBB 
# get the mapping of snps to genes
map11 <- get(load("olga_rbc_assoc/snpIDs_genes_chr11.RData"))
map11 <- sapply(map11,function(x){unique(x)})

hitsMap <- map11[is.element(names(map11),hits$geneID)]
snpRes <- data.frame(unlist(hitsMap))
len <- sapply(hitsMap,function(x){length(x)})
snpRes$gene <- rep(names(len),len)
snpRes <- merge(snpRes,hits,by.x="gene",by.y="geneID",all.x=TRUE)
colnames(snpRes)[2] <- "snpID"

sum(is.na(snpRes$singleTest.min.pval)) # 0
snpRes$singleTest.min.pval.x <- snpRes$singleTest.min.pval.y <- NULL

# make sure we aren't missing any single vars that were sig in the single SNP test in these regions
assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr11.RData")
assc <- assc[!is.na(assc$pval),]
assc <- assc[assc$pval<5e-8|is.element(assc$snpID,snpRes$snpID),]
dim(assc) # 13 17

snpRes <- merge(snpRes,assc[,c("snpID","n","Stat","pval")],by="snpID",all=TRUE)

# get the snp annot
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=11) 

snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","snpID","position","chromosome")],by="snpID",all.x=TRUE)
dim(snpRes) # 468 18

snpRes$chromosome <- snpRes$chromosome.x
snpRes$chromosome[is.na(snpRes$chromosome)] <- snpRes$chromosome.y[is.na(snpRes$chromosome)]
snpRes$chromosome.x <- snpRes$chromosome.y <- NULL

# now merge in the old SNP annot so we can get MAF
snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","exclude","MAF.study","composite.filter","quality.filter")],by="rsID",all.x=TRUE)
dim(snpRes) # 13 21

snpRes$geneName <- NA
snpRes$geneName[snpRes$gene==3043] <- "HBB"
snpRes$geneName[is.na(snpRes$geneName)] <- "notMapped"

snpRes[snpRes$geneName=="HBB",]

snpRes$rareAnalys[snpRes$geneName=="HBB"] <- TRUE

## merge in variant mapping for rare + common analysis
map11 <- get(load("olga_rbc_assoc/snpIDs_genes_inclSmallMAF_chr11.RData"))
map11 <- sapply(map11,function(x){unique(x)})

hbb <- snpRes[snpRes$geneName=="HBB",]
comVarhbb <- map11[["3043"]]
snpAnnot <- getSnpAnnotation(olgaData, chromosome=11) 
snpAnnothbb <- snpAnnot[is.element(snpAnnot$snpID,comVarhbb),]

snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
snpAnnothbb <- merge(pData(snpAnnothbb)[,c("snpID","rsID","position","chromosome")],pData(snpAnnot)[,c("rsID","exclude","MAF.study","composite.filter","quality.filter")],by="rsID",all.x=TRUE)

snpAnnothbb$commonAnalys <- TRUE

allRes <- merge(hbb[,c("rsID","MAF.study","rareAnalys","pval")],snpAnnothbb[,c("snpID","rsID","MAF.study","commonAnalys")],by="rsID",all=TRUE)

# need to merge in pval for rest of SNPs
snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
allRes <- merge(allRes,pData(snpAnnot)[,c("snpID","rsID")],by="rsID",all.x=TRUE)

assc <- getobj("../../olga_xchr_assoc/Assoc/assoc_316987_chr11.RData")
finalRes <- merge(allRes,assc[,c("snpID","pval")],by.x="snpID.y",by.y="snpID",all.x=TRUE)

finalRes$MAF.study.x <- NULL
finalRes$rareAnalys[is.na(finalRes$rareAnalys)] <- FALSE

toPr <- finalRes[,c("rsID","MAF.study.y","pval.y","rareAnalys")]
colnames(toPr) <- c("rsID","MAF","singleTest.pval","rareAnalys")
toPr

rm(list=ls())


#####
# 46. More power sims to add to dissertation

# cut sample size from 2000 to 600

# created new/updated scripts in
# /projects/users/caitlin/keats_x/new_sims folder

# all sims in new_sims/ folder have sample size of 600 = 75 iters of 8 person pedigree
# also use herit 

# qsub -q olga.q -p -50 -N pow9 -t 1-10000 batch_powerSim.sh
# which writes to powerSims_8ped_chr9_ld05_ar1/sim406_c01
# 20 variants, 8ped, auto=TRUE, effVar=c(0,0.4,0.6), c=0.1, sigx=siga=0.5
# cat sim*_res.txt >> allRes406_c01.txt
# wc -l allRes406_c01.txt ## should be 490000
# if so, rm sim406_c01*_res.txt

# powerSims_8ped_chr9_ld05_ar1/allRes406_c01.txt only has 9791 iters, not 10K

# ok, so all configurations for 8ped have been run, for chr 9 and x chr

# for null sims, just run powerSim with c=0 so all effect sizes are zero
# qsub -q olga.q -p -50 -N null -t 1-10000 batch_nullSim.sh


##
library(readr); library(dplyr)
library(tidyr); library(ggplot2)

setwd("/projects/users/caitlin/keats_x/new_sims/")

readPow <- function(fn,n=10000,prop,type,chr){
  dat <- read_delim(fn,delim=" ",col_names=FALSE,skip=1)
  dat$X1 <- NULL
  colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
  dat <- dat[!is.na(dat$model),]
  
  datNew <- dat %>%
    #mutate(iter=rep(1:n,each=12*4)) %>%
    mutate(finalpval=rep(c(rep(FALSE,11),TRUE),4*n)) %>%
    filter(finalpval==TRUE) %>%
    mutate(prop=prop) %>%
    mutate(type=type) %>%
    mutate(chr=chr)
  return(datNew)
}

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/allRes406_c02_s2x9_s2a6.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/allRes_604_c02_s2x9_s2a6.txt",prop="60/0/40",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/allRes208_c02_s2x9_s2a6.txt",prop="20/0/80",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/allRes406_c02_s2x9_s2a6.txt",prop="40/0/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/allRes604_c02_s2x9_s2a6.txt",prop="60/0/40",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/allRes208_c02_s2x9_s2a6.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allRes_c0_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims2 <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_8ped_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0,0.8)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

##### do the same for non-bt proportions now

chrX_8ped_ar1_15256 <- readPow("powerSims_8ped_chrX_ld05_ar1/allRes5356_c02_s2x9_s2a6.txt",prop="5/35/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8ped_chrX_ld05_ar1/allRes15256_c02_s2x9_s2a6.txt",prop="15/25/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8ped_chrX_ld05_ar1/allRes226_c02_s2x9_s2a6.txt",prop="20/20/60",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8ped_chr9_ld05_ar1/allRes5356_c02_s2x9_s2a6.txt",prop="5/35/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8ped_chr9_ld05_ar1/allRes15256_c02_s2x9_s2a6.txt",prop="15/25/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8ped_chr9_ld05_ar1/allRes226_c02_s2x9_s2a6.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allRes_c0_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_8ped_chr9_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims2 <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_8ped_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0,0.8)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

##### now do for 8pedFem pedigree config

chrX_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/allRes406_c02_s2x9_s2a6.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/allRes604_c02_s2x9_s2a6.txt",prop="60/0/40",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/allRes208_c02_s2x9_s2a6.txt",prop="20/0/80",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/allRes406_c02_s2x9_s2a6.txt",prop="40/0/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/allRes604_c02_s2x9_s2a6.txt",prop="60/0/40",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/allRes208_c02_s2x9_s2a6.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims2 <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_8pedFem_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0,0.8)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

##### do the same for non-bt proportions now

chrX_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/allRes5356_c02_s2x9_s2a6.txt",prop="5/35/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/allRes15256_c02_s2x9_s2a6.txt",prop="15/25/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_8pedFem_chrX_ld05_ar1/allRes226_c02_s2x9_s2a6.txt",prop="20/20/60",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/allRes5356_c02_s2x9_s2a6.txt",prop="5/35/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/allRes15256_c02_s2x9_s2a6.txt",prop="15/25/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_8pedFem_chr9_ld05_ar1/allRes226_c02_s2x9_s2a6.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_8pedFem_chrX_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_8pedFem_chr9_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims2 <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_8pedFem_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0,0.8)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()


##### finally do for set of unrelated samples

chrX_8ped_ar1_15256 <- readPow("powerSims_unrel_chrX_ld05_ar1/allRes406_c02_s2x9_s2a6.txt",prop="40/0/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_unrel_chrX_ld05_ar1/allRes604_c02_s2x9_s2a6.txt",prop="60/0/40",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_unrel_chrX_ld05_ar1/allRes208_c02_s2x9_s2a6.txt",prop="20/0/80",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_unrel_chr9_ld05_ar1/allRes406_c02_s2x9_s2a6.txt",prop="40/0/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_unrel_chr9_ld05_ar1/allRes604_c02_s2x9_s2a6.txt",prop="60/0/40",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_unrel_chr9_ld05_ar1/allRes208_c02_s2x9_s2a6.txt",prop="20/0/80",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "60/0/40"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/0/80"

dat <- read_delim("nullSims_unrel_chr9_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims2 <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="40/0/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "60/0/40"
nullSims4 <- nullSims2
nullSims4$prop <- "20/0/80"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("20/0/80","40/0/60","60/0/40"))

# plot these results now
pdf("power_unrel_btProps_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0,0.8)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

##### do the same for non-bt proportions now

chrX_8ped_ar1_15256 <- readPow("powerSims_unrel_chrX_ld05_ar1/allRes5356_c02_s2x9_s2a6.txt",prop="5/35/60",
                               type=TRUE,chr="X")
chrX_8ped_ar1_5356 <- readPow("powerSims_unrel_chrX_ld05_ar1/allRes15256_c02_s2x9_s2a6.txt",prop="15/25/60",
                              type=TRUE,chr="X")
chrX_8ped_ar1_226 <- readPow("powerSims_unrel_chrX_ld05_ar1/allRes226_c02_s2x9_s2a6.txt",prop="20/20/60",
                             type=TRUE,chr="X")
chr9_8ped_ar1_15256 <- readPow("powerSims_unrel_chr9_ld05_ar1/allRes5356_c02_s2x9_s2a6.txt",prop="5/35/60",
                               type=TRUE,chr="Autosomal")
chr9_8ped_ar1_5356 <- readPow("powerSims_unrel_chr9_ld05_ar1/allRes15256_c02_s2x9_s2a6.txt",prop="15/25/60",
                              type=TRUE,chr="Autosomal")
chr9_8ped_ar1_226 <- readPow("powerSims_unrel_chr9_ld05_ar1/allRes226_c02_s2x9_s2a6.txt",prop="20/20/60",
                             type=TRUE,chr="Autosomal")

allPow <- rbind(chrX_8ped_ar1_15256,chrX_8ped_ar1_5356,chrX_8ped_ar1_226,
                chr9_8ped_ar1_15256,chr9_8ped_ar1_5356,chr9_8ped_ar1_226)

# now read in null simulations
dat <- read_delim("nullSims_unrel_chrX_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="X") 

nullSimsa <- nullSims
nullSimsa$prop <- "15/25/60"
nullSimsb <- nullSimsa
nullSimsb$prop <- "20/20/60"

dat <- read_delim("nullSims_unrel_chr9_ld05_ar1/allRes_c00_s2x9_s2a6.txt",delim=" ",skip=1,col_names=FALSE)
dat$X1 <- NULL
colnames(dat) <- c("qstat","pvalue","rho","qmin","model")
nullSims2 <- dat %>%
  filter(qstat==0) %>%
  filter(model!="model") %>%
  mutate(finalpval=NA) %>%
  mutate(prop="5/35/60") %>%
  mutate(type=FALSE) %>% 
  mutate(chr="Autosomal") 

nullSims3 <- nullSims2
nullSims3$prop <- "15/25/60"
nullSims4 <- nullSims2
nullSims4$prop <- "20/20/60"

allRes <- rbind(allPow,nullSims,nullSims2,nullSimsa,nullSimsb,nullSims3,nullSims4)

pow <- allRes %>%
  group_by(model,prop,chr,type)
data.frame(summarize(pow,n())) # 10K for each model at each prop

alpha <- seq(from=1e-10,to=0.25,by=0.0001)
plo <- matrix(NA,nrow=48,ncol=length(alpha))
#colnames(plo) <- paste0("alpha.",alpha)
for(i in 1:length(alpha)){
  s <- summarize(pow,sum(pvalue<alpha[i])/n())
  plo[,i] <- data.frame(s)[,5]
}

plo <- data.frame(plo)
plo <- cbind(data.frame(s[,1:4]),plo)
plo$model[plo$model=="auto"] <- "MONSTER"
plo$model[plo$model=="both"] <- "KEATS-O"
plo$model[plo$model=="x"] <- "KEATS-O X only"
plo$model[plo$model=="unrel"] <- "SKATO"

allres <- plo %>%
  gather(alpha,rate,-c(model,type,prop,chr))

allres$type[allres$type==TRUE] <- "power"
allres$type[allres$type==FALSE] <- "null"
finalR <- allres %>%
  spread(type,rate) 
finalR$prop <- ordered(finalR$prop,levels=c("5/35/60","15/25/60","20/20/60"))

# plot these results now
pdf("power_unrel_trunc.pdf",width=14)
ggplot(finalR,aes(x=null,y=power,color=model)) + facet_grid(chr~prop) + 
  geom_line(size=1) + theme_bw() + scale_y_continuous(limits=c(0,0.8)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  scale_x_continuous(limits=c(0,0.1)) + xlab("False Positive Rate") + ylab("True Positive Rate")
dev.off()

rm(list=ls())


#####
# 47. Make heatmaps of LD between variants in chr 16, x chr RBC gene hits

library(GWASTools)
library(OLGApipeline); library(OLGAanalysis)
library(SNPRelate)
library(reshape2); library(ggplot2)
# snpgdsLDpair, method="corr"

setwd("/projects/users/caitlin/keats_x/keats_optimal")

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=TRUE,as.is=TRUE)
dim(geneAnnot) # 16502 13
head(geneAnnot)

hits <- geneAnnot[!is.na(geneAnnot$p.value)&-log10(geneAnnot$p.value)>6,]

# get the mapping of snps to genes
map16 <- get(load("olga_rbc_assoc/snpIDs_genes_chr16.RData"))
map16 <- sapply(map16,function(x){unique(x)})

hitsMap <- map16[is.element(names(map16),hits$geneID)]
snpRes <- data.frame(unlist(hitsMap))
len <- sapply(hitsMap,function(x){length(x)})
snpRes$gene <- rep(names(len),len)
snpRes <- merge(snpRes,hits,by.x="gene",by.y="geneID",all.x=TRUE)
colnames(snpRes)[2] <- "snpID"

# get the snp annot
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=16) 

snpRes <- merge(snpRes,pData(snpAnnot)[,c("rsID","snpID","position","chromosome")],by="snpID",all.x=TRUE)
dim(snpRes) # 95 17

snpRes$chromosome <- snpRes$chromosome.x
snpRes$chromosome[is.na(snpRes$chromosome)] <- snpRes$chromosome.y[is.na(snpRes$chromosome)]
snpRes$chromosome.x <- snpRes$chromosome.y <- NULL

# calculate pairwise LD between these variants
gdsOlga <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1/SOL_freeze1_chr-16.gds")

toCalc <- sort(unique(snpRes$snpID))
snpIDs <- read.gdsn(index.gdsn(gdsOlga,"snp.id"))
which(snpIDs==toCalc[1]) # 1767
which(snpIDs==toCalc[82]) # 5857, so need to read 5857-1767=4090
snp1 <- read.gdsn(index.gdsn(gdsOlga,"genotype"),start=c(which(snpIDs==toCalc[1]),1),count=c(4091,-1))
rownames(snp1) <- snpIDs[which(snpIDs==toCalc[1]):(which(snpIDs==toCalc[1])+4090)] 
dim(snp1) # 4091 x 849776 ie snp x sample

# loop through the snps we want pairwise LD for and calculate it
luc7lvsitfg3 <- matrix(NA,nrow=39,ncol=13)
luc7lSNPs <- snpRes$snpID[snpRes$gene==55692]
itfg3SNPs <- snpRes$snpID[snpRes$gene==83986]
for(i in 1:39){
  snpa <- snp1[rownames(snp1)==luc7lSNPs[i]]
  for(j in 1:13){
    snp2 <- snp1[rownames(snp1)==itfg3SNPs[j]]
    luc7lvsitfg3[i,j] <- snpgdsLDpair(snpa,snp2,method="corr")
  }
}
save(luc7lvsitfg3,file="luc7lvsitfg3_ld.RData")


luc7lvsrab <- matrix(NA,nrow=39,ncol=43)
rabSNPs <- snpRes$snpID[snpRes$gene==9727]
for(i in 1:39){
  snpa <- snp1[rownames(snp1)==luc7lSNPs[i]]
  for(j in 1:43){
    snp2 <- snp1[rownames(snp1)==rabSNPs[j]]
    luc7lvsrab[i,j] <- snpgdsLDpair(snpa,snp2,method="corr")
  }
}
save(luc7lvsrab,file="luc7lvsrab_ld.RData")


rabvsitfg3 <- matrix(NA,nrow=43,ncol=13)
for(i in 1:43){
  snpa <- snp1[rownames(snp1)==rabSNPs[i]]
  for(j in 1:13){
    snp2 <- snp1[rownames(snp1)==itfg3SNPs[j]]
    rabvsitfg3[i,j] <- snpgdsLDpair(snpa,snp2,method="corr")
  }
}
save(rabvsitfg3,file="rabvsitfg3.RData")

library(ggplot2); library(reshape2)
snpRes <- snpRes[order(snpRes$snpID),]
luc7lrs <- snpRes$rsID[snpRes$gene==55692]
itfg3rs <- snpRes$rsID[snpRes$gene==83986]
rownames(luc7lvsitfg3) <- luc7lrs
colnames(luc7lvsitfg3) <- itfg3rs
mlt <- melt(luc7lvsitfg3)
colnames(mlt) <- c("LUC7L","ITFG3","LD")
pdf("luc7lvsitfg3_ld.pdf")
ggplot(mlt,aes(x=LUC7L,y=ITFG3,fill=LD)) + geom_tile(color="white") +
  scale_fill_gradient(low="white",high="steelblue")
dev.off()

rabrs <- snpRes$rsID[snpRes$gene==9727]
rownames(luc7lvsrab) <- luc7lrs
colnames(luc7lvsrab) <- rabrs
mlt <- melt(luc7lvsrab)
colnames(mlt) <- c("LUC7L","RAB11FIP3","LD")
pdf("luc7lvsrab_ld.pdf")
ggplot(mlt,aes(x=LUC7L,y=RAB11FIP3,fill=LD)) + geom_tile(color="white") +
  scale_fill_gradient(low="white",high="steelblue")
dev.off()

rownames(rabvsitfg3) <- rabrs
colnames(rabvsitfg3) <- itfg3rs
mlt <- melt(rabvsitfg3)
colnames(mlt) <- c("RAB11FIP3","ITFG3","LD")
pdf("rabvsitfg3_ld.pdf")
ggplot(mlt,aes(y=ITFG3,x=RAB11FIP3,fill=LD)) + geom_tile(color="white") +
  scale_fill_gradient(low="white",high="steelblue")
dev.off()



## ok, do this for the chr 23 hits too
# get the mapping of snps to genes
map23 <- get(load("olga_rbc_assoc/snpIDs_genes_chr23.RData"))
map23 <- sapply(map23,function(x){unique(x)})

hitsMap <- map23[is.element(names(map23),hits$geneID)]
snpRes2 <- data.frame(unlist(hitsMap))
len <- sapply(hitsMap,function(x){length(x)})
snpRes2$gene <- rep(names(len),len)
snpRes2 <- merge(snpRes2,hits,by.x="gene",by.y="geneID",all.x=TRUE)
colnames(snpRes2)[2] <- "snpID"

# get the snp annot
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=23) 

snpRes2 <- merge(snpRes2,pData(snpAnnot)[,c("rsID","snpID","position","chromosome")],by="snpID",all.x=TRUE)
dim(snpRes2) # 20 17

snpRes2$chromosome <- snpRes2$chromosome.x
snpRes2$chromosome.x <- snpRes2$chromosome.y <- NULL

gdsOlga <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1/SOL_freeze1_chr-23.gds")

toCalc <- sort(unique(snpRes2$snpID))
snpIDs <- read.gdsn(index.gdsn(gdsOlga,"snp.id"))
which(snpIDs==toCalc[1]) # 931874
which(snpIDs==toCalc[20]) # 933653, so need to read 1779
sid <- which(snpIDs==toCalc[1])
snp1 <- read.gdsn(index.gdsn(gdsOlga,"genotype"),start=c(1,sid),count=c(-1,1780))
colnames(snp1) <- snpIDs[which(snpIDs==toCalc[1]):(which(snpIDs==toCalc[1])+1779)] 
dim(snp1) # 4091 x 849776 ie snp x sample

# loop through the snps we want pairwise LD for and calculate it
g6pdvsf8 <- matrix(NA,nrow=15,ncol=5)
g6pdSNPs <- snpRes2$snpID[snpRes2$gene==2539]
f8SNPs <- snpRes2$snpID[snpRes2$gene==2157]
for(i in 1:15){
  snpa <- snp1[colnames(snp1)==f8SNPs[i]]
  for(j in 1:5){
    snp2 <- snp1[colnames(snp1)==g6pdSNPs[j]]
    g6pdvsf8[i,j] <- snpgdsLDpair(snpa,snp2,method="corr")
  }
}
save(g6pdvsf8,file="g6pdvsf8_ld.RData")


snpRes2 <- snpRes2[order(snpRes2$snpID),]
g6rs <- snpRes2$rsID[snpRes2$gene==2539]
f8rs <- snpRes2$rsID[snpRes2$gene==2157]
rownames(g6pdvsf8) <- f8rs
colnames(g6pdvsf8) <- g6rs
mlt <- melt(g6pdvsf8)
colnames(mlt) <- c("F8","G6PD","LD")
pdf("g6pdvsf8_ld.pdf")
ggplot(mlt,aes(x=G6PD,y=F8,fill=LD)) + geom_tile(color="white") +
  scale_fill_gradient(low="white",high="steelblue")
dev.off()


rm(list=ls())


#####
# 48. Find GATS and AHI1 rare gene analysis p-values

# GATS is entrez id 352954, on chr 7
# AHI1 is entrez id 54806, on chr 6

library(GWASTools)
library(OLGApipeline); library(OLGAanalysis)

setwd("/projects/users/caitlin/keats_x/keats_optimal")

geneAnnot <- read.table("olga_rbc_assoc/geneList_genomewide_results.txt",header=TRUE,as.is=TRUE)
dim(geneAnnot) # 16502 13
head(geneAnnot)

geneAnnot[geneAnnot$geneID==352594,]
# chr 7, pvalue 0.840, 11 vars, rho=0

geneAnnot[geneAnnot$geneID==54806,]
# chr 6, pvalue 0.522, 62 vars, rho=0

rm(list=ls())


#####