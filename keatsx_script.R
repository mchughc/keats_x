
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
# 12. Implement burden tests for sims, too



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
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- filter(dat,!is.na(stat))

chr9_8ped_ar1 <- dat[!duplicated(dat$famSKAT),]
chr9_8ped_ar1 <- chr9_8ped_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  filter(famSKAT<=1) %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8ped_chr9ceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- filter(dat,stat=="pval")
chr9_8ped_ceu <- dat[!duplicated(dat$famSKAT),]
chr9_8ped_ceu <- chr9_8ped_ceu %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8ped_chr9mxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- filter(dat,stat=="pval")
chr9_8ped_mxl <- dat[!duplicated(dat$famSKAT),]
chr9_8ped_mxl <- chr9_8ped_mxl %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8ped") %>%
  mutate(ld="MXL")

dat <- read_delim("nullSims_8pedfem_chr9_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- filter(dat,!is.na(stat))

chr9_8pedfem_ar1 <- dat[!duplicated(dat$famSKAT),]
chr9_8pedfem_ar1 <- chr9_8pedfem_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  filter(famSKAT<=1) %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8pedfem_chr9ceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- filter(dat,stat=="pval")
chr9_8pedfem_ceu <- dat[!duplicated(dat$famSKAT),]
chr9_8pedfem_ceu <- chr9_8pedfem_ceu %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8pedfem_chr9mxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- filter(dat,stat=="pval")
chr9_8pedfem_mxl <- dat[!duplicated(dat$famSKAT),]
chr9_8pedfem_mxl <- chr9_8pedfem_mxl %>%
  mutate(chr="autosomal") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="MXL")

dat <- read_delim("nullSims_8ped_chrX_ld05_ar1/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
dat <- filter(dat,!is.na(stat))

x_8ped_ar1 <- dat[!duplicated(dat$famSKAT),]
x_8ped_ar1 <- x_8ped_ar1 %>%
  mutate(famSKAT=as.numeric(famSKAT)) %>%
  mutate(XKC_only=as.numeric(XKC_only)) %>%
  mutate(KEATS_X=as.numeric(KEATS_X)) %>%
  filter(famSKAT<=1) %>%
  mutate(chr="X") %>%
  mutate(ped="8ped") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8ped_chrXceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS_X","stat")
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
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="AR1_0.5")

dat <- read_delim("nullSims_8pedfem_chrXceu/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- filter(dat,stat=="pval")
x_8pedfem_ceu <- dat[!duplicated(dat$famSKAT),]
x_8pedfem_ceu <- x_8pedfem_ceu %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="CEU")

dat <- read_delim("nullSims_8pedfem_chrXmxl/allRes.txt",
                  delim=" ",col_names=FALSE,skip=1)
dat <- mutate(dat,X1=NULL)
colnames(dat) <- c("famSKAT","XKC_only","KEATS-X","stat")
dat <- filter(dat,stat=="pval")
x_8pedfem_mxl <- dat[!duplicated(dat$famSKAT),]
x_8pedfem_mxl <- x_8pedfem_mxl %>%
  mutate(chr="X") %>%
  mutate(ped="8pedfem") %>%
  mutate(ld="MXL")

# merge all these together
colnames(chr9_8ped_ceu)[3] <- "KEATS_X" 
colnames(chr9_8ped_mxl)[3] <- "KEATS_X" 
colnames(chr9_8pedfem_ceu)[3] <- "KEATS_X" 
colnames(chr9_8pedfem_mxl)[3] <- "KEATS_X" 
colnames(x_8pedfem_ceu)[3] <- "KEATS_X"
colnames(x_8pedfem_mxl)[3] <- "KEATS_X"

tyI_group <- rbind(chr9_8ped_ar1,chr9_8ped_ceu,chr9_8ped_mxl,
                   chr9_8pedfem_ar1,chr9_8pedfem_ceu,chr9_8pedfem_mxl,
                   x_8ped_ar1,x_8ped_ceu,x_8ped_mxl,
                   x_8pedfem_ar1,x_8pedfem_ceu,x_8pedfem_mxl)

tyI_group <- tyI_group %>% 
  gather(model,pval,-c(stat,chr,ped,ld)) %>%
  mutate(pval=as.numeric(pval)) %>%
  filter(model!="XKC_only") %>%
  group_by(chr,ped,ld,model)

n <- summarize(tyI_group,n()) # 1M each, only 10k of the AR1 iters
res <- summarize(tyI_group,
          "alpha=0.01"=sum(pval<0.01)/n(),
          "alpha=0.001"=sum(pval<0.001)/n(),
          "alpha=0.0001"=sum(pval<0.0001)/n())

# save these as xtable
library(xtable)
xtable(res)

rm(list=ls())


