library(devtools)
library(GenomicSEM)

set.seed(100)
setwd("~/genomicSEM/v2_sleepMSA/sem/")

dir_forchunk_GWAS = "~/genomicSEM/v2_sleepMSA/sem/chunks_GWAS/"
dir_forTXT = "~/genomicSEM/v2_sleepMSA/sem/txt/"

load("~/genomicSEM/v2_sleepMSA/sem/All_data_NO_GWAS.RData")


output = list.files(dir_forTXT, pattern = ".txt", full.names = T)

##### merge outputs
mGWAS = purrr::map_dfr(output, function(x){
        print(x)
        df = read.delim(x)
        df$warning = as.numeric(df$warning)
        df
})


#### subset first 25 columns are relative to F1~SNP, trailing 25 are for F2~SNP
mGWAS_F1 = mGWAS[,1:25]
mGWAS_F2 = mGWAS[,26:ncol(mGWAS)]

colnames(mGWAS_F2) = colnames(mGWAS_F1)

save(mGWAS_F1, file="~/genomicSEM/v2_sleepMSA/sem/mGWAS_F1xSNP.RData")
save(mGWAS_F2, file="~/genomicSEM/v2_sleepMSA/sem/mGWAS_F2xSNP.RData")

####### reorder based on sumstats input
all(mGWAS_F1$SNP %in% p_sumstats$SNP)
mGWAS_F1 = mGWAS_F1[match(p_sumstats$SNP,mGWAS_F1$SNP),]
####### reorder based on sumstats input
all(mGWAS_F2$SNP %in% p_sumstats$SNP)
mGWAS_F2 = mGWAS_F2[match(p_sumstats$SNP,mGWAS_F2$SNP),]

#### remove SNPs with lavaan warning 
mGWAS_F1 = mGWAS_F1[!is.na(mGWAS_F1$Z_Estimate),]
#### remove SNPs with lavaan warning 
mGWAS_F2 = mGWAS_F2[!is.na(mGWAS_F2$Z_Estimate),]



###########################################
####### calculate3 effective sample size###
###########################################
#F1 
# Because of the Cholesky model, we need to adjust with the residual heritabilities 
# (the estimates for the path loadings between the latent variable and the measured variable). 
# I use here the sqrt of the heritability(=Unstandardized estimate path loadings, i.e. first column of modelfit) 

## change 0.1993 with Unstand_EST from model output
mGWAS_F1$Neff <- ((mGWAS_F1$Z_Estimate/(mGWAS_F1$est*sqrt(0.2852*0.2852)))^2)/(2*mGWAS_F1$MAF*(1-mGWAS_F1$MAF)) 

#F2
#h2 of single traits in array. rg correlation matrix (cov2cor(LDSC$S), lambda -> array of factor loadings)
#order is --> RBD,PD,DLB,MSA
h2 = c(.1412,.0088,.1043,0.0397)
rg = cov2cor(LDSCoutput$S)
lambda = c(.76,.33,.99,.50)
numerator = sum(lambda^2 * h2) + 2*sum(upper.tri(rg) * outer(lambda, lambda)*rg)
denominator = sum(lambda^2) + 2*sum(upper.tri(outer(lambda, lambda)))
h2_F2 = numerator/denominator
## 0.2456 is h2_F2

mGWAS_F2$Neff <- ((mGWAS_F2$Z_Estimate/(mGWAS_F2$est*sqrt(.2456)))^2)/(2*mGWAS_F2$MAF*(1-mGWAS_F2$MAF))   #0.446483*0.446483=0.1993
summary(mGWAS_F2$Neff)



##### rename columns for LDscore
colnames(mGWAS_F1)[c(12,15)] = c("BETA","P")
gz1 <- gzfile("~/genomicSEM/v2_sleepMSA/sem/mGWAS_All_F1_sumstats.gz", "w")
write.table(mGWAS_F1, gz1, sep="\t", row.names = F,quote = F)
close(gz1)
##### rename columns for LDscore
colnames(mGWAS_F2)[c(12,15)] = c("BETA","P")
gz1 <- gzfile("~/genomicSEM/v2_sleepMSA/sem/mGWAS_All_F2_sumstats.gz", "w")
write.table(mGWAS_F2, gz1, sep="\t", row.names = F,quote = F)
close(gz1)