library(devtools)
library(GenomicSEM)

set.seed(100)
dir_forGWAS = "~/genomicSEM/v2_sleepMSA/GWAS/uncompressed/"
dir_forSEM = "~/genomicSEM/v2_sleepMSA/sem/"

####################################
########### MUNGE DATA #############
####################################

files=c(paste0(dir_forGWAS,"RBD_rsid_fin_Neffcorrected.tab",collapse = ""),paste0(dir_forGWAS,"nallsEtAl2019_rsid_fin_Neffcorrected.tab",collapse = ""),
        paste0(dir_forGWAS,"DLB_rsid_fin_Neffcorrected.tab",collapse = ""),
        paste0(dir_forGWAS,"MSA_rsid_fin_Neffcorrected.tab",collapse = ""))

hm3<-"~/genomicSEM/w_hm3.snplist"
trait.names<-c("RBD", "PD", "DLB", "MSA")

setwd(dir_forSEM)

munge(files=files,hm3=hm3,trait.names=trait.names)

####################################
############# LDSC #################
####################################

traits<-c(paste0(dir_forSEM,"RBD.sumstats.gz",collapse = ""),paste0(dir_forSEM,"PD.sumstats.gz",collapse = ""),
          paste0(dir_forSEM,"DLB.sumstats.gz",collapse = ""),paste0(dir_forSEM,"MSA.sumstats.gz",collapse = ""))

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(.5,.5,.5,.5)
#vector of population prevalences
population.prev<-c(.0125,.005,.01,.000072)
#the folder of LD scores
ld<-"~/genomicSEM/eur_w_ld_chr/"
#the folder of LD weights [typically the same as folder of LD scores]
wld<-"~/genomicSEM/eur_w_ld_chr/"
#name the traits
trait.names<-c("RBD","PD", "DLB", "MSA")
#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

mod_model = 'F1=~ NA*PD +ALS+DLB
F1~~1*F1
DLB ~~ 0*DLB'

#usr_model = usermodel(LDSCoutput, model=mod_model, estimation = "DWLS", CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

####################################
######## Sumstats of model #########
####################################

files=c(paste0(dir_forGWAS,"RBD_rsid_fin_Neffcorrected.tab",collapse = ""),paste0(dir_forGWAS,"nallsEtAl2019_rsid_fin_Neffcorrected.tab",collapse = ""),
        paste0(dir_forGWAS,"DLB_rsid_fin_Neffcorrected.tab",collapse = ""),
        paste0(dir_forGWAS,"MSA_rsid_fin_Neffcorrected.tab",collapse = ""))

ref="~/genomicSEM/reference.1000G.maf.0.005.txt.gz"
trait.names=c("RBD","PD", "DLB", "MSA")
se.logit=c(T,T,F,F)
OLS=c(F,F,F,F)
linprob=c(F,F,F,F)
p_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,parallel=T)


save(LDSCoutput, p_sumstats, CommonFactor_DWLS, file="~/genomicSEM/v2_sleepMSA/sem/All_data_NO_GWAS.RData")
