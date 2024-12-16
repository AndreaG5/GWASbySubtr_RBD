library(devtools)
library(GenomicSEM)

set.seed(100)
dir_forGWAS = "~/genomicSEM/v2_sleepMSA/GWAS/uncompressed/"
dir_forSEM = "~/genomicSEM/v2_sleepMSA/sem/"
setwd(dir_forSEM)

load("~/genomicSEM/v2_sleepMSA/sem/All_data_NO_GWAS.RData")


##########################
#user model
##########################

mod_model<-'F1=~ NA*RBD 
F2=~NA*RBD+PD +DLB + MSA 
F1~~1*F1
F2~~1*F2
F1~~0.1*F2
RBD~~0*RBD
a>.001
DLB~~a*DLB'


usr_model = usermodel(LDSCoutput, model=mod_model, estimation = "DWLS", CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

model<-'F1=~ NA*RBD
F2=~NA*RBD+PD +DLB + MSA 
F1~~1*F1
F2~~1*F2
RBD~~0*RBD
a>.001
DLB~~a*DLB

F1~SNP
F2~SNP

F1~~0.1*F2

SNP~~SNP'

total_rows <- nrow(p_sumstats)
chunk_size <- 10000
num_chunks <- ceiling(total_rows / chunk_size)
for (i in 1:num_chunks) {
        files_li = list.files("~/genomicSEM/v2_sleepMSA/sem/chunks/", pattern="RData", full.names = T)
        file_index = grep(paste0("_",i, ".RData$", collapse=""), files_li)
        
        load(files_li[file_index])
        res = userGWAS(covstruc=LDSCoutput,SNPs=chunk,estimation="DWLS",model=model,sub =c("F1~SNP","F2~SNP"),
                       Q_SNP = T, parallel= TRUE,smooth_check=TRUE, toler = 1e-26, cores=4)
        path_name = paste0("~/genomicSEM/v2_sleepMSA/sem/chunks_GWAS01/chunk_GWAS_", i, ".RData", collapse = "")
        txt_name = paste0("~/genomicSEM/v2_sleepMSA/sem/txt01/chunk_GWAS_", i, ".txt", collapse = "")
        save(res, file = path_name)
        write.table(res, txt_name, quote = F, sep = "\t", col.names = T, row.names = T)   
        
        
}

