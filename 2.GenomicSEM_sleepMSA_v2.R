library(devtools)
library(GenomicSEM)

set.seed(100)
dir_forGWAS = "~/genomicSEM/v2_sleepMSA/GWAS/uncompressed/"
dir_forSEM = "~/genomicSEM/v2_sleepMSA/sem/"


load("~/genomicSEM/v2_sleepMSA/sem/All_data_NO_GWAS.RData")


total_rows <- nrow(p_sumstats)
chunk_size <- 10000
num_chunks <- ceiling(total_rows / chunk_size)
for (i in 1:num_chunks) {
        if(i<num_chunks){
                if(i==1){
                        end_index <- i*chunk_size
                        chunk <- p_sumstats[i:end_index, ]
                        
                        save(chunk, file=paste0("~/genomicSEM/v2_sleepMSA/sem/chunks/chunk_", i,".RData", collapse = ""))   
                }else{
                        end_index <- i*chunk_size
                        ini_index = (i-1)*10000+1
                        chunk <- p_sumstats[ini_index:end_index, ]
                        
                        save(chunk, file=paste0("~/genomicSEM/v2_sleepMSA/sem/chunks/chunk_", i,".RData", collapse = ""))
                }    
        }else{
                end_index <- nrow(p_sumstats)
                ini_index = (i-1)*10000+1
                chunk <- p_sumstats[ini_index:end_index, ]
                save(chunk, file=paste0("~/genomicSEM/v2_sleepMSA/sem/chunks/chunk_", i,".RData", collapse = ""))     
        }
        
}
