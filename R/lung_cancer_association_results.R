library(data.table)
library(dplyr)

###

model <- fread('/home/labs/dnalab/share/lims/xduran/Lung/data.model')
logistic <- fread('/home/labs/dnalab/share/lims/xduran/Lung/data.assoc.logistic')
adjusted <- fread('/home/labs/dnalab/share/lims/xduran/Lung/data.assoc.logistic.adjusted')

lung_bim <- read.table('inst/extdata/lung/chr12_imputed.bim', sep = '\t', stringsAsFactors = FALSE)

View(model[1:10000])

model.p <- model %>%
  filter(
    P < 0.001
  )

View(lung_bim[1:100,])

snp.split <- do.call(rbind, strsplit(model.p$SNP, ':'))

###


path = "/home/labs/dnalab/igalvan/Documents/dnalabuser1/gpfs/work/pr94bu/di73ceq/IMPPC_Lung/GWImp_COMPSs/test_progres_new_compss/old_gwimpcompss/outputs/associations_covTreatment/IMPPC_Lung_progres_for_1kgphase3/summary/"

setwd(path)

assoc_results = fread("filtered_by_all_results_1kgphase3_chr_1_to_22.txt")

assoc_results$OR = exp(assoc_results$frequentist_add_beta_1)

assoc_results$chr_position = do.call("paste",list(assoc_results$chr,assoc_results$position,sep=":")) 


assoc_results_filtered = assoc_results %>% filter(info>0.7) %>% filter(OR<10 & OR>0.1) %>% 
  filter(frequentist_add_pvalue<1e-3) %>% select(rs_id_all,chr,position,frequentist_add_pvalue)


# afegim els p-values corregits 


assoc_results_corrected = fread("corrected_pvalues_1kgphase3_chr_1_to_22.txt")


assoc_results_corrected$chr_position = do.call("paste",list(assoc_results_corrected$CHR,
                                                            assoc_results_corrected$BP,
                                                            sep=":")) 

colnames(assoc_results_corrected)[3] = "P_adjusted"


assoc_results_filtered_corrected = left_join(assoc_results_filtered,assoc_results_corrected %>% 
                                    dplyr::select(chr_position,P_adjusted))



head(assoc_results_filtered_corrected)
