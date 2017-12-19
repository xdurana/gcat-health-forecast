library(data.table)
library(dplyr)

source('R/pathway.R')

gwas_pathways <- function(gwas.assoc, ds.name) {
  
  gwas.assoc.filtered <- gwas.assoc %>%
    filter(
      P_adjusted < 10^-4 &
        all_maf > 0.01 &
        OR < 20
    ) %>%
    rename(
      name=chr_position,
      chrom=CHR,
      chromStart=BP,
      ref=alleleA,
      alt=alleleB
    ) %>%
    transform(
      chromEnd=chromStart + 1
    ) %>%
    select(
      name,
      chrom,
      chromStart,
      chromEnd,
      P_adjusted,
      all_maf,
      OR,
      ref,
      alt
    )
  
  report_pathway(gwas.assoc.filtered, ds.name)
}

run <- function() {
  
  gwas_pathways(fread('inst/extdata/lung/assoc_results_whole_genome_impute_07.txt'), 'gwas_impute_07')
  gwas_pathways(fread('inst/extdata/lung/assoc_results_whole_genome.txt'), 'gwas')
}