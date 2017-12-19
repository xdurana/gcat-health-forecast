options(java.parameters = "-Xmx8000m")

library(seq2pathway)
library(qqman)
library(dplyr)
library(tidyr)
library(xlsx)
library(Vennerable)
library(data.table)

hf_to_snp_dat <- function(chr) {

  print(sprintf("Transforming chr%s", chr))
  
  directory <- sprintf('inst/extdata/lung/lung.cancer_results/chr%s', chr)
  variants <- fread(file.path(directory, 'general_features_info.csv'))
  variants_filtered <- variants %>%
    filter(stability > 5) %>%
    transform(
      chr=chr
    )
  
  lung_bim <- fread(sprintf('/home/labs/dnalab/share/lims/xduran/Lung/all_genome/chr%s/chr%s_imputed.bim', chr, chr))
  
  lung.pos <- variants_filtered
  lung.pos$name <- lung.pos$names
  lung.pos$name <- gsub("b'", "", as.character(lung.pos$name))
  lung.pos$name <- do.call(rbind, strsplit(lung.pos$name, '_'))[,1]
  
  lung.pos <- lung.pos %>%
    merge(lung_bim, all.x = TRUE, by.x = 'name', by.y = 'V2')

  snp_dat <- lung.pos %>%
    transform(
      chrom=chr,
      chromStart=as.integer(V4),
      chromEnd=as.integer(V4) + 1,
      stability=as.integer(as.character(stability)),
      scaled_importances='.',
      ref=V5,
      alt=V6
    ) %>%
    dplyr::select(
      name,
      chrom,
      chromStart,
      chromEnd,
      stability,
      scaled_importances,
      ref,
      alt
    ) %>%
    filter(
      stability > 0
    )
  
  print(sprintf("Transforming chr%s stability [%s, %s])", chr, min(variants$stability), max(variants$stability)))
  print(sprintf("Transforming chr%s ok", chr))
  
  snp_dat  
}

function() {
  variants <- do.call(rbind, lapply(c(1:22), hf_to_snp_dat))
  variants_50 <- variants %>%
    filter(stability >= 40)
  report_pathway(variants_50, 'hf')
}