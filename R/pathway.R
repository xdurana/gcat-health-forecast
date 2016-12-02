options(java.parameters = "-Xmx8000m")

library(seq2pathway)
library(dplyr)
library(tidyr)
library(xlsx)

get_vcf <- function() {
  
  lung <- read.xlsx2('inst/extdata/lung/genomic_general_features_info.xlsx', sheetIndex = 1)
  snp_dat <- describe_snp_bim(lung)
  
  write.csv(snp_dat, 'output/lung/snp_dat.csv', row.names = FALSE, quote = FALSE)

  vcf <- snp_dat %>%
    mutate(
      CHROM=chrom,
      POS=chromStart,
      ID=paste(chrom, chromStart, sep = ':'),
      REF=ref,
      ALT=alt,
      QUAL='.',
      FILTER='.',
      INFO='.'
    ) %>%
    dplyr::select(
      CHROM,
      POS,
      ID,
      REF,
      ALT,
      QUAL,
      FILTER,
      INFO
    )
  
  write.table(vcf, 'output/lung/lung.vcf', col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ' ')
  vcf
}

get_vep <- function() {
  
  vep <- read.table("inst/extdata/lung/WDWXwGbATmbKuqpM.txt", sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  
  snp_dat_select <- snp_dat %>%
    mutate(
      Uploaded_variation=paste(chrom, chromStart, sep = ':')
    ) %>%
    dplyr::select(
      Uploaded_variation,
      stability,
      scaled_coefficients
    )
  
  vep <- vep %>%
    merge(
      snp_dat_select,
      all.x = TRUE
    ) %>%
    arrange(
      desc(scaled_coefficients),
      desc(stability)
    )
  
  vep_dat <- vep %>%
    mutate(
      name=Existing_variation,
      chrom=as.integer(as.character(do.call(rbind, strsplit(Uploaded_variation, ':'))[,1])),
      chromStart=as.integer(as.character(do.call(rbind, strsplit(Uploaded_variation, ':'))[,2])),
      chromEnd=chromStart+1
    ) %>%
    dplyr::select(
      name,
      chrom,
      chromStart,
      chromEnd
    ) %>%
    unique()
  
  gene_dat <- runseq2gene(
    inputfile = snp_dat,
    genome = "hg19",
    adjacent = TRUE,
    SNP = TRUE,
    search_radius = 150000,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )
  
  genes <- gene_dat$seq2gene_FullResult
  write.csv2(genes, 'output/lung/genes.csv', row.names = FALSE, quote = TRUE)
  
  gene_list <- as.data.frame(table(genes$gene_name))
  write.csv2(gene_list, 'output/lung/gene_list.csv', row.names = FALSE, quote = TRUE)
  
    
  ### PATHWAY
  
  path_dat <- runseq2pathway(
    inputfile = snp_dat,
    genome = "hg19",
    adjacent = FALSE,
    SNP = TRUE,
    search_radius = 150000,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )
  
  genes <- pat_dat$seq2gene_result$seq2gene_FullResult %>%
    merge(snp_dat)
  
  pathways <- pat_dat$gene2pathway_result.FET %>%
    with(
      rbind(
        GO_BP,
        GO_CC,
        GO_MF
      )
    )

  write.xlsx2(vep, 'output/lung/lung_vep.xlsx')
  write.xlsx2(genes, 'output/lung/lung_genes.xlsx')
  write.xlsx2(pathways, 'output/lung/lung_pathways.xlsx')
}

#' @title Prepare for runseq2pathway
#' @param lung
describe_snp_bim <- function(lung) {
  
  lung_bim <- read.table('inst/extdata/lung/chr12_imputed.bim', sep = '\t', stringsAsFactors = FALSE)

  lung.pos <- lung
  lung.pos$name <- lung.pos$names
  lung.pos$name <- gsub("b'", "", as.character(lung.pos$name))
  lung.pos$name <- do.call(rbind, strsplit(lung.pos$name, '_'))[,1]
  
  lung.pos <- lung.pos %>%
    merge(lung_bim, all.x = TRUE, by.x = 'name', by.y = 'V2')
  
  snp_dat <- lung.pos %>%
    transform(
      chrom=12,
      chromStart=as.integer(V4),
      chromEnd=as.integer(V4) + 1,
      stability=as.integer(as.character(stability)),
      ref=V5,
      alt=V6
    ) %>%
    dplyr::select(
      name,
      chrom,
      chromStart,
      chromEnd,
      stability,
      scaled_coefficients,
      ref,
      alt
    ) %>%
    filter(
      stability > 0
    )
  
  snp_dat
}

top20 <- function(vep, n) {
  
  vep_unique <- vep %>%
    dplyr::select(
      Uploaded_variation,
      stability,
      scaled_coefficients,
      SYMBOL,
      Gene
    ) %>%
    unique()
  
  df <- (as.data.frame(table(vep_unique[1:n,]$SYMBOL)) %>% arrange(desc(Freq)))[1:20,]
  df
}