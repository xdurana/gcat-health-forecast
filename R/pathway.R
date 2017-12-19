options(java.parameters = "-Xmx8000m")

library(seq2pathway)
library(qqman)
library(dplyr)
library(tidyr)
library(xlsx)
library(Vennerable)
library(data.table)

get_vcf <- function() {
  
  lung <- read.xlsx2('inst/extdata/lung/genomic_general_features_info.xlsx', sheetIndex = 1)
  snp_dat <- describe_snp_bim(lung, '12')
  
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

  ### GET GENES
  
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
  write.csv2(genes, 'output/lung/genes_hf.csv', row.names = FALSE, quote = TRUE)
  
  gene_list <- as.data.frame(table(genes$gene_name))
  write.csv2(gene_list, 'output/lung/gene_hf_list.csv', row.names = FALSE, quote = TRUE)

  ### FILTER SNP
  
  snp_dat_filtered <- snp_dat[as.numeric(as.character(snp_dat$stability)) > 5,]
  #snp_dat_filtered <- snp_dat[as.numeric(as.character(snp_dat$scaled_coefficients)) > 0.45,]
  
  ### GET GENES
  
  gene_dat <- runseq2gene(
    inputfile = snp_dat_filtered,
    genome = "hg19",
    adjacent = TRUE,
    SNP = TRUE,
    search_radius = 150000,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )

  genes <- gene_dat$seq2gene_FullResult
  write.csv2(genes, 'output/lung/genes_hf_filtered.csv', row.names = FALSE, quote = TRUE)
  
  gene_list <- as.data.frame(table(genes$gene_name))
  write.csv2(gene_list, 'output/lung/gene_hf_filtered_list.csv', row.names = FALSE, quote = TRUE)
  
  genes$gene_name <- as.character(genes$gene_name)
  
  selected_variants <- genes[genes$gene_name %in% c('VWF', 'KRAS', 'KSR2', 'LUM', 'KERA'),]
  write.csv2(selected_variants, 'output/lung/selected_variants.csv', row.names = FALSE, quote = TRUE)

  ### PATHWAY
  
  path_dat <- runseq2pathway(
    inputfile = snp_dat_filtered,
    genome = "hg19",
    adjacent = FALSE,
    SNP = TRUE,
    search_radius = 150000,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )

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
describe_snp_bim <- function(lung, chr) {
  
  dir_bim <- sprintf('/home/labs/dnalab/share/lims/xduran/Lung/all_genome/chr%s', chr)
  lung_bim <- fread(dir_bim)
  
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

#' @title 
lung_gcat_all <- function() {
  
  all_snp <- read.csv('inst/extdata/lung/assoc_results_whole_genome.txt', sep = ' ', stringsAsFactors = FALSE)
  all_snp_dat <- all_snp %>%
    dplyr::rename(
      chrom=CHR,
      chromStart=BP
    ) %>%
    transform(
      name=sprintf("%s:%s:%s:%s", chrom, chromStart, alleleA, alleleB),
      chromEnd=chromStart + 1
    ) %>%
    dplyr::select(
      name,
      chrom,
      chromStart,
      chromEnd,
      P_adjusted
    )

  vcf <- all_snp %>%
    mutate(
      CHROM=CHR,
      POS=BP,
      ID=paste(CHR, POS, sep = ':'),
      REF=alleleA,
      ALT=alleleB,
      QUAL=P_adjusted,
      FILTER='PASS',
      INFO='.',
      FORMAT='.'
    ) %>%
    dplyr::select(
      CHROM,
      POS,
      ID,
      REF,
      ALT,
      QUAL,
      FILTER,
      INFO,
      FORMAT
    )  
  
  write.table(vcf, 'output/lung/lung.vcf', col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  all_snp_dat_4 <- all_snp_dat %>% filter(P_adjusted < 10^-4)
  all_snp_dat_5 <- all_snp_dat %>% filter(P_adjusted < 10^-5)

  report_pathway(all_snp_dat_4, 150000)
}

report_pathway <- function(snps, ds.name, radius=150000) {

  path_dat <- runseq2pathway(
    inputfile = snps,
    genome = "hg19",
    adjacent = FALSE,
    SNP = FALSE,
    search_radius = radius,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )
  
  genes <- path_dat$seq2gene_result$seq2gene_FullResult
  gene_list <- as.data.frame(table(genes$gene_name))
  genes$gene_name <- as.character(genes$gene_name)

  pathways <- path_dat$gene2pathway_result.FET %>%
    with(
      rbind(
        GO_BP,
        GO_CC,
        GO_MF
      )
    ) %>%
    arrange(
      FDR
    )
  
  directory <- sprintf('output/lung/%s', ds.name)
  
  write.csv2(genes, file.path(directory, 'genes.csv'), row.names = FALSE, quote = TRUE)
  write.csv2(gene_list, file.path(directory, 'gene_list.csv'), row.names = FALSE, quote = TRUE)
  write.xlsx2(pathways, file.path(directory, 'pathways.xlsx'))
}

lung_gcat_genes <- function() {
  
  chr12 <- read.csv('inst/extdata/lung/chr12_association_results.txt', sep = ' ', stringsAsFactors = FALSE)
  
  snp_dat_gcat <- chr12 %>%
    dplyr::rename(
      name=rs_id_all,
      chrom=chr,
      chromStart=position
    ) %>%
    transform(
      chromEnd=chromStart + 1
    ) %>%
    dplyr::select(
      name,
      chrom,
      chromStart,
      chromEnd,
      frequentist_add_pvalue,
      P_adjusted
    )

  gene_dat <- runseq2gene(
    inputfile = snp_dat_gcat,
    genome = "hg19",
    adjacent = TRUE,
    SNP = TRUE,
    search_radius = 150000,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )

  genes <- gene_dat$seq2gene_FullResult
  write.csv2(genes, 'output/lung/genes_gcat.csv', row.names = FALSE, quote = TRUE)
  
  gene_list <- as.data.frame(table(genes$gene_name))
  write.csv2(gene_list, 'output/lung/gene_gcat_list.csv', row.names = FALSE, quote = TRUE)
  
  genes$gene_name <- as.character(genes$gene_name)
  
  selected_variants <- genes[genes$gene_name %in% c('VWF', 'KRAS', 'KSR2', 'LUM', 'KERA'),]
  write.csv2(selected_variants, 'output/lung/selected_variants.csv', row.names = FALSE, quote = TRUE)
}

lung_gcat_genes_10_5 <- function() {
  
  chr12 <- read.csv('inst/extdata/lung/chr12_association_results.txt', sep = ' ', stringsAsFactors = FALSE) %>%
    filter(
      frequentist_add_pvalue <= 10^-5
    )
  
  snp_dat_gcat <- chr12 %>%
    dplyr::rename(
      name=rs_id_all,
      chrom=chr,
      chromStart=position
    ) %>%
    transform(
      chromEnd=chromStart + 1
    ) %>%
    dplyr::select(
      name,
      chrom,
      chromStart,
      chromEnd,
      frequentist_add_pvalue,
      P_adjusted
    )
  
  gene_dat <- runseq2gene(
    inputfile = snp_dat_gcat,
    genome = "hg19",
    adjacent = TRUE,
    SNP = TRUE,
    search_radius = 150000,
    PromoterStop = FALSE,
    NearestTwoDirection = TRUE
  )
  
  genes <- gene_dat$seq2gene_FullResult
  write.csv2(genes, 'output/lung/genes_gcat_10_5.csv', row.names = FALSE, quote = TRUE)
  
  gene_list <- as.data.frame(table(genes$gene_name))
  write.csv2(gene_list, 'output/lung/gene_gcat_list_10_5.csv', row.names = FALSE, quote = TRUE)
  
  genes$gene_name <- as.character(genes$gene_name)
  
  selected_variants <- genes[genes$gene_name %in% c('VWF', 'KRAS', 'KSR2'),]
  write.csv2(selected_variants, 'output/lung/selected_variants.csv', row.names = FALSE, quote = TRUE)
}

gwas_catalog_genes <- function() {
  
  gwas_genes <- read.csv('inst/extdata/lung/gwas-association-downloaded_2016-12-12-lung cancer.tsv', sep = '\t', stringsAsFactors = FALSE) %>%
    filter(
      CHR_ID == 12
    )
  
  gwas_genes_ids <- unlist(strsplit(gwas_genes$MAPPED_GENE, split = '\\ - ')) %>%
    strsplit(split = ',') %>%
    unlist() %>%
    unique()
  
  write.csv2(gwas_genes_ids, 'output/lung/gwas_catalog_lung_genes.csv', row.names = FALSE, col.names = FALSE)
}

clinvar_genes <- function() {
  clinvar_genes <- read.csv(file.path(directory, 'inst/extdata/lung/clinvar_result_lung_cancer.txt'), sep = '\t', stringsAsFactors = FALSE) %>%
    filter(
      Chromosome == 12
    )
  
  clinvar_genes_ids <- unlist(strsplit(clinvar_genes$Gene.s., split = '\\|')) %>%
    unique()
}

disgenet_ling <- function() {
  disgenet_lung <- read.xlsx2(file.path(directory, 'inst/extdata/lung/disgenet_lung_cancer.xls'), sheetIndex = 1)
}
  
compare <- function() {
  
  genes_gcat <- read.csv2('output/lung/gene_gcat_list.csv', stringsAsFactors = FALSE)
  genes_hf <- read.csv2('output/lung/gene_hf_list.csv', stringsAsFactors = FALSE)
  gwas_genes_ids <- read.csv2('output/lung/gwas_catalog_lung_genes.csv', stringsAsFactors = FALSE)
  
  vlung <- Venn(
    SetNames = c("GWAS", "HF", "Reported"),
    list(
      genes_gcat$Var1,
      genes_hf$Var1,
      gwas_genes_ids$x
    )
  )
  
  vlung@IntersectionSets$`111`
  
  plot(vlung, doWeights = FALSE, type = "circles")
  
}