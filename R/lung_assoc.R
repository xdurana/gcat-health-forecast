library(data.table)
library(dplyr)

lung_bim <- do.call(rbind, lapply(c(1:22), function(chr) {
  fread(sprintf('/home/labs/dnalab/share/lims/xduran/Lung/all_genome/chr%s/chr%s_imputed.bim', chr, chr))
}))

assoc <- fread('/home/labs/dnalab/share/lims/R/health-forecast/inst/extdata/lung/assoc_results_whole_genome.txt')
assoc <- fread('/home/labs/dnalab/share/lims/R/health-forecast/inst/extdata/lung/assoc_results_whole_genome_impute_07.txt')

lung_bim <- lung_bim %>%
  rename(
    CHR=V1,
    BP=V4
  )

assoc <- filter(
  P_adjusted < 10^-4 &
  all_maf > 0.01
) %>%
  View()


assoc_bim <- merge(assoc, lung_bim) %>%
  filter(
    P_adjusted < 10^-4
  )

v2 <- unique(assoc_bim$V2)

snps <- '/home/labs/dnalab/share/lims/R/health-forecast/inst/extdata/lung/assoc_results_whole_genome_snps.txt'

write.table(
  v2,
  snps,
  row.names = FALSE,
  quote = FALSE)

lungdir <- '/home/labs/dnalab/share/lims/xduran/Lung'

lapply(
  c(1:22),
  function(chr) {
    system(
      sprintf(
        'plink --bfile %s/all_genome/chr%s/chr%s_imputed --extract %s --make-bed --out %s/output/chr%s_assoc --noweb', lungdir, chr, chr, snps, lungdir, chr))
  }
)

lapply(
  c(1:22),
  function(chr) {
    system(
      sprintf(
        'plink --bfile %s/all_genome/chr%s/chr%s_imputed --max-maf 0.1 --make-bed --out %s/output/chr%s_assoc --noweb', lungdir, chr, chr, lungdir, chr))
  }
)

system('plink --merge-list allfiles.txt --make-bed --out output/lung --noweb')

###
# RESULTS

assoc.lamplink <- fread('/imppc/labs/dnalab/xduran/fim/output/lung/assoc_0.1/assoc_ld.lamplink')

assoc.lamplink.tuples <- assoc.lamplink %>%
  transform(
    tuples = rowSums(assoc.lamplink[,10:17])
  ) %>%
  filter(
    tuples > 0
  )

### COUNT GENES AND VARIANTS

lst <- fread('/home/labs/dnalab/share/lims/R/health-forecast/output/lung/10^-4/genes_gcat.csv')

genelst <- c('ATF1', 'PAX7', 'TBX3', 'IRX5', 'IRX3', 'CERS5', 'ATF1', 'PAX7', 'PROX1', 'TBX3', 'IRX5', 'IRX3', 'CERS5',
  'PROX1', 'TBX3', 'TASP1', 'EBF2', 'ING5',
  'PTPN14', 'TBX3', 'TLE3', 'IRX5', 'SALL3', 'ESF1', 'SMYD2', 'EBF2', 'IRX3', 'ING5',
  'ADM', 'EIF4EBP2', 'PDE4D', 'RAPGEF2', 'PCLO',
  'FCGR2A', 'FCGR3B', 'FCGR2C', 'FCGR2B', 'FCGR3A',
  'AGT', 'ANGPT1', 'GATA3', 'PROX1', 'NRP1',
  'HHEX', 'PITX2', 'TLE3', 'TLE4', 'FZD4', 'PYGO1', 'WWOX', 'CXXC4', 'NKD2', 'RSPO2',
  'GATA3', 'CGA', 'GHRL', 'TBX3', 'FZD4') %>% unique()

snplst <- lst[lst$gene_name %in% genelst,]$name %>% unique()

snplst %>% length()
genelst %>% length()
