library(biomaRt)

to_chr_position <- function(snp_ids) {
  snp_mart = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
  snp_attributes = c("refsnp_id", "chr_name", "chrom_start")
  snp_locations = getBM(
    attributes=snp_attributes,
    filters="snp_filter",
    values=snp_ids,
    mart=snp_mart
  ) %>%
    transform(
      chrom_end=chrom_start+1
    )
  snp_locations
}

mymart <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")

filters <- listFilters(mymart)
attributes <- listAttributes(mymart)

getBM(
  attributes = c("ensembl_gene_id"),
  filters = c("phenotype_description"),
  values = c("DIABETES MELLITUS INSULIN-DEPENDENT 10"),
  mart = mymart
)