library(data.table)

genes_hf_40 <- fread('output/lung/hf_40/gene_list.csv')
genes_hf_45 <- fread('output/lung/hf_45/gene_list.csv')
genes_hf_50 <- fread('output/lung/hf_50/gene_list.csv')
genes_gwas <- fread('output/lung/gwas/gene_list.csv')

vlung <- Venn(
  SetNames = c("GWAS", "HF_45", "HF_50"),
  list(
    genes_gwas$Var1,
    genes_hf_45$Var1,
    genes_hf_50$Var1
  )
)

vlung <- Venn(
  SetNames = c("GWAS", "HF_40", "HF_50"),
  list(
    genes_gwas$Var1,
    genes_hf_40$Var1,
    genes_hf_50$Var1
  )
)

vlung@IntersectionSets$`111`
  
plot(vlung, doWeights = FALSE, type = "circles")
