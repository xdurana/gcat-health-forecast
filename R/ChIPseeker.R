library(ChIPseeker)

function() {
  compKEGG <- compareCluster(
    geneCluster = genes,
    fun = "enrichKEGG",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  plot(
    compKEGG,
    showCategory = 15,
    title = "KEGG Pathway Enrichment Analysis"
  )
}