library(topGO)

function() {
  geneList <- genes.summary$Var1[1:10]
  sampleGOdata <- new(
    "topGOdata",
    description = "Simple session", ontology = "BP",
    allGenes = geneList, geneSel = topDiffGenes,
    nodeSize = 10,
    annot = annFUN.db, affyLib = affyLib
  )
}