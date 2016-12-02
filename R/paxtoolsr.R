library(paxtoolsr)
library(igraph)

function() {
  browseVignettes("paxtoolsr")
  
  gene <- genes.summary$Var1[1:10]
  t1 <- graphPc(source = gene, kind = "neighborhood", format = "BINARY_SIF", verbose = TRUE)
  t2 <- t1[which(t1[, 2] == "controls-state-change-of"), ]
  g <- graph.edgelist(as.matrix(t2[, c(1, 3)]), directed = FALSE)
  plot(g, layout = layout.fruchterman.reingold)
}