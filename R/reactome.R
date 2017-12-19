library(ReactomePA)

data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)

## [1] "4312"  "8318"  "10874" "55143" "55388" "991"
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

dotplot(x, showCategory=15)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)

cnetplot(x, categorySize="pvalue", foldChange=geneList)


x <- enrichPathway(gene = gene_list$Var1, pvalueCutoff = 0.05, readable = TRUE)
