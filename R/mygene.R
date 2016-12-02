library(mygene)

function() {
  genes.summary$Var1[1:10]
  res <- queryMany(genes.summary$Var1[1:10], scopes='symbol', fields=c('entrezgene', 'go'), species='human')
  res[1, 'go.BP'][[1]]
}