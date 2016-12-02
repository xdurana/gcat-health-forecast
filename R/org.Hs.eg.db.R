library(org.Hs.eg.db)

function() {
  genesym <- genes.summary$Var1[1:10]
  geneid <- select(
    org.Hs.eg.db,
    keys=genesym,
    keytype="SYMBOL",
    columns="ENTREZID"
  )
  
  geneid
  
  select(
    org.Hs.eg.db,
    geneid$ENTREZID,
    c(
      "PFAM",
      "GO"
    )
  )
}