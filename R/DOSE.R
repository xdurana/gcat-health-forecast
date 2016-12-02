library(DOSE)

a <- c("DOID:14095", "DOID:5844", "DOID:2044", "DOID:8432", "DOID:9146",
       "DOID:10588", "DOID:3209", "DOID:848", "DOID:3341", "DOID:252")
b <- c("DOID:9409", "DOID:2491", "DOID:4467", "DOID:3498", "DOID:11256")
doSim(a[1], b[1], measure="Wang")
s <- doSim(a, b, measure="Wang")
s

simplot(s,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=5,
        font.size=14, xlab="", ylab="")


g1 <- c("84842", "2524", "10590", "3070", "91746")
g2 <- c("84289", "6045", "56999", "9869")

geneSim(g1[1], g2[1], measure="Wang", combine="BMA")
gs <- geneSim(g1, g2, measure="Wang", combine="BMA")
gs

simplot(gs,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=5,
        font.size=14, xlab="", ylab="")

g3 <- c("57491", "6296", "51438", "5504", "27319", "1643")
clusters <- list(a=g1, b=g2, c=g3)
mclusterSim(clusters, measure="Wang", combine="BMA")


library(DOSE)
data(geneList)
y <- gseDO(geneList,
           nPerm         = 100, 
           minGSSize     = 120,
           pvalueCutoff  = 0.2, 
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(y, 3)

ncg <- gseNCG(geneList,
              nPerm         = 100, 
              minGSSize     = 120,
              pvalueCutoff  = 0.2, 
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3)

cnetplot(ncg, categorySize="pvalue", foldChange=geneList)

enrichMap(y, n=20)
