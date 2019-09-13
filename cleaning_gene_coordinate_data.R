## Cleaning gene coordinate data

## Import gene coordinates, and remove extraneous columns
#Downloaded as ftp://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/ensGene.txt.gz
gene.ids <- read.table("~/Downloads/ensGene.txt")
gene.ids <- subset(gene.ids,V3 %in% c("chr2L","chr2R","chr3L","chr3R","chrX"))
gene.ids <- subset(gene.ids,!duplicated(V13))
gene.ids$V3 <- as.character(gene.ids$V3)
gene.ids$V3 <- as.factor(gene.ids$V3)
levels(gene.ids$V3) <- 1:5
write.table(gene.ids[,c(13,3,5,6)],"ensgenes.txt",quote=F,row.names=F,col.names=F)