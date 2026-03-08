#This script performs Functional Enrichment Analysis using the fgsea package
#The results table generated from the Advanced RNA-seq analysis is used as input
#The functional data comes from the Molecular Signatures Database (as used for the microarray analysis as well)

library(fgsea)
library(dplyr)
library(ggplot2)

results <-read.csv("result_treatment.csv",header =T, row.names = 1)

#fix the rownames
head(results)

results$full_name <- rownames(results)

results$ensembl_gene_id <- sub("_.*$", "", results$full_name)

#grab the hallmark gene set from bioinf.wehi.edu.au
#generally it is better to cache download results, not download each time, but this file is small!
#You may also wish to check for updates to the annotation.
download.file("https://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata", destfile = "mouse_H_v5p2.rdata")
load("mouse_H_v5p2.rdata")

#We need to map to EntrezIDs- the annotation uses EntrezID
head(Mm.H)

#
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
#listAttributes(mart)
ens2entrez <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), mart=mart)
head(ens2entrez)

#Map gene names to entrezgene_id
#if NA remove (we only need to consider mappable genes)

results2 <- inner_join(results, ens2entrez, by=join_by("ensembl_gene_id"))

head(results2,n=40)

#sanity check one or two of the results eg
#ENSMUSG00000046101  303.147520    -0.26407455 0.1563231 -1.68928724 0.0911643983 0.277178884        240697
ens2entrez[which(ens2entrez$ensembl_gene_id  == "ENSMUSG00000046101"),]

#remove any that are NA
results2 <-results2[is.na(results2$entrezgene_id)==FALSE,]
dim(results2)

head(results2)

#This analysis takes a ranked list of numbers with entrez IDs as names
#We need to choose how to order the hits
#Here I order by stat value- we could use bidirectional ordering using FC or some combination of FC and FDR
results2 <-results2[order(results2$stat,decreasing=T),]

results2 <-results2[,c("stat","entrezgene_id")]

#remove non-unique row names- we don't mind deleting these in this case
#it is wise to check the genes that are duplicated to see if they should be excluded
results2 <-results2[!duplicated(results2$entrezgene_id),]


rownames(results2)<-results2[,"entrezgene_id"]
results2["entrezgene_id"] <-NULL


#it might want  a named vector
rn <-rownames(results2)
colnames(results2) <-NULL
results2 <-results2[,1]
names(results2)<-rn

head(results2)

#time to perform the enrichment analysis
results2  <- results2[!(is.na(names(results2)))]
results2  <- results2[!(is.na(results2))]
#results2 <-log2(results2)

fgseaRes <- fgsea(Mm.H, results2, minSize=25, maxSize = 500)

fgseaRes <-fgseaRes[order(fgseaRes$NES,decreasing=T),]
head(fgseaRes, 50)
png("FAM_enrch.png")
plotEnrichment(Mm.H[["HALLMARK_FATTY_ACID_METABOLISM"]],results2) +labs(title="HALLMARK_FATTY_ACID_METABOLISM")
dev.off()

#etc
#More analysis approaches can be found
#https://bioconductor.org/packages/release/bioc/html/fgsea.html

#plotGseatable, adapted from: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
toppathwaysup <- fgseaRes[NES>0][head(order(padj), n=10), pathway]
toppathwaysdown <- fgseaRes[NES<0][head(order(padj), n=10), pathway]
toppathways <- c(toppathwaysup, rev(toppathwaysdown))
png("top_enrch.png")
plotGseaTable(Mm.H[toppathways], results2, fgseaRes, gseaParam=0.5)+labs(title="Top pathways")
dev.off()
collapsedPathways <- collapsePathways(fgseaRes[order(padj)][padj < 0.05], 
                                      Mm.H, results2)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
png("topcolps_enrch.png")
plotGseaTable(Mm.H[mainPathways], results2, fgseaRes, 
              gseaParam = 0.5)+labs(title="Top collapsed pathways")
dev.off()

#barplot of top pathways
idx_pw <- match(toppathways, fgseaRes$pathway)
df <- data_frame(x=toppathways, y=fgseaRes$NES[idx_pw])
df <- df[order(df$y), ]
png("top_pwbar.png")
ggplot(df, aes(x=x, y=y))+
  geom_col()+
  labs(
    title = "Top pathways (up and down)",
    x="pathways",
    y="NES"
  )+
  coord_flip()
dev.off()

