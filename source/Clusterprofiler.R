################################
#Project: UVB_UA RNA-seq project
#Package: Clusterprofiler
#Author: Ran Yin
#Date: 12-6-2018
################################

#Require org.Hs.eg.db
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")

library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
library(readr)

colnames(wk2_UA)
#input wk2 data
wk2_UA <- read_csv("data/wk2 UA.csv")

#selec up or down regulated genes from wk2
wk2_up <- unique(subset(wk2_UA,
                        log2FoldChange >=2))

wk2_down <- unique(subset(wk2_UA,
                        log2FoldChange <=-2))

#Write CSV file for clusterprofiler
write.csv(wk2_up, file = "./data/For Clusterprofiler/wk2up.csv")
write.csv(wk2_down, file = "./data/For Clusterprofiler/wk2down.csv")







x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)

data(gcSample)
hg <- gcSample[[1]]
head(hg)

eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)

bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)


