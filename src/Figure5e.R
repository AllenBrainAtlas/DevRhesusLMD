# Init workspace
if (!grepl("reports$", getwd())) setwd("reports")

# Load libraries
library(pheatmap)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
source("../src/fReorderFactorLevels.R")
source("../src/mdmr_r.beta.1.0.1.r")


# Load macaque starting data
load(file="../cache/nhp_PrePost_StartingData.RData")
keep.samples <- which(
  samplePrePost$subregion %in% c("ACG", "V1") 
  # samplePrePost$toplevel == "NCX" 
  # samplePrePost$data_set == "pre"
)
# keep.samples <- which(samplePrePost$data_set == "pre")
expr.subset <- exprPrePost[, keep.samples]

gibbs <- read.csv("../lib/Karaca2015/gibbs_genes.csv")
gibbs <- gibbs[gibbs$Source == "Karaca2015", -grep("Source", colnames(gibbs))]
mgene.id <- match(gibbs$Gene, probes$macaque_genesymbol)
hgene.id <- match(gibbs$Gene, probes$human_genesymbol)
hgene.id[is.na(hgene.id)] <- mgene.id[is.na(hgene.id)]
mgene.sym <- probes$macaque_genesymbol[hgene.id]
hgene.sym <- mgene.sym
for (i in 1:length(mgene.sym)) {
  hum.sym <- probes$human_genesymbol[probes$macaque_genesymbol == mgene.sym[i]]
  if (! is.na(hum.sym[1])) {
    hgene.sym[i] <- hum.sym
  }
}

gibbs$Gene2 <- hgene.sym
gibbs <- na.omit(gibbs)
gibbs.annot <- apply(gibbs[, ! grepl("Gene", colnames(gibbs))], 2, 
                     function(x) tapply(x, gibbs$Gene2, max))
gibbs.annot <- as.data.frame(gibbs.annot[, apply(gibbs.annot, 2, sd) > 0])
gibbs.annot$CC_abnormal <- apply(gibbs.annot[, c("CC_agenesis", "CC_dysgenesis", 
                                                 "CC_hypoplasia")], 1, 
                                 function(x) ifelse(any(x == 1), 1, 0))

keep.genes <- mgene.sym[match(rownames(gibbs.annot), hgene.sym)]
cor1 <- cor(t(expr.subset[keep.genes, ]))
colnames(cor1) <- rownames(cor1) <- rownames(gibbs.annot)
gibbs.annot2 <- gibbs.annot[, c(4, 12, 7, 5, 13, 8, 3, 2, 1)]
colnames(gibbs.annot2) <- rev(c("Candidate gene", "Primary microcephaly", 
                                "Microcephaly", "CC agenesis", "CC abnormality", 
                                "Cortical atrophy", "Cerebral atrophy", 
                                "Brain malformations", "Seizures"))

# Heatmap
micro.pal <- brewer.pal(9, "Reds")
atrophy.pal <- brewer.pal(9, "Oranges")
cc.pal <- brewer.pal(9, "Blues")
mri.pal <- brewer.pal(9, "Greens")
seiz.pal <- brewer.pal(9, "Purples")
annot.colors = list(`Candidate gene` = c("grey95", "black"),
                    `Primary microcephaly` = c(micro.pal[1], micro.pal[7]),
                    `Microcephaly` = c(micro.pal[1], micro.pal[3]),
                    `Cortical atrophy` = c(atrophy.pal[1], atrophy.pal[6]),
                    `Cerebral atrophy` = c(atrophy.pal[1], atrophy.pal[3]),
                    `CC agenesis` = c(cc.pal[1], cc.pal[6]),
                    `CC abnormality` = c(cc.pal[1], cc.pal[4]),
                    `Brain malformations` = c(mri.pal[1], mri.pal[5]),
                    `Seizures` = c(seiz.pal[1], seiz.pal[5]))

pheatmap(cor1, border_color = NA,
         annotation_row = gibbs.annot2, annotation_legend = FALSE, 
         annotation_colors = annot.colors, show_colnames = FALSE,
         treeheight_row = 0, cutree_rows = 4, cutree_cols = 4)

# Change orientation
pheatmap(cor1, border_color = NA,
         annotation_row = gibbs.annot2, annotation_legend = FALSE, 
         annotation_colors = annot.colors, 
         # show_colnames = FALSE,
         treeheight_col = 0, cutree_rows = 4, cutree_cols = 4)


# Phenotype var explained (MDMR) ########
dist1 <- sqrt(2 * (1 - cor1))
kPerm <- 1000
mdmr1 <- MDMR(DMATRIX=dist1, DATA=gibbs.annot, PERMS=kPerm, UNIVARIATE=TRUE) 
mdmr1 <- mdmr1$UNIVARIATE
mdmr1 <- mdmr1[order(-mdmr1[,"PVE"]),]
mdmr1
