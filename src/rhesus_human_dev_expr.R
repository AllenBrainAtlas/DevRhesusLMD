# Load libraries 
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(boot)

# Load functions 
source(file="../src/fReorderFactorLevels.R")
source(file="../src/fConvertPcdtoEventScore.R")

##### Load data #####
# Load macaque starting data
load(file="../cache/nhp_PrePost_StartingData.RData")

# Load human starting data
load(file="../lib/human/brainspan_longitudinal.RData", verbose = TRUE)

# Load human starting data (BrainCloud dlPFC)
load(file="../lib/human/human_BrainCloud.RData", verbose = TRUE)


#######
# Macaque
#######
# Keep subset of macaque samples
m.region <- "ACG"  # ACG / V1
keep.samples <- which(samplePrePost$subregion == m.region)
# keep.samples <- which(samplePrePost$region == m.region)
sample.subset <- droplevels(samplePrePost[keep.samples, ])
expr.subset <- exprPrePost[, keep.samples]

age.logpcd <- factor(sample.subset$age_log2pcd)

# Calculate mean expression by age (macaque)
expr.mean <- t(apply(expr.subset, 1, function(x)
  tapply(x, age.logpcd, mean)))
age.levels <- as.numeric(levels(age.logpcd))

# Relabel macaque genes with human orthologs
rownames(expr.mean) <- probes$human_genesymbol

# Calc most variable genes (macaque)
expr.var <- apply(expr.mean, 1, sd)
expr.range <- apply(expr.mean, 1, function(x) max(x) - min(x))

########
# Human - BrainSpan
########
# Calculate mean expression by age (human)
h.region <- "MFC"  # MFC / V1C
h.age <- numAge2[age2]
h.keep.samples <- which(reg2 == h.region)
# h.keep.samples <- which(reg2 == h.region & h.age <= 6000) # Remove later ages
h.age.logpcd <- factor(round(log2(h.age[h.keep.samples]), 2))
dat2.subset <- dat2[, h.keep.samples]
h.expr.mean <- t(apply(dat2.subset, 1, function(x) tapply(x, h.age.logpcd, mean)))
# Remove genes without variation
h.sd <- apply(h.expr.mean, 1, sd)
h.expr.mean <- h.expr.mean[h.sd > 0, ]
h.age.levels <- as.numeric(levels(h.age.logpcd))

########
# Human - BrainCloud
########
# Calc mean expression by binned age
h2.region <- "DFC"
h2.keep.samples <- which(samp.dat$age_yrs < 41)
samp.dat.subset <- droplevels(samp.dat[h2.keep.samples, ])
agepcd.cat <- round(log2(tapply(samp.dat.subset$age_pcd, 
                                samp.dat.subset$age_cat, mean)), 2)
h2.age.logpcd <- factor(agepcd.cat[samp.dat.subset$age_cat])
expr.dat3.subset <- expr.dat3[, h2.keep.samples]
h2.expr.mean <- t(apply(expr.dat3.subset, 1, 
                        function(x) tapply(x, h2.age.logpcd, mean)))
h2.age.levels <- as.numeric(levels(h2.age.logpcd))


# Get number of genes
num.genes <- nrow(expr.mean)
h.num.genes <- nrow(h.expr.mean)
h2.num.genes <- nrow(h2.expr.mean)


# Keep orthologs and match gene order
ortho.genes <- Reduce(intersect, list(rownames(expr.mean), 
                                      rownames(h.expr.mean), 
                                      rownames(h2.expr.mean)))
ortho.genes <- sort(ortho.genes)
keep.genes <- match(ortho.genes, rownames(expr.mean))
expr.mean.subset <- expr.mean[keep.genes, ]

h.keep.genes <- match(ortho.genes, rownames(h.expr.mean))
h.expr.mean.subset <- h.expr.mean[h.keep.genes, ]

h2.keep.genes <- match(ortho.genes, rownames(h2.expr.mean))
h2.expr.mean.subset <- h2.expr.mean[h2.keep.genes, ]

# Calc expr z-scores
expr.mean.subsetz <- t(apply(expr.mean.subset, 1, scale))
colnames(expr.mean.subsetz) <- colnames(expr.mean.subset)

h.expr.mean.subsetz <- t(apply(h.expr.mean.subset, 1, scale))
colnames(h.expr.mean.subsetz) <- colnames(h.expr.mean.subset)

h2.expr.mean.subsetz <- t(apply(h2.expr.mean.subset, 1, scale))
colnames(h2.expr.mean.subsetz) <- colnames(h2.expr.mean.subset)

# Store data in lists
species.expr <- list(expr.mean.subset, 
                     h.expr.mean.subset, 
                     h2.expr.mean.subset)
names(species.expr) <- c("macaque", "human", "human_bc")

species.exprz <- list(expr.mean.subsetz, 
                      h.expr.mean.subsetz, 
                      h2.expr.mean.subset)
names(species.exprz) <- c("macaque", "human", "human_bc")


##### Fit dev expr trend #####
# Combine species data
dev.expr <- rbind(data.frame(species="macaque", melt(expr.mean.subset), 
                             exprz=melt(expr.mean.subsetz)$value), 
                  data.frame(species="human", melt(h.expr.mean.subset), 
                             exprz=melt(h.expr.mean.subsetz)$value), 
                  data.frame(species="human_bc", melt(h2.expr.mean.subset), 
                             exprz=melt(h2.expr.mean.subsetz)$value))
colnames(dev.expr) <- c("species", "gene", "log2pcd", "expr", "exprz")

# Calc Translating Time escore
dev.expr$escore <- apply(dev.expr, 1, 
                         function(x) ConvertPcdtoEventScore(as.numeric(x["log2pcd"]), 
                                                            x["species"]))
dev.expr$escore <- round(dev.expr$escore, 2)

# Calc spline fit
# Select escores that are near all 3 species
escore1 <- c(0.27, 0.36, 0.46, 0.54, 0.75, 0.85, 0.96, 1.32)  # All species
# escore1 <- c(0.25, 0.38, 0.45, 0.49, 0.58, 0.76, 0.91, 1.28)  # Human vs. rhesus
dev.pred <- list()
dev.pred.se <- list()
for (species1 in c("macaque", "human", "human_bc")) {
  species.dev <- subset(dev.expr, species==species1)
  pred.df <- NULL
  pred.se.df <- NULL
  for (i in 1:length(ortho.genes)) {
    # Select rows corresponding to each gene
    gene.row <- seq(i, nrow(species.dev), by=length(ortho.genes))
    gene.dev <- species.dev[gene.row, ]
    
    # Fit expr trend
    gene.lo <- loess(expr ~ escore, gene.dev, degree=1)
    expr.pred <- predict(gene.lo, data.frame(escore = escore1), se = TRUE)
    pred.df <- rbind(pred.df, expr.pred$fit)
    pred.se.df <- rbind(pred.se.df, expr.pred$se.fit)
  }
  rownames(pred.df) <- ortho.genes
  rownames(pred.se.df) <- ortho.genes
  dev.pred[[species1]] <- pred.df
  dev.pred.se[[species1]] <- pred.se.df
}
save(dev.expr, dev.pred, dev.pred.se, 
     file = "../cache/dev_expr_species/dev.pred_hrh.RData")



# Calc dev expr correlations
cor.pairs <- list(c("human", "human_bc"), c("human", "macaque"), 
                  c("human_bc", "macaque"))
cor.all <- matrix(NA, length(ortho.genes), length(cor.pairs), 
                  dimnames = list(ortho.genes, c("h2h", "hrh", "h2rh")))
for (i in 1:ncol(cor.all)) {
  pair1 <- dev.pred[[cor.pairs[[i]][1]]]
  pair2 <- dev.pred[[cor.pairs[[i]][2]]]
  cor1 <- sapply(ortho.genes, function(x) cor(pair1[x, ], pair2[x, ], 
                                              use = "pairwise"))
  cor.all[, i] <- cor1
}
save(cor.all, file = "../cache/dev_expr_species/cor.all_hrh.RData")
