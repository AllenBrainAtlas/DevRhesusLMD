# >>>>>>>>>>>>>>>>>>>>>>
# * Initialize workspace ####
# >>>>>>>>>>>>>>>>>>>>>>

# Load libraries
library("reshape2")
library("dplyr")
library("relaimpo")


# >>>>>>>>>>>>>>>>>>>>>>>>
# * Define functions ####
# >>>>>>>>>>>>>>>>>>>>>>>>

source(file="src/fFindSampleSubsets.R")
source(file="src/fReorderFactorLevels.R")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# * Load data files #######
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Load starting data: expr, samples, and probes
load(file="cache/nhp_PrePost_StartingData.RData")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# * Calc laminar spatial enrichment (V1 L1-WM, L4 combined) #######
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Combine MZ/L1 and IZ/WM
samplePrePost$layer_dev <- gsub("MZ", "L1", samplePrePost$layer_dev)
samplePrePost$layer_dev <- gsub("IZ", "WM", samplePrePost$layer_dev)

# Select groups of samples (V1 L1-WM)
struc.subset <- with(samplePrePost, 
                     which(subregion == "V1" &
                             layer_dev %in% c("WM", paste0("L", 1:6)) & 
                             ! age %in% paste0("E", c(40, 50, 70, 80))))

# Perform independent laminar analysis of each sample group
exprPrePost.subset <- exprPrePost[, struc.subset]
samplePrePost.subset <- droplevels(samplePrePost[struc.subset, ])

laminar.enrichment.all <- vector("list", nlevels(samplePrePost.subset$age))
names(laminar.enrichment.all) <- levels(samplePrePost.subset$age)

for (timepoint in levels(samplePrePost.subset$age)) {
  age.subset <- which(samplePrePost.subset$age == timepoint)
  coef1 <- data.frame()
  for (gene in probes$macaque_genesymbol) {
    gene.subset <- which(probes$macaque_genesymbol == gene)
    expr1 <- exprPrePost.subset[gene.subset, age.subset]
    expr1.scaled <- scale(expr1)
    
    layer1 <- factor(samplePrePost.subset$layer_dev[age.subset])
    layer1.levels <- levels(layer1)
    layer1.num <- length(layer1.levels)
    
    # Create effects contrast formula (grand mean deviation coding)
    lm1.vi.all <- matrix(NA, layer1.num, layer1.num)
    
    # Due to multi-collinearity of templates, leave each template out
    for (contr.layer in 1:layer1.num) {
      p1 <- paste0("layer_dev == '", layer1.levels[-contr.layer], "'")
      p2 <- paste0("layer_dev == '", layer1.levels[contr.layer], "'")
      p3 <- paste0("scale(I((", p1, ") - (", p2, ")))")
      p4 <- paste(p3, collapse = " + ")
      
      contr.layer.formula <- as.formula(paste0("expr1.scaled ~ ", p4))
      lm1 <- lm(contr.layer.formula, 
                data = samplePrePost.subset[age.subset, ])
      
      # Estimate proportion of variance explained by each structure
      # lmg method preferred but slower than car transformation (cf Bi 2012)
      lm1.vi <- calc.relimp(lm1, type = "car", rela = TRUE)
      lm1.R2 <- lm1.vi$R2
      lm1.vi.sign <- lm1.vi$car * sign(summary(lm1)$coef[-c(1), "Estimate"])
      lm1.vi.all[-contr.layer, contr.layer] <- lm1.vi.sign
    }
    # Find average variable importance across sets of templates
    lm1.vi.mean <- apply(lm1.vi.all, 1, mean, na.rm=TRUE)
    lm1.vi.r2 <- lm1.vi.mean / sum(abs(lm1.vi.mean)) * lm1.R2
    
    # Clean up layer names
    names(lm1.vi.r2) <- layer1.levels
    coef1 <- rbind(coef1, data.frame(timepoint, gene, lm1.R2, t(lm1.vi.r2)))
  }
  laminar.enrichment.all[[timepoint]] <- coef1
}

# Update rownames with gene symbols
for (timepoint in levels(samplePrePost.subset$age)) { 
  rownames(laminar.enrichment.all[[timepoint]]) <- 
    laminar.enrichment.all[[timepoint]]$gene
}

# Save laminar enrichment data ####
laminar.folder <- "cache/laminar_genes/"
laminar.file <- paste0("laminar.enrichment.all_", struc.name, "_L1-WM_dev.RData")
save(laminar.enrichment.all, file = paste0(laminar.folder, laminar.file))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# * Find V1 laminar genes by age/layer (L1-WM only, L4 combined) ####
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# List of structures for analysis
struc.names <- "NCX_Ocx_V1"
laminar.folder <- "cache/laminar_genes/"
laminar.summary.all <- data.frame()

for (struc.name in struc.names) {
  # Load laminar enrichment information for structure
  laminar.file <- paste0("laminar.enrichment.all_", struc.name, "_L1-WM_dev.RData")
  load(file = paste0(laminar.folder, laminar.file))
  
  # Find laminar genes for each age and layer
  for (age1 in names(laminar.enrichment.all)) {
    laminar.enrichment <- laminar.enrichment.all[[age1]]
    layer.col <- ! colnames(laminar.enrichment) %in% 
      c("timepoint", "gene", "lm1.R2", "toplevel")
    layer.R2 <- laminar.enrichment[, layer.col]
    
    # For each gene, find partial R^2 of 2nd most enriched layer
    penult.layer.R2 <- apply(layer.R2, 1, function(x) x[order(-x)[2]])
    
    # Find laminar enriched and specific genes
    enriched.threshold <- 0.25
    specific.threshold <- 0.1
    for (layer1 in colnames(layer.R2)) {
      gene.index <- which(layer.R2[, layer1] > enriched.threshold)
      laminar.enriched.genes <- laminar.enrichment$gene[gene.index]
      laminar.enriched <- gene.index > 0
      
      # Gene is specific to layer if no other layer explains more than
      # threshold of variance
      laminar.specific <- abs(penult.layer.R2[gene.index]) < specific.threshold
      
      # Keep same number of columns if no laminar genes
      if (length(gene.index) == 0) {
        laminar.enriched.genes <- NA
        laminar.enriched <- NA
        laminar.specific <- NA
      }
      
      # Create table summary of laminar genes
      laminar.summary <- data.frame(region = struc.name, age1, layer1, 
                                    gene = laminar.enriched.genes,
                                    enriched = laminar.enriched, 
                                    specific = laminar.specific)
      laminar.summary.all <- rbind(laminar.summary.all, laminar.summary)
    }
  }
}  # Repeat laminar summary analysis for next structure

# Reorder age/layer factors
laminar.summary.all$age1 <- 
  ReorderFactorLevels(laminar.summary.all$age1, 
                      level.order = levels(samplePrePost$age), 
                      ordered = TRUE)
laminar.summary.all$layer1 <- 
  ReorderFactorLevels(laminar.summary.all$layer1, 
                      level.order = levels(samplePrePost$layer), 
                      ordered = FALSE)

# Rename columns
colnames(laminar.summary.all) <- gsub("1", "", colnames(laminar.summary.all))

# Add macaque Entrez ids
match.genes <- match(laminar.summary.all$gene, probes$macaque_genesymbol)
laminar.summary.all$geneid <- probes$macaque_entrezid[match.genes]

# Save table of laminar genes
write.csv(laminar.summary.all, file = "analysis/laminar_genes/laminar.summary.all_L1-WM_dev.csv", row.names = FALSE)


# Get lists of genes that have maximal laminar enrichment
max.laminar.genes <- NULL
for (timepoint in names(laminar.enrichment.all)) {
  gene.list <- NULL
  layer.names <- colnames(laminar.enrichment.all[[timepoint]])
  layer.names <- layer.names[! layer.names %in% c("timepoint", "gene", "lm1.R2")]
  for (layer in layer.names) {
    gene.list <- cbind(gene.list, laminar.enrichment.all[[timepoint]][, "gene"][order(laminar.enrichment.all[[timepoint]][, layer], decreasing = TRUE)])
  }
  colnames(gene.list) <- paste(timepoint, layer.names, sep = "_")
  max.laminar.genes <- cbind(max.laminar.genes, gene.list)
}

write.table(max.laminar.genes, file = "analysis/laminar_genes/V1(L1-WM)_max_laminar_genes_byage.txt", quote = FALSE, row.names = FALSE, sep = "\t")


