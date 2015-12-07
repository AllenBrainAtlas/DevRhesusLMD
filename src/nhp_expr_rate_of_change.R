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
source(file = "src/fReorderFactorLevels.R")

ConcatenateListItems <- function(list1, from = 1, to = length(list1), 
                                 sep = "_") {
  paste(list1[from:to], collapse = sep) 
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>
# Calculate mean expr for all regions ####
# >>>>>>>>>>>>>>>>>>>>>>>>>>>

# Load data
load(file = "cache/nhp_PrePost_StartingData.RData")

# Calculate mean expr for all structures ####
agg.gene.subsets <- list(1:5000, 5001:10000, 10001:nrow(exprPrePost))
kInfoVars <- c("struc_age", "age_log2pcd")
expr.summary.bystruc.long <- data.frame()  # Init var
for (i.subset in agg.gene.subsets) {
  agg.dat <- data.frame(samplePrePost[, kInfoVars], 
                        t(exprPrePost[i.subset, ]))
  
  agg.dat.long <- melt(agg.dat, id.vars = kInfoVars)
  colnames(agg.dat.long) <- sub("variable", "gene", colnames(agg.dat.long))
  colnames(agg.dat.long) <- sub("value", "expr", colnames(agg.dat.long))
  
  # Calc mean expr by structure
  agg.dat.long.mean <- dcast(agg.dat.long, gene ~ struc_age, mean, value.var = "expr")
  expr.summary.bystruc.long <- rbind(expr.summary.bystruc.long, 
                                     melt(agg.dat.long.mean, id.vars = "gene"))
}

# Reorder data frame
colnames(expr.summary.bystruc.long)[2:3] <- c("struc_age", "expr")
expr.summary.bystruc.long <- arrange(expr.summary.bystruc.long, struc_age, gene)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>
# Calculate mean expr for cortex (layer_dev) ####
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
kInfoVars <- c("subregion", "layer_dev", "age")
keep.samples <- which(samplePrePost$toplevel == "NCX")

agg.dat <- data.frame(samplePrePost[keep.samples, kInfoVars], 
                      t(exprPrePost[, keep.samples]))

agg.dat.long <- melt(agg.dat, id.vars = kInfoVars)
colnames(agg.dat.long) <- sub("variable", "gene", colnames(agg.dat.long))
colnames(agg.dat.long) <- sub("value", "expr", colnames(agg.dat.long))

region.grouping <- "area"
if (region.grouping == "layer") {
  # Calc mean expr by layer
  agg.dat.long.mean <- dcast(agg.dat.long, gene ~ subregion + layer_dev + age,
                           mean, value.var = "expr")
} else if (region.grouping == "area") {
  # Calc mean expr by area
  agg.dat.long.mean <- dcast(agg.dat.long, gene ~ subregion + age,
                             mean, value.var = "expr")  
}

expr.summary.bystruc.long <- melt(agg.dat.long.mean, id.vars = "gene")

# Reorder data frame
colnames(expr.summary.bystruc.long)[2:3] <- c("struc_age", "expr")
expr.summary.bystruc.long <- arrange(expr.summary.bystruc.long, struc_age, gene)



# >>>>>>>>>>>>>>>>>>>>>>>>>>>
# Calculate mean expr for toplevel ####
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
region.grouping <- "toplevel"
kInfoVars <- c("toplevel", "age")

agg.dat <- data.frame(samplePrePost[, kInfoVars], t(exprPrePost))

agg.dat.long <- melt(agg.dat, id.vars = kInfoVars)
colnames(agg.dat.long) <- sub("variable", "gene", colnames(agg.dat.long))
colnames(agg.dat.long) <- sub("value", "expr", colnames(agg.dat.long))

# Calc mean expr by region
agg.dat.long.mean <- dcast(agg.dat.long, gene ~ toplevel + age,
                             mean, value.var = "expr")

expr.summary.bystruc.long <- melt(agg.dat.long.mean, id.vars = "gene")

# Reorder data frame
colnames(expr.summary.bystruc.long)[2:3] <- c("struc_age", "expr")
expr.summary.bystruc.long <- arrange(expr.summary.bystruc.long, struc_age, gene)



# >>>>>>>>>>>>>>>>>>>>>>>>>>>
# Calculate rates of expression change ####
# >>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split struc_age variable into structure and age
sample.name <- as.character(expr.summary.bystruc.long$struc_age)
sample.name.parsed <- strsplit(sample.name, split = "_")
sample.struc <- sapply(sample.name.parsed, function(x)
  ConcatenateListItems(x, from = 1, to = (length(x) - 1)))
sample.age <- sapply(sample.name.parsed, function(x)
  ConcatenateListItems(x, from = length(x), to = length(x)))
sample.age.log2pcd <- factor(sample.age, levels = levels(samplePrePost$age))
levels(sample.age.log2pcd) <- levels(as.factor(samplePrePost$age_log2pcd))
sample.age.log2pcd <- as.numeric(as.character(sample.age.log2pcd))

# Update table of mean expression by structure
expr.summary.bystruc.long$struc <- as.factor(sample.struc)
expr.summary.bystruc.long$age <- sample.age
expr.summary.bystruc.long$age_log2pcd <- sample.age.log2pcd

# Calc Event Score (cf. Translating Time, Workman et al. 2013)
event.score <- (log(2^expr.summary.bystruc.long$age_log2pcd) - 3.27) / 2.413
expr.summary.bystruc.long$event_score <- event.score

# Order samples by struc/gene/age
expr.summary.bystruc.long <- arrange(expr.summary.bystruc.long,
                                     struc, gene, age_log2pcd)

# Init vars
expr.summary.bystruc.long$expr_diff <- NA
expr.summary.bystruc.long$expr_diff_dir <- NA
expr.summary.bystruc.long$expr_diff_rate <- NA
expr.summary.bystruc.long$expr_diff_rate_eventscore <- NA

for(struc1 in levels(expr.summary.bystruc.long$struc)) {
  subset1 <- which(expr.summary.bystruc.long$struc == struc1)
  expr.summary.bystruc.long.subset <- expr.summary.bystruc.long[subset1, ]
  
  num.ages <- nlevels(factor(expr.summary.bystruc.long.subset$age))
  
  first.ages <- which(expr.summary.bystruc.long.subset$age_log2pcd == 
                        min(expr.summary.bystruc.long.subset$age_log2pcd))
  
  last.ages <- which(expr.summary.bystruc.long.subset$age_log2pcd == 
                       max(expr.summary.bystruc.long.subset$age_log2pcd))
  
  expr.summary.bystruc.long.subset$expr_diff[-first.ages] <- 
    expr.summary.bystruc.long.subset$expr[-first.ages] -
    expr.summary.bystruc.long.subset$expr[-last.ages]
  
  expr.summary.bystruc.long.subset$expr_diff_dir <- 
    ifelse(expr.summary.bystruc.long.subset$expr_diff >= 0, 
           "Doublings", "Halvings")
  
  # Calc doublings/halvings per year
  expr.summary.bystruc.long.subset$expr_diff_rate[-first.ages] <- 
    expr.summary.bystruc.long.subset$expr_diff[-first.ages] / 
    (2 ^ expr.summary.bystruc.long.subset$age_log2pcd[-first.ages] - 
       2 ^ expr.summary.bystruc.long.subset$age_log2pcd[-last.ages]) * 365
  
  # Calc doublings/halvings per unit event score (Translating Time)
  expr.summary.bystruc.long.subset$expr_diff_rate_eventscore[-first.ages] <- 
    expr.summary.bystruc.long.subset$expr_diff[-first.ages] / 
    (expr.summary.bystruc.long.subset$event_score[-first.ages] - 
       expr.summary.bystruc.long.subset$event_score[-last.ages])
  
  # Update mean expression by struc/age table
  expr.summary.bystruc.long[subset1, ] <- expr.summary.bystruc.long.subset
}

# Save expression summary
if (region.grouping == "layer") {
  save(expr.summary.bystruc.long, file = "cache/expr_dynamics/expr.summary.bystruc.long_NCXlayer_dev.RData")
} else if (region.grouping == "area") {
  save(expr.summary.bystruc.long, file = "cache/expr_dynamics/expr.summary.bystruc.long_NCXarea.RData")
} else if (region.grouping == "toplevel") {
  save(expr.summary.bystruc.long, file = "cache/expr_dynamics/expr.summary.bystruc.long_toplevel.RData")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Convert to wide format (like laminar analysis)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
load("cache/expr_dynamics/expr.summary.bystruc.long_NCXlayer_dev.RData")

# Extract cortical area/layer info
ctx.region <- gsub("_.*$", "", expr.summary.bystruc.long$struc)
ctx.layer <- gsub("^.*_", "", expr.summary.bystruc.long$struc)
expr.summary.bystruc.long$area <- ctx.region
expr.summary.bystruc.long$layer <- ctx.layer

# Reorder factor levels
age.order <- c("E40", "E50", "E70", "E80", "E90", "E120", 
               "0M", "3M", "12M", "48M")
age1 <- ReorderFactorLevels(expr.summary.bystruc.long$age, 
                            level.order = age.order, ordered = TRUE)
expr.summary.bystruc.long$age <- age1

layer.order <- rev(c("Hem", "WM", "VZi", "VZo", "VZ", "SZi", "SZo", "SZ", 
                     "IZ", "IFZ", "TMZ", "OFZ", "ICD", "SP", "CPi", "L6", "L5", 
                     "CPo", "CP", "L4", "L4Cb", "L4Ca", "L4B", "L4A", 
                     "L3", "L2-3", "L2.3", "L2", "L1", "MZ"))
layer1 <- ReorderFactorLevels(expr.summary.bystruc.long$layer, 
                              level.order = layer.order, ordered = FALSE)
expr.summary.bystruc.long$layer <- layer1

# Keep ACG, S1, V1 with expr rate info
keep.layers <- c("VZ", "SZ", "IZ", "SP", "L6", "L5", "L4", "L3", "L2", "MZ")
keep.struc <- with(expr.summary.bystruc.long, 
                   which(area %in% c("ACG", "S1", "V1") & 
                           layer %in% keep.layers & 
                           !is.na(expr_diff_rate)))
expr.summary.bystruc.long <- droplevels(expr.summary.bystruc.long[keep.struc, ])

# Convert to wide format
for (region1 in c("ACG", "V1")) {
  ages <- levels(factor(subset(expr.summary.bystruc.long, area == region1)$age))
  expr.rate.all <- vector("list", length(ages))
  names(expr.rate.all) <- ages
  for (age1 in ages) {
    expr.summary.subset <- droplevels(subset(expr.summary.bystruc.long, 
                                             area == region1 & age == age1))
    expr.summary.subsetw <- dcast(expr.summary.subset, gene ~ layer, 
                                  value.var="expr_diff_rate")
    # Reorder genes to match probes table
    gene.order <- match(probes$macaque_genesymbol, expr.summary.subsetw$gene)
    expr.rate.all[[age1]] <- expr.summary.subsetw[gene.order, ]
  }
  # Save rates of change by cortical area
  save(expr.rate.all, file=paste0("cache/expr_dynamics/expr.rate.all_", 
                                  region1, ".RData"))
}
