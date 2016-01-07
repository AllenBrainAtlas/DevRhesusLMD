# Init workspace
if (!grepl("reports$", getwd())) setwd("reports")

# Load libraries
library(limma)
library(lmerTest)

# Load macaque starting data
load(file="../cache/nhp_PrePost_StartingData.RData")

keep.samples <- which(samplePrePost$data_set == "pre")
expr.subset <- exprPrePost[, keep.samples]
samp.subset <- droplevels(samplePrePost[keep.samples, ])

# Label genes with chromosome
ncbiInfo <- read.csv(file = "../lib/rhesus/NCBI_macaque_gene_info_strand_06Feb2014.csv",
                     row.names=2, stringsAsFactors = FALSE)
# Update Y chromosome info
chrY = scan(file = "../lib/rhesus/rhesus_ychr_genes_Hughes2012.csv", "character")
# Pseudo-autosomal genes
par.genes <- read.csv(file = "../lib/rhesus/PAR_genes.csv")
ncbiInfo$chromosome[ncbiInfo$Symbol %in% chrY] <- "Y"
chr = as.character(ncbiInfo$chromosome)
chr[is.na(chr)] <- ncbiInfo$chromosome[is.na(chr)] <- "Un"
names(chr) <- rownames(ncbiInfo)
chr <- chr[as.character(probes$macaque_entrezid)]
chr[is.na(chr)] <- "Un"
probes$chromosome <- chr
probes$chromosome[probes$human_genesymbol %in% par.genes$Approved.Symbol] <- "PAR"

# Sex DEX of all genes (not permuted)
try(lme.all <- read.csv(file = "../analysis/sex_dex/lme_sexdex.csv", 
                        row.names = 1))

if (!exists(lme.all)) {
  lme.all <- t(apply(expr.subset, 1, function(x) {
    lme1 <- lmer(x ~ sex + (1 | toplevel) + (1 | age/donor_name), 
                 data = samp.subset)
    # Report fixed effect and approx p-value
    c(fixef(lme1)["sexM"], pval = anova(lme1)$P)
  }))
  
  write.csv(lme.all, "../analysis/sex_dex/lme_sexdex.csv", row.names = FALSE)
  
  
  # Calc mean expr of genes across prenatal samples ######
  mean.expr <- apply(expr.subset, 1, mean)
}


# Permute donor sex within age (64 possible perms)
try(load(file = "../analysis/sex_dex/lme_sexdex_perm.RData"))

# VERY SLOW - RUN ONLY IF NECESSARY
if (!exists(lme.perm)) {
  tbl1 <- table(samp.subset$donor_name, samp.subset$age, samp.subset$sex)
  donor.m <- apply(tbl1[, , "M"], 2, function(x) rownames(tbl1[, , "M"])[which(x > 0)])
  donor.f <- apply(tbl1[, , "F"], 2, function(x) rownames(tbl1[, , "F"])[which(x > 0)])
  
  opt1 <- c(1, 2)
  perm <- 0
  donor.perm.all <- data.frame()
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        for (l in 1:2) {
          for (m in 1:2) {
            for (n in 1:2) {
              perm <- perm + 1
              choice1 <- c(i, j, k, l, m, n)
              ages <- levels(samp.subset$age)
              for (idx in 1:6) {
                age1 <- ages[idx]
                male1 <- cbind(perm, c(donor.m[1, age1], 
                                       donor.f[choice1[idx], age1]), "M")
                female1 <- cbind(perm, c(donor.m[2, age1], 
                                         donor.f[opt1[opt1 != choice1[idx]], age1]), "F")
                donor.perm.all <- rbind(donor.perm.all, male1, female1)
              }
            }
          }
        }
      }
    }
  }
  colnames(donor.perm.all) <- c("perm", "donor", "sex")
  
  # Calc linear mixed model with permuted sex
  lme.perm <- vector("list", 2)
  names(lme.perm) <- c("sexM", "pval")
  for (p in p.start:p.end) {
    samp.subset$sex_perm <- samp.subset$sex
    donor.perm <- subset(donor.perm.all, perm == p)
    for (i in 1:nrow(donor.perm)) {
      donor1 <- as.character(donor.perm$donor[i])
      sex1 <- as.character(donor.perm$sex[i])
      samp.subset$sex_perm[samp.subset$donor_name == donor1] <- sex1
    }
    print(table(samp.subset$sex_perm))
    
    # Run LME
    lme.genes <- t(apply(expr.subset, 1, function(x) {
      lme1 <- lmer(x ~ sex_perm + (1 | toplevel) + (1 | age/donor_name), 
                   data = samp.subset)
      # Report fixed effect and approx p-value
      c(fixef(lme1)["sex_permM"], pval = anova(lme1)$P)
    }))
    
    sex.perm <- sapply(lme.genes, function(x) x["sex_permM"])
    pval <- sapply(lme.genes, function(x) x["pval"])
    if (length(lme.perm[["sexM"]]) == 0) {
      lme.perm[["sexM"]] <- sex.perm
      lme.perm[["pval"]] <- pval
    } else {
      lme.perm[["sexM"]] <- cbind(lme.perm[["sexM"]], sex.perm)
      lme.perm[["pval"]] <- cbind(lme.perm[["pval"]], pval)
    }
  }
  
  save(lme.perm, file = "../analysis/sex_dex/lme_sexdex_perm.RData")
  
}


# Extended Data Figure 3
lme.all$pval[is.na(lme.all$pval) | lme.all$pval < 0] <- 1
lme.all$gene <- rownames(lme.all)
lme.all <- lme.all[order(lme.all$M.F_log2ratio), ]

lme.perm[["sexM"]][is.na(lme.perm[["sexM"]])] <- 0
lme.perm.p <- apply(lme.perm[["sexM"]], 2, sort)

perm <- apply(lme.perm.p[, apply(lme.perm.p, 2, sd) > 0], 1, median)
perm.qtop <- apply(lme.perm.p[, apply(lme.perm.p, 2, sd) > 0], 1, quantile, 0.975)
perm.qbottom <- apply(lme.perm.p[, apply(lme.perm.p, 2, sd) > 0], 1, quantile, 0.025)
orig <- lme.all$M.F_log2ratio

# pdf(file = "../analysis/sex_dex/lme_sexdex_qqplot.pdf", width = 5, height = 6)
plot(perm, orig, 
     col = c("red", "blue")[ifelse(lme.all$M.F_log2ratio < 0, 1, 2)], las = 1,
     ylim = c(-2.5, 3.5), 
     xlab = "Expected sex differential expression - log2(M:F ratio)", 
     ylab = "Observed sex differential expression - log2(M:F ratio)")
abline(a = 0, b = 1)
lines(perm, perm.qtop, lty = "dashed")
lines(perm, perm.qbottom, lty = "dashed")
gene.lab <- c(12440, 12441)
text(perm[gene.lab], orig[gene.lab], labels = c("EIF1AY", "LOC693361"),   
     pos = 2, cex = 0.7)
# dev.off()
