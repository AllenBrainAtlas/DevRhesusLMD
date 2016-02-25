# Load functions
library(limma)
library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(RColorBrewer)


#### Extended Data Figure 7a ####

# Load data (human, human2, rhesus)
load(file = "../cache/dev_expr_species/dev.pred_hrh.RData")
load(file = "../cache/dev_expr_species/cor.all_hrh.RData")

# Plot correlation as a function of gene variability
sd.mean <- apply(dev.pred[["macaque"]], 1, sd, na.rm = TRUE)
sd.sets <- seq(1, length(sd.mean), by = 100)
sd.range <- sort(sd.mean)[sd.sets]

# cons.compare <- colnames(cor.boot)[seq(2, 13, by = 3)]
cons.compare <- colnames(cor.all)[-3]
cons.prop <- matrix(NA, length(sd.range), length(cons.compare), 
                    dimnames = list(sd.range, cons.compare))
for (cons1 in cons.compare) {
  cor1 <- cor.all[, cons1]
  #   prop1 <- sapply(sd.range, function(x) sum(cor1[sd.mean >= x] > 0.5) / 
  #                     sum(sd.mean >= x))
  prop1 <- sapply(sd.range, function(x) median(cor1[sd.mean >= x]))
  cons.prop[, cons1] <- prop1
}

top.n.genes <- length(sd.mean) - sd.sets + 1
plot(top.n.genes, cons.prop[, 1], ylim = c(0.8, 1), type = "n", las = 1,
     xlab = "Top N most variable genes", 
     ylab = "Median correlation")
for (i in 1:ncol(cons.prop)) {
  lines(top.n.genes, cons.prop[, i], col = i, lwd = 4)
}
legend("bottomleft", fill = 1:2, cex = 0.7, bty = "n", 
       legend = c("Human vs. Human 2", "Human vs. Rhesus"))



#### Extended Data Figure 7b ####

# source(file = "../src/rhesus_human_dev_expr.R")  # Run to generate files loaded below
load(file = "../cache/dev_expr_species/dev.expr_dev.pred_ACG.RData")
load(file = "../cache/dev_expr_species/cor.all_hrhrat.RData")

# Define variable genes
sd1 <- apply(dev.pred[["macaque"]], 1, sd)
var.genes <- which(sd1 > quantile(sd1, 0.5))

# Plot conserved genes vs. correlation threshold
cons.trend <- NULL
q.list <- c(seq(0, 0.9, by = 0.1), 0.99)
par(mfrow = c(1, 4))
for (q1 in q.list) {
  var.genes <- which(sd1 > quantile(sd1, q1))
  cor.subset <- cor.all[var.genes, ]

  thresh.list <- seq(0, 0.9, by = 0.1)
  venn.df <- NULL
  for (cor.thresh1 in thresh.list) {
    venn.matrix <- cor.subset
    venn.matrix[1:length(venn.matrix)] <- 0
    venn.matrix[cor.subset > cor.thresh1] <- 1
    venn1 <- vennCounts(venn.matrix)[, "Counts"]
    venn.cnt <- c(venn1[c(3, 2, 5)], sum(venn1[c(4, 6:8)]))
    venn.prop <- venn.cnt / sum(venn.cnt)
    venn.df <- rbind(venn.df, venn.prop)
  }
  colnames(venn.df) <- c("Human-specific", "Rhesus-specific", 
                         "Primate-specific", "Conserved")
  rownames(venn.df) <- thresh.list
  barplot(t(venn.df), border = NA, col = brewer.pal(ncol(venn.df), "Set1"), 
          xlab = "Correlation threshold", ylab = "Proportion of genes", 
          main = paste("Top", length(var.genes), "most variable genes"), las = 1)
  if (q1 == 0.99) {
    legend("topright", fill = rev(brewer.pal(ncol(venn.df), "Set1")), 
           legend = rev(colnames(venn.df)))
  }
  
  # Store values for cor threshold of 0.5
  cons.trend <- rbind(cons.trend, cbind(length(var.genes), venn.df))
}
colnames(cons.trend)[1] <- "num_genes"

cons.trendl <- melt(as.data.frame(cons.trend), id = "num_genes")

ggplot(cons.trendl, aes(x = num_genes, y = value, col = variable)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line", size = 1) +
  xlab("Top N most variable genes") +
  ylab("Proportion of correlated genes") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank()) +
  guides(color = guide_legend(reverse = TRUE))

# Mean +- SD proportion of genes in each evolutionary category
with(subset(cons.trendl, num_genes != 42), tapply(value, variable, mean))
with(subset(cons.trendl, num_genes != 42), tapply(value, variable, sd))
