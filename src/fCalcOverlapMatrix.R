# Load libraries
library(OrderedList)
source(file = "../src/fcompareListsFast.R")

CalcOverlapMatrix <- function (max_genes, alphas=0.01151) {
  gene_list_comparison_summary <- vector("list", 2)
  names(gene_list_comparison_summary) <- c("score", "overlap")
  for (m1 in names(gene_list_comparison_summary)) {
    matrix1 <- matrix(0, nrow=ncol(max_genes), ncol=ncol(max_genes), 
                      dimnames=list(colnames(max_genes), colnames(max_genes)))
    gene_list_comparison_summary[[m1]] <- matrix1
  }
  
  for (layer1 in 1:ncol(max_genes)) {
    for (layer2 in 1:ncol(max_genes)) {
      # Compare only top n ranked genes (n=100 ~ alpha=0.1151)
      gene_list_comparison <- compareListsFast(max_genes[,layer1], 
                                               max_genes[,layer2],
                                               alphas=alphas, two.sided=FALSE, 
                                               no.reverse=TRUE)
      gene_list_overlap <- getOverlap(gene_list_comparison, percent=1)
      score1 <- gene_list_overlap$score
      overlap1 <-  gene_list_overlap$direction * 
        length(gene_list_overlap$intersect)
      gene_list_comparison_summary[["score"]][layer1,layer2] <- score1
      gene_list_comparison_summary[["overlap"]][layer1,layer2] <- overlap1
    }
  }
  return(gene_list_comparison_summary)
}
