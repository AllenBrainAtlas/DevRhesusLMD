# Load libraries
library(GOstats)
library(org.Hs.eg.db)

# Load macaque starting data
if (! exists("probes")) {
  load(file="../cache/nhp_PrePost_StartingData.RData")
}

# Define background set of genes
entrezUniverse <- unique(na.omit(probes$human_entrezid))

# Define function
CalcGoEnrich <- function(gene.ids, id.type="entrez_id", 
                         hgCutoff = 0.05, onto="BP", 
                          cond = FALSE) {
  # Check if need to lookup Entrez IDs
  if (id.type != "entrez_id") {
    gene.ids <- probes$human_entrezid[match(gene.ids, 
                                            probes$macaque_genesymbol)]
    gene.ids <- na.omit(unique(gene.ids))
  }
  
  # Check if enough gene ids entered
  if (length(gene.ids) < 10 | all(is.na(gene.ids))) {
    return(NA)
  } else {
    params <- new("GOHyperGParams", 
                  geneIds = gene.ids, 
                  universeGeneIds = entrezUniverse, 
                  annotation = "org.Hs.eg.db", # Human entrez id to GO mapping
                  ontology = onto,   # BP, CC or MF
                  pvalueCutoff = hgCutoff, 
                  conditional = cond, 
                  testDirection = "over")
    go1 <- hyperGTest(params)
    go1.p <- pvalues(go1)
    go1.p <- go1.p[go1.p < hgCutoff]  # Return GO ids with significant pvalues
    return(go1.p)
    # sigCategories(go1)
    # go1.df <- summary(go1)
    # return(go1.df)
  }
}
