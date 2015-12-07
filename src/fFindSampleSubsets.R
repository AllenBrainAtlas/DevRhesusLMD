FindSampleSubsets <- function(sample.index, sample.names) {
  sample.subset.list <- vector("list", length(sample.names))
  names(sample.subset.list) <- sample.names
  
  # For each structure, return row indices for corresponding samples
  for (sample.name in sample.names) {
    # Select subset based on name search (reg expr)
    sample.name.search <- paste0("^", sample.name, "_")
    sample.subset <- grep(sample.name.search, sample.index)
    
    # Remove GE from BG subset
    if (sample.name == "BG") {
      GE.subset <- grep("^BG_GE_", sample.index[sample.subset])
      sample.subset <- sample.subset[-GE.subset]
    }
    
    # Add cortical hem to hippocampus subset
    if (sample.name == "HP") { 
      Hem.subset <- grep("^NCX_Tcx_Hem_", sample.index)
      sample.subset <- c(sample.subset, Hem.subset)
    }
    # Save sample subset
    sample.subset.list[[sample.name]] <- sample.subset
  }
  return(sample.subset.list)
}