ReorderFactorLevels <- function(var1, level.order, ordered = FALSE) {
  # Reorder all or a subset of factor levels of a vector
  # Args:
  #   var1: charactor or factor vector with levels to order
  #   level.order: ordering for a subset of existing factor levels
  #   ordered: create ordered factor type?
  #
  # Returns:
  #   factor vector with reordered levels
  if (!"factor" %in% class(var1)) { var1 <- as.factor(var1) }
  current.levels <- levels(var1)
  level.order <- level.order[level.order %in% current.levels]
  missing.levels <- current.levels[! current.levels %in% level.order]
  new.levels <- c(level.order, missing.levels)
  var1.reordered <- factor(var1, ordered = ordered, levels = new.levels)
  return(var1.reordered)
}
