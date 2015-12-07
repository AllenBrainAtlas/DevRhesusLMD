library(OrderedList)

compareListsFast <- function(ID.List1, ID.List2, mapping = NULL, 
                             two.sided = TRUE, B = 1000, alphas = NULL, 
                             invar.q = 0.5, min.weight = 1e-05, 
                             no.reverse = FALSE) {
  res <- list()
  if (!is.null(mapping)) {
    res$mapping <- list()
    tmp <- ID.List1 %in% mapping[, 1]
    if (any(!tmp)) 
      cat(sum(!tmp), " of ", length(tmp), " elements in first list not found in mapping\n")
    res$mapping$missed1 <- sum(!tmp)
    res$mapping$rawlen1 <- length(tmp)
    tmp <- ID.List2 %in% mapping[, 2]
    if (any(!tmp)) 
      cat(sum(!tmp), " of ", length(tmp), " elements in second list not found in mapping\n")
    res$mapping$missed2 <- sum(!tmp)
    res$mapping$rawlen2 <- length(tmp)
    mapping <- mapping[mapping[, 1] %in% ID.List1, , drop = FALSE]
    mapping <- mapping[mapping[, 2] %in% ID.List2, , drop = FALSE]
    n <- nrow(mapping)
    mapIDs <- apply(mapping, 1, paste, collapse = "/")
    oo <- order(match(mapping[, 1], ID.List1))
    ID.List1 <- mapIDs[oo]
    oo <- order(match(mapping[, 2], ID.List2))
    ID.List2 <- mapIDs[oo]
  }
  tmp <- sum(!(ID.List1 %in% ID.List2))
  if (tmp > 0) 
    stop(tmp, " element(s) of first list not found in second")
  tmp <- sum(!(ID.List2 %in% ID.List1))
  if (tmp > 0) 
    stop(tmp, " element(s) of second list not found in first")
  n <- length(ID.List2)
  Ranks.List1 <- match(ID.List1, ID.List2)
  Ranks.List2 <- 1:n
  if (is.null(alphas)) {
    nn <- c(100, 150, 200, 300, 400, 500, 750, 1000, 1500, 
            2000, 2500)
    alphas <- -log(min.weight)/nn
  }
  else {
    nn <- floor(-log(min.weight)/alphas)
  }
  select <- nn < n
  alphas <- alphas[select]
  nn <- nn[select]
  nalphas <- length(alphas)
  res$n <- n
  res$call <- list(B = B, alphas = alphas, invar.q = invar.q, 
                   two.sided = two.sided, min.weight = min.weight, no.reverse = no.reverse)
  res$nn <- nn
  res$scores <- numeric(nalphas)
  res$revScores <- numeric(nalphas)
  res$pvalues <- numeric(nalphas)
  res$revPvalues <- numeric(nalphas)
  class(res) <- "listComparison"
  res$overlaps <- overlap(ID.List1, ID.List2, max(nn))
  if (no.reverse) {
    res$revOverlaps <- rep(0, length(res$overlaps))
  }
  else {
    res$revOverlaps <- overlap(ID.List1, rev(ID.List2), max(nn))
  }
  if (two.sided) {
    res$overlaps <- c(res$overlaps, rev(overlap(rev(ID.List1), 
                                                rev(ID.List2), max(nn))))
    if (no.reverse) {
      res$revOverlaps <- rep(res$revOverlaps, 2)
    }
    else {
      res$revOverlaps <- c(res$revOverlaps, rev(overlap(rev(ID.List1), 
                                                        ID.List2, max(nn))))
    }
  }
  bases <- exp(-alphas)
  res$scores <- scoreRankings(Ranks.List1, Ranks.List2, nn, 
                              bases, two.sided)
  if (no.reverse) {
    res$revScores <- rep(0, length(res$scores))
  }
  else {
    res$revScores <- scoreRankings(Ranks.List1, n + 1 - Ranks.List2, 
                                   nn, bases, two.sided)
  }
  for (i in 1:nalphas) {
    res$pvalues[i] <- 1
    if (no.reverse) {
      res$revPvalues[i] <- 1
    }
    else {
      res$revPvalues[i] <- 1
    }
  }
  res$ID.List1 <- ID.List1
  res$ID.List2 <- ID.List2
  return(res)
}
