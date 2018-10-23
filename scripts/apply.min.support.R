apply.min.support <- function (phy, tol = 0.5) 
  # Adapted di2multi function from the ape package to plot polytomies
  # based on numeric node support values
  # (di2multi does this based on edge lengths)
  # Needs adjustment for unrooted trees as currently skips the first edge
{
    phy$node.label[phy$node.label==""] = "100"
  if (is.null(phy$edge.length)) 
    stop("the tree has no branch length")
  if (is.na(as.numeric(phy$node.label[2])))
    stop("node labels can't be converted to numeric values")
  if (is.null(phy$node.label))
    stop("the tree has no node labels")
  ind <- which(phy$edge[, 2] > length(phy$tip.label))[as.numeric(phy$node.label[2:length(phy$node.label)]) < tol]
  n <- length(ind)
  if (!n) 
    return(phy)
  foo <- function(ancestor, des2del) {
    wh <- which(phy$edge[, 1] == des2del)
    for (k in wh) {
      if (phy$edge[k, 2] %in% node2del) 
        foo(ancestor, phy$edge[k, 2])
      else phy$edge[k, 1] <<- ancestor
    }
  }
  node2del <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node2del) 
      next
    foo(anc[i], node2del[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
  phy
}
