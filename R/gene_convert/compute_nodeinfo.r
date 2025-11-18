#' Compute some useful metrics for PPI nodes
#' @param g An igraph object containing PPI information.
#' @param weight_attr Character; The attribute of weight in this PPI object.
#' @importFrom igraph E V edge_attr_names edge_attr
#' @importFrom igraph degree strength betweenness closeness
#' @importFrom igraph eigen_centrality page_rank coreness
#' @importFrom igraph transitivity eccentricity articulation_points
#' @return The input \code{igraph} object with additional vertex
#' @export

compute_nodeinfo <- function(g, weight_attr = "score") {
  stopifnot(inherits(g, "igraph"))

  ## PPI score
  w <- NULL
  if (!is.null(weight_attr) && weight_attr %in% igraph::edge_attr_names(g)) {
    w <- igraph::edge_attr(g, "score")
  }

  w_dist <- NULL
  if (!is.null(w)) {
    w_dist <- 1 / pmax(w, .Machine$double.eps)
  }

  ## compute node info
  message("Calculating degree / strength ...")
  V(g)$degree <- igraph::degree(g, mode = "all")
  V(g)$strength <- if (!is.null(w)) igraph::strength(g, mode = "all", weights = w) else NA_real_

  message("Calculating betweenness (unweighted) ...")
  V(g)$betweenness <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)

  if (!is.null(w_dist)) {
    message("Calculating betweenness (weighted) ...")
    V(g)$betweenness_w <- igraph::betweenness(g, directed = FALSE,
                                              weights = w_dist,
                                              normalized = TRUE)
  } else {
    V(g)$betweenness_w <- NA_real_
  }

  message("Calculating closeness (unweighted) ...")
  V(g)$closeness <- igraph::closeness(g, normalized = TRUE)

  if (!is.null(w_dist)) {
    message("Calculating closeness (weighted) ...")
    V(g)$closeness_w <- igraph::closeness(g, normalized = TRUE,
                                          weights = w_dist)
  } else {
    V(g)$closeness_w <- NA_real_
  }

  message("Calculating eigenvector centrality ...")
  V(g)$eigen_centrality <- igraph::eigen_centrality(g, weights = w)$vector

  message("Calculating PageRank ...")
  V(g)$pagerank <- igraph::page_rank(g, weights = w)$vector

  message("Calculating coreness (k-core) ...")
  V(g)$coreness <- igraph::coreness(g, mode = "all")

  message("Calculating local clustering coefficient ...")
  V(g)$clustering_coef <- igraph::transitivity(g, type = "local", isolates = "zero")

  message("Calculating eccentricity ...")
  V(g)$eccentricity <- igraph::eccentricity(g)

  message("Detecting articulation points ...")
  art <- igraph::articulation_points(g)
  V(g)$is_articulation <- FALSE
  V(g)$is_articulation[art] <- TRUE

  ## compute MCC
  message("Calculating MCC ...")
  g <- compute_MCC(g)

  return(g)
}
