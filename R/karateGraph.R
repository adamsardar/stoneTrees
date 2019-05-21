#' Zachary's famous karate graph
#'
#' An undirected graph describing a social network of friendships between 34 members of a karate club at a US university in the 1970's. It serves as a small graph to test and demonstrate Steiner tree problems.
#'
#' @docType data
#' @usage data(karateGraph)
#'
#' @format An igraph object with 34 nodes and 78 edges
#' \describe{
#'   \item{name}{Unique identifiers for each vertex}
#'   \item{isTerminal}{A boolean flag detailing as to whether a vertex is a fixed node or not}
#' }
#' @source \url{https://networkdata.ics.uci.edu/data.php?id=105}
#' @references W. W. Zachary, An information flow model for conflict and fission in small groups, Journal of Anthropological Research 33, 452-473 (1977)
"karateGraph"
