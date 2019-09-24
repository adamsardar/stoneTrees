#' gene42 example from SteinLib
#'
#' A directed graph representing a Steiner arborescence problem (directed Steiner problems) from a genetic context. This corresponds to gene42 from the steinLib dataset.
#'
#' Note: it is only advised to use this, rather difficult, Steiner tree problem in conjunction with an advanced MILP solver like CPLEX.
#'
#' @docType data
#' @usage data(gene42_igraph)
#'
#' @format An igraph object with 335 nodes and 456 edges
#' \describe{
#'   \item{name}{Unique identifiers for each vertex - note, not an index despite being character numbers}
#'   \item{isTerminal}{A boolean flag detailing as to whether a vertex is a fixed node or not}
#'   \item{edgeWeight}{A belief score between 0 and 1 indicating how strong the evidence is for an interaction}
#' }
#' @source \url{http://steinlib.zib.de/showset.php?GENE}
"gene42_igraph"
