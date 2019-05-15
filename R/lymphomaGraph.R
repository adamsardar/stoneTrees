#' Lymphoma example from SteinLib
#' 
#' An udirected graph representing a Maximum-Weight Connected Subgraph (MWCS) problem from a genetic context
#'
#' @docType data
#' @usage data(lymphomaGraph)
#'
#' @format An igraph object with 2034 nodes and 7756 edges
#' \describe{
#'   \item{name}{Unique identifiers for each vertex - note, not an index despite being character numbers}
#'   \item{nodeScore}{A score or weight on the node. Produced using the BioNet algroithm using diffuse large B-cell lymphoma data}
#' }
#' @source \url{http://bioconductor.org/packages/release/data/experiment/html/DLBCL.html}
"lymphomaGraph"