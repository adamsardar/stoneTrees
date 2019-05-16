#' Lymphoma example from DIMAC11 challenge
#'
#' An udirected graph representing a Maximum-Weight Connected Subgraph (MWCS) problem from a genetic context
#' Instance used in the DIMACS11 Steiner Tree problems challengs (2015).
#'
#' @docType data
#' @usage data(lymphomaGraph)
#'
#' @format An igraph object with 2034 nodes and 7756 edges
#' \describe{
#'   \item{name}{Unique identifiers for each vertex - note, not an index despite being character numbers}
#'   \item{nodeScore}{A score or weight on the node. Produced using the BioNet algroithm using diffuse large B-cell lymphoma data}
#' }
#' @source \url{http://dimacs11.zib.de/downloads.html}
"lymphomaGraph"
