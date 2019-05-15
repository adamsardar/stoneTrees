globalVariables(names=c(".","name"))
#' @title create a module of the set of all-pairs shortest paths (APSP)
#' @description a convenience function for generating APSP for a set of seeds - not that this is ALL paths of the shortest length, not just one
#' @param seeds A set of seen nodes. These should correspond to the 'name' field of the searchGraph
#' @param searchGraph An igraph object representation of the network within which to search from shortest-path pairs
#' @param omitNA Logical flag to dictate if NA values should be stripped from input (default: TRUE)
#' @return An igraph representation of the APSP sampling
#' @importFrom ensurer ensure
#' @importFrom stringr str_c
#' @importFrom stats na.omit
#' @export
calculateAPSPmodule <- function(seeds,searchGraph,omitNA=TRUE){
  
  #Capture search graph name
  searchGraphName <- deparse(substitute(searchGraph))
  
  searchGraph %<>% validateSearchGraph
  
  ifelse(omitNA,
         seeds %<>% na.omit(seeds),
         if(any(is.na(V(searchGraph)$name))) warning("NA nodes found in input seed list but 'omitNA' set to FALSE - this can lead to unstable beahviour") )
  
  seedsNotInGaph <- setdiff(seeds,V(searchGraph)$name)
  seedsInGaph <- intersect(seeds,V(searchGraph)$name)
  
  if(length(seedsNotInGaph) == length(seeds)){stop("None of the seeds were found in graph!")}
  if(length(seedsNotInGaph) > 0){ warning(str_c(c("Not all seeds found in graph: ",seedsNotInGaph)))}
  
  APSP_vertices <- all_shortest_paths(searchGraph,
                                      from = V(searchGraph)[name %in% seedsInGaph],
                                      to = V(searchGraph)[name %in% seedsInGaph],
                                      mode="all") %>% .$res %>% unlist %>% as.integer %>% unique %>% na.omit
  process_APSP_module <- induced_subgraph(searchGraph,V(searchGraph)[APSP_vertices])
  
  V(process_APSP_module)$isSeed <- FALSE
  V(process_APSP_module)[name %in% seeds]$isSeed <- TRUE
  
  process_APSP_module %<>% set_graph_attr("SamplingMethod","All-Pairs Shortest Paths Of Seed Nodes")
  process_APSP_module %<>% set_graph_attr("SearchNetwork",searchGraphName)    
  
  return(process_APSP_module)
}