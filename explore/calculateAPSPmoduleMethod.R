library(purrr)
library(igraph)
library(doParallel)
library(foreach)
library(magrittr)

calculateAPSPmoduleETXMethod <- function(seeds,
                                         searchGraph,
                                         omitNA = TRUE,
                                         mode = "all",
                                         useDepth = FALSE,
                                         maxNetworkSize = 1000,
                                         depth = 2){
  
  
  if(omitNA){ seeds %<>% na.omit(seeds) }else{
    
    if(any(is.na(seeds))){ warning("NA nodes found in input seed list but 'omitNA' set to FALSE - this can lead to
                                 unstable behaviour") } }
  
  seedsNotInGraph <- setdiff(seeds, V(searchGraph)$name)
  
  seedsInGraph <- intersect(seeds, V(searchGraph)$name)
  
  if (length(seedsNotInGraph) == length(seeds)) {
    stop("None of the seeds were found in graph!")
  }
  
  if (length(seedsNotInGraph) > 0) {
    warning("Not all seeds found in graph: ", seedsNotInGraph)
  }
  
  
  
  APSP_Paths <- foreach(seed = iter(seedsInGraph))%dopar%{
    
    all_shortest_paths(searchGraph, from = seed, to = V(searchGraph)[name %in% seedsInGraph], mode = mode) %>%
      .$res } %>% unlist(recursive = FALSE)
  
  APSP_PathsLength <- map_dbl(APSP_Paths, ~{ length(.x)-1 } )
  
  APSP_PathsLengthMin <- min(APSP_PathsLength)
  
  APSP_PathsLengthMax <- max(APSP_PathsLength)
  
  if(!useDepth){
    
    APSP_vertices <- map(APSP_PathsLengthMin:APSP_PathsLengthMax, ~{  APSP_Paths[APSP_PathsLength <= .x] %>%
        unlist %>%
        as.integer %>%
        unique %>%
        na.omit } )
    
    maxDepth <- c(APSP_PathsLengthMin:APSP_PathsLengthMax)[max(which(map_int(APSP_vertices, length) <= maxNetworkSize))]
    
    APSP_verticesFiltered <- APSP_vertices[max(which(map_int(APSP_vertices, length) < maxNetworkSize))]
    
    process_APSP_module <- induced_subgraph(searchGraph, V(searchGraph)[APSP_verticesFiltered[[1]]])
    
    process_APSP_module %<>% set_graph_attr("maxDepth", as.double(maxDepth))
    
    process_APSP_module %<>% set_graph_attr("maxNetworkSize", as.double(maxNetworkSize))
    
  }else{ APSP_verticesFiltered <- APSP_Paths[APSP_PathsLength <= depth] %>%
    unlist %>%
    as.integer %>%
    unique %>%
    na.omit
  
  process_APSP_module <- induced_subgraph(searchGraph, V(searchGraph)[APSP_verticesFiltered])
  process_APSP_module %<>% set_graph_attr("maxDepth", as.double(depth))
  
  process_APSP_module %<>% set_graph_attr("maxNetworkSize", NULL)
  
  }
  
  V(process_APSP_module)$isSeed <- FALSE
  
  V(process_APSP_module)[name %in% seeds]$isSeed <- TRUE
  
  process_APSP_module %<>% set_graph_attr("samplingMethod", "APSP")
  
  if(!is_connected(process_APSP_module)){ warning( "APSP_module is disconnected, you might want to focus on the biggest connected component of that module" ) }
  
  return(process_APSP_module)
  
}
