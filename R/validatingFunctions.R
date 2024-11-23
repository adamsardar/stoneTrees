#' @import igraph
check_network = function(network2validate, singleWeakComponent = TRUE, isDirected = FALSE){

  stopifnot("Input network must be an igraph object" = is.igraph(network2validate))

  check_bool(singleWeakComponent)
  stopifnot("Input network must a single connected component (consider using igraph::decompose?)" =
             length(decompose(network2validate, mode = "weak")) == 1 | !singleWeakComponent)

  check_bool(isDirected)
  stopifnot("Expecting is.directed() output of graph to be equal to isDirected" = 
            is.directed(network2validate) == isDirected)

  if("nodeScore" %in% vertex_attr_names(network2validate)){

    stopifnot("nodeScore values for input graph must be numeric vectors" = is.numeric(V(network2validate)$nodeScore),
              "nodeScore values must have no missing (NA) values" = all(!is.na(V(network2validate)$nodeScore)) )
  }

  if("isTerminal" %in% vertex_attr_names(network2validate)){
    
    stopifnot("isTerminal values for input graph must be logical vectors" = is.logical(V(network2validate)$isTerminal),
              "isTerminal values must have no missing (NA) values" = all(!is.na(V(network2validate)$isTerminal)) )
  }

  invisible(network2validate)
}