globalVariables(c("condensedNode","potentialTerminal","segmentID",".vertexID"))

#' Collapse neighbouring terminal and potential terminal nodes together
#'
#' Part of the Steiner Tree solving method involves enumerating all vs all node separators of nodes in differing
#' components. This scales very, very badly and for some runs can actually be the major computational bottleneck. It's
#' also unnecessary as neighbouring prize nodes and terminals will always be in the resultant graph.
#'
#' This function groups together potential terminals and adds their relevant node attributes together sensibly. For use
#' with steiner solvers in the stoneTrees package.
#'
#' @param graphToCondense igraph network that we wish to contract
#' @param condensedNodeSep (optional) Specify the value to separate node names with
#'
#' @return A graph with neighbouring terminal and potential terminal nodes collapsed to the same node
#' @seealso uncondenseGraph
condenseSearchGraph = function(graphToCondense, condensedNodeSep = ";"){

  check_network(graphToCondense, singleWeakComponent = FALSE, isDirected = FALSE)
  check_string(condensedNodeSep, allow_empty = FALSE)
  stopifnot("nodeNameSep must not be in any node names already" = !any(V(graphToCondense)$name %like% condensedNodeSep))

  V(graphToCondense)$.vertexID = 1:vcount(graphToCondense)

  originalGraph = graphToCondense

  if(!"isTerminal" %in% vertex_attr_names(graphToCondense)){ V(graphToCondense)$isTerminal = FALSE }
  if(!"nodeScore" %in% vertex_attr_names(graphToCondense)){ V(graphToCondense)$nodeScore = -1 }

  V(graphToCondense)$condensedNode = FALSE

  if(length(V(graphToCondense)[nodeScore > 0 | isTerminal]) == 0){stop("No potential terminals in input graph!")}

  V(graphToCondense)$potentialTerminal = FALSE
  V(graphToCondense)[nodeScore > 0 | isTerminal]$potentialTerminal = TRUE

  #Induce a subgraph of potential terminals, decompose into the connected components and then assign segmentIDs to the blocks
  potentialTerminalGraph = induced_subgraph(graphToCondense, V(graphToCondense)[potentialTerminal])

  nodeDT = as_data_frame(graphToCondense, what = "vertices") %>% data.table(key = ".vertexID")
  nodeDT[, segmentID := .vertexID]

  lapply(decompose(potentialTerminalGraph, mode = "weak", min.vertices = 2),
         function(g){ nodeDT[V(g)$.vertexID, segmentID := min(segmentID)]; return(invisible(nodeDT[]))})

  collapseNodeName = function(x){return(paste0(x, collapse = condensedNodeSep))}
  isCollapsed = function(x){return(length(x) > 1)}

  condensedGraph = contract(graphToCondense,
                              mapping = nodeDT[V(graphToCondense)$.vertexID, as.integer(as.factor(segmentID)) ], #Create a integer sequence along the segmentID variable
                              vertex.attr.comb = list(nodeScore = sum,
                                                      name = collapseNodeName,
                                                      isTerminal = any,
                                                      condensedNode = isCollapsed,
                                                      .vertexID = collapseNodeName,
                                                      "ignore") )

  condensedGraph %<>% set_graph_attr("originalGraph", originalGraph) # Keep the original graph present as an attribute. Useful in the uncondensation.
  condensedGraph %<>% set_graph_attr("nodeNameSep", condensedNodeSep)

  return(condensedGraph)
}

#' Decondense a graph following a procedure to group together prize nodes. If the input was not condensed, then we
#' just return it back
#'
#' @param condensedGraph An igraph representation of the condensed graph (where prize nodes). The result of calling 'condenseGraph()'
#'
#' @return uncondensed graph
#' @seealso condenseSearchGraph
uncondenseGraph = function(condensedGraph){

  if(! 'condensedNode' %in% vertex_attr_names(condensedGraph)){
    #No condensation occured, so we'll just return the original graph
    return(condensedGraph)
  }

  stopifnot("condensedGraph must be an igraph object resultant from running 'condenseSearchGraph'" = 
    'condensedNode' %in% vertex_attr_names(condensedGraph))
  stopifnot("condensedGraph must be an igraph object resultant from running 'condenseSearchGraph'" = 
    'originalGraph' %in% graph_attr_names(condensedGraph))
  stopifnot("condensedGraph must be an igraph object resultant from running 'condenseSearchGraph'" = 
    'nodeNameSep' %in% graph_attr_names(condensedGraph))

  check_logical(is.logical(V(condensedGraph)$condensedNode))

  nodeNameSep = graph_attr(condensedGraph)$nodeNameSep
  originalGraph = graph_attr(condensedGraph,"originalGraph")

  condensedGraphVertexIDs = V(condensedGraph)[condensedNode == TRUE]$.vertexID %>%
    str_split(nodeNameSep) %>%
    unlist %>%
    c(V(condensedGraph)[condensedNode == FALSE]$.vertexID) %>%
    as.integer

  condensedGraph %<>% delete_graph_attr("originalGraph")

  graph_attr(originalGraph) = graph_attr(condensedGraph)

  returnGraph = induced_subgraph(originalGraph,V(originalGraph)[.vertexID %in% condensedGraphVertexIDs])

  return(  delete_vertex_attr(returnGraph, ".vertexID"))
}
