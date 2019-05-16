#' Collapse neighbouring terminal and potential terminal nodes together
#'
#' Part of the Steiner Tree solving method involves enumerating all vs all node separators of nodes in differing
#' components. This scales very, very badly and for some runs can actually be the major computational bottleneck. It's
#' also unnecessary as neighbouring prize nodes and terminals will always be in the resultant graph.
#'
#' @param graphToCondense igraph network that we wish to contract
#' @param condensedNodeSep (optional) Specify the value to separate node names with
#'
#' @return A graph with neighbouring terminal and potential terminal nodes collapsed to the same node
#' @seealso uncondenseGraph
#' @importFrom stringr str_c
#' @importFrom ensurer ensure
condenseSearchGraph <- function(graphToCondense, condensedNodeSep = ";"){

  validateIsNetwork(graphToCondense)

  condensedNodeSep %<>% ensure(is.character,
                          all(! V(graphToCondense)$name %like% .),
                          err_desc = "nodeNameSep must not be in any node names")

  veretexAttr <- vertex_attr_names(graphToCondense)
  originalGraph <- graphToCondense
  graphToCondense %<>% set_graph_attr("originalGraph", originalGraph)

  if(!"isTerminal" %in% veretexAttr){ V(graphToCondense)$isTerminal <- FALSE }
  if(!"nodeScore" %in% veretexAttr){ V(graphToCondense)$nodeScore <- -1 }

  positiveScoringNodeGraph <- induced_subgraph(graphToCondense, V(graphToCondense)[nodeScore > 0 | isTerminal])

  if(vcount(positiveScoringNodeGraph) == 0){stop("No prize terminals in input graph!")}

  graphSegments <- decompose.graph(positiveScoringNodeGraph)
  V(graphToCondense)$condensedNode <- FALSE

  for(graphSegment in graphSegments){

    if(vcount(graphSegment) == 1){ next}

    nNodes <- vcount(graphToCondense)
    segmentNodes <- V(graphSegment)

    segmentProperties <- list(name = vertex_attr(graphSegment,name="name") %>% str_c(collapse = condensedNodeSep),
                              condensedNode = TRUE)

    segmentProperties$isTerminal <- any(V(graphSegment)$isTerminal)

    totalSegmentScore <- sum(V(graphSegment)$nodeScore)
    segmentProperties$nodeScore <- totalSegmentScore

    #Grab all neighbours of nodes to condense
    neighbouringNodes <- ego(graphToCondense, order=1, nodes = V(graphSegment)$name) %>% unlist %>% unique %>%
      setdiff( ego(graphToCondense, order=0, nodes = V(graphSegment)$name) %>% unlist)

    graphToCondense %<>% add_vertices(1,attr = segmentProperties)

    #Create node pairs so as to add as edges
    neighbourPairs <- rep(nNodes+1,2*length(neighbouringNodes))
    neighbourPairs[seq(2,2*length(neighbouringNodes),by=2)] <- neighbouringNodes

    graphToCondense %<>% add_edges(neighbourPairs)
    graphToCondense %<>% delete_vertices(V(graphSegment)$name)
  }

  return(graphToCondense)
}

globalVariables(c("condensedNode"))
#' Decondense a graph following a procedure to group together prize nodes. If the input was not condensed, then we
#' just return it back
#'
#' @param condensedGraph An igraph representation of the condensed graph (where prize nodes )
#' @param originalGraph (optional) The original graph that the condensation procedure was run on. Should probably be left
#' empty as the function will evaluate and draw in the original network
#' @param nodeNameSep (optional) If a different node-name separator was used in the graph condensation step
#'
#' @return uncondensed graph
#' @seealso condenseSearchGraph
#' @importFrom ensurer ensure
uncondenseGraph <- function(condensedGraph,originalGraph = NULL,nodeNameSep = ";"){

  if(! 'condensedNode' %in% list.vertex.attributes(condensedGraph)){
    #No condensation occured, so we'll just return the original graph
    return(condensedGraph)
  }

  condensedGraph %<>%
    #validateSearchGraph %>%
                      ensure(
                          'condensedNode' %in% list.vertex.attributes(.),
                          'originalGraph' %in% list.graph.attributes(.),
                          all(is.logical(V(.)$condensedNode)),
                          err_desc = "condensedGraph must be an igraph object resultant from running 'condenseSearchGraph'")

  nodeNameSep %<>% ensure(is.character,
                          all(V(condensedGraph)[condensedNode == TRUE]$name %like% .),
                          all(! V(condensedGraph)[condensedNode == FALSE]$name %like% .),
                          err_desc = "nodeNameSep must be only in condensed node names and not in any other node")

  if(is.null(originalGraph)){
    originalGraph <- graph_attr(condensedGraph,"originalGraph")
  }

  condensedGraphNodeNames <- V(condensedGraph)[condensedNode == TRUE]$name %>%
    str_split(nodeNameSep) %>%
    unlist %>%
    c(V(condensedGraph)[condensedNode == FALSE]$name)

  condensedGraph %<>% delete_graph_attr("originalGraph")
  graph_attr(originalGraph) <- graph_attr(condensedGraph)

  returnGraph <- induced_subgraph(originalGraph,V(originalGraph)[name %in% condensedGraphNodeNames])

  return(returnGraph)
}
