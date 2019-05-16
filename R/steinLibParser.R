globalVariables(c("edgeVerbosity","edgeInfo","rowIndex","terminalVerbosity","terminalInfo"))
#' @title extract out the edge/arc and terminal information from a STP file
#' @description A function used to read in SteinLib formatted files. It is not the intention of the author for general users to call this directly.
#' @param fileName The file name of the Steiner problem in SteinLib format
#' @return list A list of two components: 'terminalNodes', which details node terminal information and 'edges', which
#' @references \url{http://steinlib.zib.de/format.php}
#' @importFrom stringr str_match str_extract str_split
#' @importFrom stats na.exclude
parseSTPfile <- function(fileName){

  file.connection <- file(fileName,open="r")
  file.lines <- readLines(file.connection)
  close(file.connection)

  nEdges <- file.lines %>% str_match('Edges\\s+\\d+') %>% na.omit() %>% str_extract('\\d+') %>% as.integer()
  nArcs <- file.lines %>% str_match('Arcs\\s+\\d+') %>% na.omit() %>% str_extract('\\d+') %>% as.integer()
  nNodes <- file.lines %>% str_match('Nodes\\s+\\d+') %>% na.omit() %>% str_extract('\\d+') %>% as.integer()
  nTerminals <- file.lines %>% str_match('Terminals\\s+\\d+') %>% na.omit() %>% str_extract('\\d+') %>% as.integer()

  if(length(nEdges) > 0 & length(nArcs) == 0){nArcs <- 0}
  if(length(nArcs) > 0 & length(nEdges) == 0){nEdges <- 0}

  if(length(nEdges) > 1){ stop('MStTP parser failed. There should only be one line that declares the number of edges')  }
  if(length(nArcs) > 1){ stop('MStTP parser failed. There should only be one line that declares the number of arcs')  }
  if(length(nNodes) > 1){ stop('MStTP parser failed. There should only be one line that declares the number of nodes')  }
  if(length(nTerminals) > 1){ stop('MStTP parser failed. There should only be one line that declares the number of terminals')  }

  graphEdgesDT <- file.lines %>% str_match('^[AE]?\\s+\\d+\\s+\\d+\\s*\\d*') %>% na.exclude() %>% str_split('\\s+') %>% data.table(edgeInfo = .)
  graphEdgesDT[,edgeVerbosity := length(edgeInfo[[1]]),by=1:nrow(graphEdgesDT)]

  if(graphEdgesDT[,unique(edgeVerbosity)] %>% length > 1){ stop('Differing amounts of information in each edge line!')}

  if(graphEdgesDT[,unique(edgeVerbosity)] == 3){
    #For cases where there are no edgeweights

    graphEdgesDT[,`:=`(edgeType=edgeInfo[[1]][1],
                       start = as.integer(edgeInfo[[1]][2]),
                       end = as.integer(edgeInfo[[1]][3])
                      )
    , by = 1:nrow(graphEdgesDT)]

  }else if(graphEdgesDT[,unique(edgeVerbosity)] == 4){
    #For cases where there are edgeweights

    graphEdgesDT[,`:=`(edgeType=edgeInfo[[1]][1],
                       start = as.integer(edgeInfo[[1]][2]),
                       end = as.integer(edgeInfo[[1]][3]),
                       weight = as.numeric(edgeInfo[[1]][4])
                      )
                , by = 1:nrow(graphEdgesDT)]

  }else{
    stop('Parsing of edges failed!')
  }

  graphEdgesDT[,edgeInfo:=NULL]
  graphEdgesDT[,edgeVerbosity:=NULL]

  #There are some odd times that people who produce files don't use the edges and arcs field properly. So we'll just check the sum
  if(nrow(graphEdgesDT) != nEdges+nArcs){stop('Edge count wrong!')}


  # Now lets pass the terminal information - note the optional P in the rejex for terminals that must be present
  terminalSet <- grep('^TP?\\s+',file.lines,perl=TRUE,value=TRUE) %>% str_split('\\s+') %>% data.table(terminalInfo = .)
  terminalSet[,rowIndex := seq(1,nrow(terminalSet))]
  terminalSet[,terminalVerbosity := length(terminalInfo[[1]]),by=rowIndex]

  if(terminalSet[,unique(terminalVerbosity)] %>% length > 1){stop("Differing amounts of data in each entry.")}

  if(terminalSet[,unique(terminalVerbosity)] == 2){

    terminalSet[,`:=`(nodeType = terminalInfo[[1]][1],
                      nodeIndex = as.integer(terminalInfo[[1]][2])),
                by = 1:nrow(terminalSet)]
  }else if(terminalSet[,unique(terminalVerbosity)] == 3){
    #We have not only node information, but also node score info

    terminalSet[,`:=`(nodeType = terminalInfo[[1]][1],
                      nodeIndex = as.integer(terminalInfo[[1]][2]),
                      nodeScore = as.numeric(terminalInfo[[1]][3])),
                by = 1:nrow(terminalSet)]
  }else{
    stop("Unable to parse the terminal information in steinLib formatted file")
  }

  terminalSet[,terminalInfo := NULL]
  terminalSet[,terminalVerbosity := NULL]
  terminalSet[,rowIndex := NULL]

  if(nrow(terminalSet) != nTerminals){stop('Terminal count wrong!')}

  return(  list(terminalNodes = terminalSet, edges = graphEdgesDT))
}

globalVariables(c("edgeType","weight","nodeIndex"))
#' @title create an igraph representation of a MStT problem
#' @description The SteinLib resource provides a reference collection of Steiner Tree problems and solutions. Very useful for validation. This
#' function allows you to read in graphs relating to the minimum Steiner Tree problem.
#' @param MStTPstpFile The URI of the Steiner problem in SteinLib format
#' @return steinerProblem An igraph object with vertex names and a per vertex flag 'isTerminal' which specifies if a node is a terminal or not.
#'   Note that upon creation of the igraph object, node names (as defined in the STP file) and nodeIDs (as produced by igraph).
#' @references \url{http://steinlib.zib.de/format.php}
#' @export
readMStTPgraph <- function(MStTPstpFile){

  stpGraphInformation <- parseSTPfile(MStTPstpFile)

  edgesDT <- stpGraphInformation$edges
  #We're gonna create a digraph by default, so we'll add reverse arcs back if we have and edges (which are undirected)

  if(ncol(edgesDT) == 4){
    edgesDT %<>% rbind(edgesDT[edgeType == 'E',.(start = end,end=start,edgeType,weight)])
    MStTP_igraph <- graph.data.frame(edgesDT[,.(start,end,edgeWeight=weight)])

  }else if(ncol(edgesDT) == 3){
    edgesDT %<>% rbind(edgesDT[edgeType == 'E',.(start = end,end=start,edgeType)])
    MStTP_igraph <- graph.data.frame(edgesDT[,.(start,end)])

  }else{  stop('Parsing of MStTP file failed on edges\n') }

  if(edgesDT[,(length(unique(edgeType)) == 1 & unique(edgeType) == 'E')]){

    MStTP_igraph %<>% as.undirected
    #Note, this will strip out edgeWeight attributes (thanks igraph ...), but by the definition of SteinLib format,
    #there shouldn't be any edgeWeight info for E type (undirected) edges
  }

  terminalDT <- stpGraphInformation$terminalNodes

  #It's very possible for there to be nodes in the terminal list that do not have edges.
  nodes2add <- terminalDT[! nodeIndex %in% V(MStTP_igraph)$name, nodeIndex]
  for(node2add in nodes2add){
    MStTP_igraph %<>% add(vertex(node2add))
  }

  V(MStTP_igraph)$isTerminal <- FALSE
  V(MStTP_igraph)[name %in% stpGraphInformation$terminalNodes$nodeIndex ]$isTerminal <- TRUE
  #Here we assume that any terminal listed in a file is a so-called fixed terminal - i.e. it must be considered in a candidate solution

  if(length(which(V(MStTP_igraph)$isTerminal )) >= 0.3*vcount(MStTP_igraph) ){warning('Over 30% of your graph nodes are fixed terminals. Are you sure that this is what you want?')}

  return(MStTP_igraph)
}

globalVariables(c("edgeType","weight","nodeIndex","nodeScore","start","end"))
#' @title create an igraph representation of a MWCS problem
#' @description The SteinLib resource provides a reference collection of Steiner Tree problems and solutions. Very useful for validation. This function allows
#' you to read in graphs relating to maximum-weight connected subgraph problems.
#' @param MWCSstpFile The URI of the Steiner problem in SteinLib format
#' @return steinerProblem An igraph object with vertex names and a per vertex flag 'isTerminal' which specifies if a node is a terminal or not.
#'   Note that upon creation of the igraph object, node names (as defined in the STP file) and nodeIDs (as produced by igraph).
#' @references \url{http://steinlib.zib.de/format.php}
#' @export
readMWCSgraph <- function(MWCSstpFile){

  stpGraphInformation <- parseSTPfile(MWCSstpFile)

  edgesDT <- stpGraphInformation$edges
  #We're gonna create a digraph by default, so we'll add reverse arcs back if we have and edges (which are undirected)

  if(ncol(edgesDT) == 4){
    edgesDT %<>% rbind(edgesDT[edgeType == 'E',.(start = end,end=start,edgeType,weight)])
    MWCS_igraph <- graph.data.frame(edgesDT[,.(start,end,edgeWeight=weight)])

  }else if(ncol(edgesDT) == 3){
    edgesDT %<>% rbind(edgesDT[edgeType == 'E',.(start = end,end=start,edgeType)])
    MWCS_igraph <- graph.data.frame(edgesDT[,.(start,end)])

  }else{  stop('Parsing of MStTP file failed on edges\n') }

  #If we only have undirected edges, we should make the graph undirected
  if(edgesDT[,(length(unique(edgeType)) == 1 & unique(edgeType) == 'E')]){

    MWCS_igraph %<>% as.undirected
    #Note, this will strip out edgeWeight attributes (thanks igraph ...), but by the definition of SteinLib format,
    #there shouldn't be any edgeWeight info for E type (undirected) edges
  }

  terminalDT <- stpGraphInformation$terminalNodes

  if(nrow(terminalDT) < vcount(MWCS_igraph)){warning('It is recommended that you provide a node score per vertex for the MWCS problem')}

  #It's very possible for there to be nodes in the terminal list that do not have edges.
  nodes2add <- terminalDT[! nodeIndex %in% V(MWCS_igraph)$name, nodeIndex]
  for(node2add in nodes2add){
    MWCS_igraph %<>% add(vertex(node2add))
  }

  if(!'nodeScore' %in% colnames(terminalDT)){  stop('Terminal scores must be present for some of the nodes in the MWCS problem (by definition!)') }

  terminalDT[,nodeIndex := as.character(nodeIndex)]
  setkey(terminalDT,'nodeIndex')
  V(MWCS_igraph)$nodeScore <- as.numeric(NA)
  V(MWCS_igraph)[terminalDT[V(MWCS_igraph)$name,nodeIndex]]$nodeScore <- terminalDT[V(MWCS_igraph)$name,nodeScore]

  return(MWCS_igraph)
}

globalVariables(c("index","edgeType","edgeWeight","nodeType","indexA","indexB","terminalPrize"))
#' @title create an igraph representation of a MStT problem
#' @description The heinz program uses a slightly modified version of steinLib format, with no documentation. This function
#'   writes out a MWCS problem into an appropriate form for use with heinz
#' @param file.name The file name of the output file
#' @param network.igraph Graph object to be written out
#' @param remark SteinLib format has a remark field (rarely used), but feel free to add something in here
#' @param problem SteinLib format has a problem field (rarely used), but feel free to add something in here
#' @param name SteinLib format has a name field (rarely used), but feel free to add something in here
#' @return node2ID A data.table of nodeName to nodeIndex mappings
#' @references \url{http://steinlib.zib.de/format.php}
#' @references \url{https://github.com/ls-cwi/heinz}
writeMStTPfile_heinzFormat <-function(network.igraph,file.name, remark ='Nothing of note',problem = 'MStTP', name = 'Autogenerated Name'){

  network.igraph %<>% as.undirected
  terminal.nodes <- V(network.igraph)[isTerminal == TRUE]$name

  node2ID <- data.table( node = unique(V(network.igraph)$name))
  node2ID[,index := 1:nrow(node2ID)]

  edges2ID <- get.edgelist(network.igraph) %>% data.table
  setnames(edges2ID,c('V1'),c('node'))
  edges2ID <- merge(node2ID,edges2ID, by = 'node')
  setnames(edges2ID, c('node','index','V2'),c('nodeA','indexA','node'))
  edges2ID <- merge(node2ID,edges2ID, by = 'node')
  setnames(edges2ID, c('node','index'),c('nodeB','indexB'))

  edge.type <- 'E'
  edges2ID[,edgeType := edge.type]
  #To do - allow for variable edge weight
  if(is.element('edgeWeight', list.edge.attributes(network.igraph))){

    edges2ID[,edgeWeight := E(network.igraph)$edgeWeight]
  }else{

    edges2ID[,edgeWeight := 0]
  }


  node2ID[,nodeType := 'T']

  #god I hate the STP file format. But it's a standard, so best not to rock the boat on this one.
  sink(file = file.name)
  cat("33D32945 STP File, STP Format Version 1.0\n\n")

  #Comment Section
  cat("SECTION Comment\n")
  cat(paste(c('Name "',  name,'"',"\n"),collapse = ''))
  cat(paste(c('Date "',  as.character(Sys.Date()),'"',"\n"),collapse = ''))
  cat(paste(c('Creator "',system('whoami', intern=TRUE),'"',"\n"),collapse = ''))
  cat(paste(c('Remark "',remark,'"',"\n"),collapse = ''))
  cat(paste(c('Problem "',problem,'"',"\n"),collapse = ''))
  cat("END\n\n")

  cat("SECTION Graph\n")
  cat(paste(c('Nodes',nrow(node2ID),"\n"),collapse = ' '))
  if(is.directed(network.igraph)){

    cat(paste(c('Arcs',nrow(edges2ID),"\n"),collapse = ' '))
    write.table(edges2ID[order(indexA)][,.(edgeType,indexA,indexB,edgeWeight)],file='', quote = FALSE, row.names=FALSE, col.names=FALSE)
  }else{

    cat(paste(c('Edges',nrow(edges2ID),"\n"),collapse = ' '))
    write.table(edges2ID[order(indexA)][,.(edgeType,indexA,indexB)],file='', quote = FALSE, row.names=FALSE, col.names=FALSE)
  }

  cat("END\n\n")

  node2ID[,terminalPrize := -1]
  node2ID[node %in% terminal.nodes, terminalPrize := 100]

  cat("SECTION Terminals\n")
  cat(paste(c('Terminals',nrow(node2ID),"\n"),collapse = ' '))
  write.table(node2ID[, .(nodeType,index,terminalPrize)],file='', quote = FALSE, row.names=FALSE, col.names=FALSE)
  cat("END\n\n")

  cat("EOF\n")
  sink()

  return(node2ID[,.(node,index,terminalPrize)])
}

