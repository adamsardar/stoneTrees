globalVariables(c("score1","nodeA","nodeB","node","V1","V2"))
#' @title use the heinz algorithm to find a solution to the maximum weight connected subgraph problem
#' @description heinz ("heavy induced subgraphs") computes solutions to the maximum-weight connected subgraph problem.
#' @param searchGraph_igraph The igraph representation of the graph that we look to search over with the MWCS solver
#' @param induceSubgraph Flag to dictate if you would prefer a tree or the induced subgraph of the solution
#' @param n.threads Specifies the maximum number of threads that the heinz binary can use
#' @param verbosity A flag that specifies if you want information on what the function is doing
#' @param timeLimit The maximum number of seconds that the heinz binary should be allowed to run for (default: 600 seconds)
#' @return igraphModule An igraph representation of the solution (with induced edges from the original search graph if requested)
#' @references An algorithmic framework for the exact solution of the prize-collecting Steiner tree problem. Mathematical Programming, Series B, 105(2-3):427-449, 2006, by I. Ljubic, R. Weiskircher, U. Pferschy, G. Klau, P. Mutzel, and M. Fischetti.
#' @references M. T. Dittrich, G. W. Klau, A. Rosenwald, T. Dandekar and T. Mueller (2008) Identifying functional modules in protein-protein interaction networks: an integrated exact approach. (ISMB2008) Bioinformatics 24: 13. i223-i231 Jul.
#' @references D. Beisser, G. W. Klau, T. Dandekar, T. Mueller and M. Dittrich (2010) BioNet: an R-package for the Functional Analysis of Biological Networks. Bioinformatics 26:08. 1129-1130 Apr.
#' @importFrom utils write.table
#' @export
findHeinzMWCS <- function(searchGraph_igraph, n.threads = 2, induceSubgraph = FALSE, verbose = TRUE, timeLimit = 600){

  interactomeName <- deparse(substitute(searchGraph_igraph))

  #Heinz doesn't like undirected graphs - clearly there's some CPLEX constraint that can end up in a loop.
  searchGraph_UndirectedIgraph <- searchGraph_igraph %>% simplify %>% as.undirected

  edgelist <- searchGraph_UndirectedIgraph %>% get.edgelist %>% data.table
  edgelist[,score1 := 0]
  setnames(edgelist, c('V1','V2'),c('nodeA','nodeB'))

  nodelist <- data.table(node = V(searchGraph_UndirectedIgraph)$name, score1 = V(searchGraph_UndirectedIgraph)$nodeScore)
  nodelist %<>% .[,.(score1 = max(score1)),by='node']

  edgelist %<>% .[nodeA %in% nodelist$node | nodeB %in% nodelist$node]
  nodelist %<>% .[node %in% edgelist$nodeA | node %in% edgelist$nodeB]

  setnames(edgelist, c('nodeA'), c('#nodeA'))
  setnames(nodelist, c('node'), c('#node'))

  MWCS.filename <- tempfile(pattern = "MWCSinput")
  MWCS.result <- tempfile(pattern = "MWCSresult")

  nodelist %>% write.table(file = paste(c(MWCS.filename,'nodes'),collapse = '.'), quote = FALSE, row.names = FALSE)
  edgelist %>% unique.data.frame %>% write.table(file = paste(c(MWCS.filename,'edges'),collapse = '.'), quote = FALSE, row.names = FALSE)

  if(verbose){ message("Temp input file stub is: ",MWCS.filename) }

  system(paste(c("heinz -n ",MWCS.filename,".nodes -e ",MWCS.filename,".edges -m ",n.threads," -v ",verbosity," -o ",MWCS.result," -t ",timeLimit),collapse=''),
         ignore.stdout = TRUE, ignore.stderr = !verbose )

  if(verbose){ message("Temp output file is: ",MWCS.result) }

  if(! file.exists(MWCS.result)){
    warning("No solution found by heinz within allowed time. Returning empty graph.")
    return(graph.empty())
  }

  cplex.result <- fread(MWCS.result, header = FALSE)[! V1 %like% '#']
  MWCS.nodes <- cplex.result[V2 != 'NaN', V1]

  MWCS.module <- induced.subgraph(searchGraph_igraph,vids = MWCS.nodes)

  file.remove(paste(c(MWCS.filename,'nodes'), collapse = '.'))
  file.remove(paste(c(MWCS.filename,'edges'), collapse = '.'))
  file.remove(MWCS.result)

  MWCS.module %<>% set_graph_attr("SearchNetwork",interactomeName)

  if(induceSubgraph){

    return(MWCS.module)
  }else{

    return( minimum.spanning.tree(MWCS.module) )
  }
}


globalVariables("isTerminal")
#' @title use the heinz algorithm to find a solution to the minimum steiner problem
#' @description heinz ("heavy induced subgraphs") computes solutions to the minimum Steiner Tree Problem.
#' @inheritParams findHeinzMWCS
#' @param searchGraph.igraph The igraph representation of the graph that we look to search over with the MStTP solver. Note, this graph should have vertex attributes of 'isTerminal'
#' @return igraphModule An igraph representation of the solution (with induced edges from the original search graph if requested)
#' @family minimumSteinerTrees
#' @references An algorithmic framework for the exact solution of the prize-collecting Steiner tree problem. Mathematical Programming, Series B, 105(2-3):427-449, 2006, by I. Ljubic, R. Weiskircher, U. Pferschy, G. Klau, P. Mutzel, and M. Fischetti.
#' @references M. T. Dittrich, G. W. Klau, A. Rosenwald, T. Dandekar and T. Mueller (2008) Identifying functional modules in protein-protein interaction networks: an integrated exact approach. (ISMB2008) Bioinformatics 24: 13. i223-i231 Jul.
#' @references D. Beisser, G. W. Klau, T. Dandekar, T. Mueller and M. Dittrich (2010) BioNet: an R-package for the Functional Analysis of Biological Networks. Bioinformatics 26:08. 1129-1130 Apr.
#' @references \url{http://dimacs11.cs.princeton.edu/contest/challenge-results.pdf}
#' @export
findHeinzMStTP <- function(searchGraph.igraph, n.threads = 2, induceSubgraph = FALSE, verbose = TRUE, timeLimit=600){

  interactomeName <- deparse(substitute(searchGraph.igraph))

  terminalNodes <- V(searchGraph.igraph)[isTerminal == TRUE]$name

  V(searchGraph.igraph)$nodeScore <- -1
  V(searchGraph.igraph)[name %in% terminalNodes]$nodeScore <- 1000

  #It turns out that converting the problem to the MWCS and solving that is actually better than the dedicated MStTP solver in heinz
  sol_MSTtP <- findHeinzMWCS(searchGraph.igraph,n.threads = n.threads,
                induceSubgraph = induceSubgraph,
                verbose = verbose,
                timeLimit = timeLimit)

  sol_MSTtP %<>% set_graph_attr("SearchNetwork",interactomeName)

  return(sol_MSTtP)
}
