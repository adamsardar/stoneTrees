#' Collect degenerate and sub-optimal solutions to Steiner problems (MStTP or MWCS) with uniform or no edge weights.
#'
#' Rather than find just a single solution to a MStTP/MWCS, one can populate a solution pool with multiple degenerate/tolerable solutions.
#'
#' This class is derived from *nodeCentricSteinerTreeProblem*: all methods available in the superclass are available here. The difference is that after
#' each acceptable solution is found, the solution is a.) stored in a solution pool and b.) used to generate a 'novelty' constraint on future solutions.
#'
#' @docType class
#' @format R6Class \code{subOptimalSteinerProblem} Construct an object representation of a multiple-solution Steiner tree/minimum weight connected subgraph (MWCS) problem
#'
#' @section methods:
#'
#' Alongisde those for *nodeCentricSteinerTreeProblem*
#' \describe{
#'    \item{\code{new(network, solverChoice = chooseSolver(), verbose = TRUE, presolveGraph = TRUE, solverTimeLimit = 300, solverTrace = as.integer(verbose), solutionTolerance = 0)}}{Constructor for the subOptimalSteinerProblem class. Alongside the arguments for the super-class constructor, there is also 'solutionTolerence', which instructs the object as to the gap between optimal and observed solution that is acceptable.}
#'    \item{\code{identifyMultipleSteinerSolutions(maxItr = 10)}}{Add solutions to the solution pool. maxItr is an argument dictating the number of runs through the obtimsation procedure.}
#'    \item{\code{getSolutionPoolGraphs(collapseSols = TRUE)}}{Either return a list of solutions within tolerence (collapseSols = FALSE) or pool all solutions together and return a single graph (collapseSols = TRUE, defaults)}
#'    \item{\code{getSolutionPoolScores()}}{Compute the scores of the solutions in the solution pool. These are in the same order as the list of graphs returned by $getSolutionPoolGraphs(FALSE)}
#'    \item{\code{getOptimumScore()}}{Returns the optimum score from solutions in the solution pool}
#'    \item{\code{getSolutionTolerance()}}{Retreive the tolerence that permits a solution to be added to the solution pool in future calls to $identifyMultipleSteinerSolutions()}
#'    \item{\code{setSolutionTolerance(x)}}{Alter the tolerence that permits a solution to be added to the solution pool in future calls to $identifyMultipleSteinerSolutions()}
#' }
#'
#' @examples
#' library(igraph)
#'
#'  # Maximum-Weight Connected Subgraph (MWCS) - find sub-optimal solutions
#'
#'  ## Vertex attribute details node costs/prizes
#'  head(V(lymphomaGraph)$nodeScore)
#'
#'  lymphoma_multiMWCS <- subOptimalSteinerProblem$new(lymphomaGraph, solutionTolerance = 0.5)
#'
#'  #Populate the solution pool with multiple solutions - notice the
#'  lymphoma_multiMWCS$identifyMultipleSteinerSolutions()
#'
#'  lymphoma_multiMWCS$getSolutionPoolGraphs(collapseSols = FALSE)
#'
#'  lymphoma_multiMWCS$getSolutionPoolScores()
#'
#'  #All solution scores are within tolerance
#'  diff(range(lymphoma_multiMWCS$getSolutionPoolScores()))
#' @references Fischetti M, Leitner M, Ljubić I, Luipersbeck M, Monaci M, Resch M, et al. Thinning out Steiner trees: a node-based model for uniform edge costs. Math Program Comput. dimacs11.cs.princeton.edu; 2017;9: 203–229.
#' @references Beisser D, Klau GW, Dandekar T, Müller T, Dittrich MT. BioNet: An R-Package for the functional analysis of biological networks. Bioinformatics. 2010;26: 1129–1130.
#' @references \url{https://en.wikipedia.org/wiki/Steiner_tree_problem}
#' @family SteinerProblemSolver
#' @seealso nodeCentricSteinerTreeProblem
#' @importFrom sets set set_union
#' @export
subOptimalSteinerProblem <- R6Class("subOptimalSteinerProblem",
                                    inherit = nodeCentricSteinerTreeProblem,
  public = list(
    
    #Overide
    initialize = function(network, solverChoice = chooseSolver(),
                          verbose = TRUE, presolveGraph = TRUE,
                          solverTimeLimit = 300, solutionTolerance = 0,
                          solverTrace = as.integer(verbose)){
      
      super$initialize(network = network,
                       solverChoice = solverChoice,
                       verbose = verbose,
                       presolveGraph = presolveGraph,
                       solverTrace = solverTrace,
                       solverTimeLimit = solverTimeLimit)
      
      self$setSolutionTolerance(solutionTolerance + 1E-10)  # Add epsilon
      
      private$setNoveltyConstraints()
      
      return(invisible(self))
    },

    getSolutionPool = function(){
      
      return(private$solutionIndciesPool)
    },
    
    getSolutionPoolGraphs = function(collapseSols = TRUE){
      
      if(collapseSols){
        
        #Ensure that the solution pool is up to date when we induce the subgraph. Since we are using a set, there is no cost to this
        return( uncondenseGraph(induced.subgraph(private$searchGraph, V(private$searchGraph)[unique(unlist( self$getSolutionPool()))])) )
      }else{
        
        return( self$getSolutionPool() %>%
                  as.list %>%
                  lapply( function(indices){ induced.subgraph(private$searchGraph, V(private$searchGraph)[indices])}) ) %>%
          lapply(uncondenseGraph)
      }
    },
    
    getSolutionPoolScores = function(){
      
      return( self$getSolutionPool() %>%
                as.list %>%
                sapply( function(indices){ super$getNodeDT()[.nodeID %in% indices, sum(nodeScore)] }) %>%
                unlist )
    },
    
    getOptimumScore = function(){ return( max( self$getSolutionPoolScores(), na.rm = TRUE) ) },
    
    getNoveltyConstraints = function(){return( private$novelSolutionsConstraint )},
    
    getSolutionTolerance = function(){return(private$tolerance)},
    
    setSolutionTolerance = function(x){ private$tolerance <- validateSinglePositiveSemiDefiniteNumeric(x) ; return(invisible(self))},
    
    identifyMultipleSteinerSolutions = function(maxItr = 10){
      
      validateSingleInteger(maxItr)
      
      self$findSingleSteinerSolution()
      private$solutionIndciesPool <- set_union(self$getSolutionPool(), sets::set(private$currentSolutionIndices) )
      
      multiSteinerItr <- 1
      
      while(multiSteinerItr <= maxItr){
        
        private$setNoveltyConstraints()
        
        super$solve()
        multiSteinerItr %<>% add(1)
        
        #add solution graph if connected, else add connectivity constraints and resolve
        if( super$isSolutionConnected() ){
          
          #If the absolute difference between scores is within tolerance, then add to pool
          if(  abs(super$getCurrentSolutionScore() - self$getOptimumScore()) <= private$tolerance ){
            
            private$solutionIndciesPool <- set_union(self$getSolutionPool(), sets::set(private$currentSolutionIndices) )
          }else{
            
            message("Next feasible solution is outside of solution tolerance! Consider increasing it with $setSolutionTolerance(x) method?")
            break()
          }
          
        }else{ super$addConnectivityConstraints() }
      }
      
      return( invisible(self) )
    }
    
  ),
  
  private = list(
    
    # Overide the superclass
    generateConstraintMatrix = function(){
      
      return(rbind(private$fixedTerminalConstraints$variables,
                   private$nodeDegreeConstraints$variables,
                   private$twoCycleConstraints$variables,
                   private$connectivityConstraints$variables,
                   private$novelSolutionsConstraint$variables) )
    },
    
    # Overide the superclass
    generateConstraintRHS = function(){
      
      return(c(private$fixedTerminalConstraints$rhs,
               private$nodeDegreeConstraints$rhs,
               private$twoCycleConstraints$rhs,
               private$connectivityConstraints$rhs,
               private$novelSolutionsConstraint$rhs))
    },
    
    # Overide the superclass
    generateConstraintDirections = function(){
      
      return(c(private$fixedTerminalConstraints$directions,
               private$nodeDegreeConstraints$directions,
               private$twoCycleConstraints$directions,
               private$connectivityConstraints$directions,
               private$novelSolutionsConstraint$directions))
    },
    
    
    # Add a constraint that we cannot have a soluton that we have already seen
    # This constraint is not from the original paper, but it is quite simple
    # For each solution, sum_i y_i > 0 for i !in a solution
    setNoveltyConstraints = function(){
      
      noveltyConstraintsList <-
        private$solutionIndciesPool %>%
        as.list %>%
        lapply(function(solIndices){
          noveltyConstraint <- Matrix(1, nrow = 1, ncol = vcount(private$searchGraph), sparse = TRUE)
          noveltyConstraint[solIndices] <- 0
          return(noveltyConstraint)})
      
      if(private$verbosity) message("Adding ", length(noveltyConstraintsList)," novelty constraint(s) ...")
      
      #Deal with empty solution pools - add a matrix with no rows but the correct columns
      noveltyConstraintsList %<>% c(list(Matrix(nrow = 0, ncol = vcount(private$searchGraph))))
      
      noveltyConstraintsMatrix <- Reduce(rbind, noveltyConstraintsList)
      
      private$novelSolutionsConstraint <- list(
        variables = noveltyConstraintsMatrix,
        directions = rep(">=",nrow(noveltyConstraintsMatrix)) ,
        rhs = rep(1,nrow(noveltyConstraintsMatrix)))
    },
    
    solutionIndciesPool = sets::set(),
    
    novelSolutionsConstraint = list(),
    
    tolerance = numeric()
  )

)
