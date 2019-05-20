

#' @importFrom sets set set_union
nodeCentricSteinerForestProblem <- R6Class("nodeCentricSteinerForestProblem",
                                           inherit = nodeCentricSteinerTreeProblem,
                                           public = list(

                                             #Overide
                                             initialize = function(network, solverChoice = chooseSolver(), verbose = TRUE, solverTimeLimit = 600){

                                               super$initialize(network,
                                                                solverChoice = solverChoice,
                                                                verbose = verbose,
                                                                presolveGraph = FALSE,
                                                                solverTimeLimit = solverTimeLimit)

                                               if(super$getNodeDT()[isTerminal == TRUE, nrow(.SD)] <= 2) stop("Steiner Forest routines require at least 3 fixed terminals (preferably many more!)")

                                               return(invisible(self))
                                             },

                                             sampleMultipleBootstrapSteinerSolutions = function(nBootstraps = 100){

                                               # solve normal steiner tree - this produces a bunch of connectivity constraints and
                                               # will also likely ensure that the solution is connected
                                               self$findSingleSteinerSolution()
                                               private$solutionIndciesPool <- set_union(self$getSolutionPool(), sets::set(private$currentSolutionIndices) )

                                               for(i in 1:nBootstraps){

                                                 private$resampleFixedTerminals()

                                                 super$solve()

                                                 #add solution graph if connected, else add connectivity constraints and resolve
                                                 if( super$isSolutionConnected() ){

                                                   # TODO keep a solution pool and test for a convergence in solutions
                                                   private$solutionIndciesPool <- set_union(self$getSolutionPool(), sets::set(private$currentSolutionIndices) )
                                                 }else{

                                                   super$addConnectivityConstraints()
                                                 }
                                               }

                                               return(invisible(self))
                                             },

                                             getSolutionPool = function(){

                                               return(private$solutionIndciesPool)
                                             },

                                             getSolutionPoolGraphs = function(collapseSols = TRUE){

                                               if(collapseSols){

                                                 #Ensure that the solution pool is up to date when we induce the subgraph. Since we are using a set, there is no cost to this
                                                 return( induced.subgraph(private$searchGraph, V(private$searchGraph)[unique(unlist( self$getSolutionPool()))]))
                                               }else{

                                                 return( self$getSolutionPool() %>%
                                                   as.list %>%
                                                   lapply( function(indices){ induced.subgraph(private$searchGraph, V(private$searchGraph)[indices])}) )
                                               }
                                             },

                                             #Overide: This overides the parent classes method and freshly regenerates the Steiner solution afresh each time
                                             findSingleSteinerSolution = function(maxItr = 20){

                                               private$fixedTerminalIndices <- super$getNodeDT()[isTerminal == TRUE, .nodeID]
                                               private$currentSolutionIndices <- integer()

                                               return(super$findSingleSteinerSolution(maxItr = maxItr))
                                             }

                                           ),
                                           private = list(

                                             resampleFixedTerminals = function(pSuccess = 0.5){

                                               if(private$verbosity){message("Bootstrap sampling seeds/fixed terminals ...")}

                                               #Flush the set of existing terminals
                                               private$fixedTerminalIndices <- integer()

                                               while(length(private$fixedTerminalIndices) == 0){

                                                 #Looks complicated, but really it just resamples the isTerminal/fixedTerminal nodeIDs
                                                 private$fixedTerminalIndices <- super$getNodeDT()[isTerminal == TRUE, .SD[sample(c(TRUE,FALSE), length(isTerminal), replace = TRUE, prob = c(pSuccess,1-pSuccess)), .nodeID]]
                                               }

                                               #Regenerate the fixed terminal constraints now that we have a different
                                               super$addFixedTerminalConstraints()

                                               return(invisible(self))
                                             },

                                             #This will be a set of integer sets
                                             solutionIndciesPool = sets::set()
                                           )
)
