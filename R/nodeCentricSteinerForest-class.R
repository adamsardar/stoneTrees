

nodeCentricSteinerForestProblem <- R6Class("nodeCentricSteinerForestProblem",
                                           inherit = nodeCentricSteinerTreeProblem,
                                           public = list(


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
                                               private$solutionPool %<>% c( list(private$solutionGraph) )

                                               for(i in 1:nBootstraps){

                                                 private$resampleFixedTerminals()
                                                 private$solutionPool %<>% c( list(private$solutionGraph) )

                                                 super$solve()

                                                 #add solution graph if connected, else add connectivity constraints and resolve

                                                 if(is.connected(private$solutionGraph)){

                                                   # TODO keep a solution pool and test for a convergence in solutions
                                                   private$solutionPool %<>% c( list(private$solutionGraph) )
                                                 }else{

                                                   private$addConnectivityConstraints()
                                                 }
                                               }


                                               return(invisible(self))
                                             },

                                             getSolutionPool = function(){

                                               return(private$solutionPool)
                                             },

                                             #TODO have a dedicated addToSolutionPool method which checks previous solutions by hash

                                             #This overides the parent classes method and freshly regenerates the Steiner solution afresh each time
                                             findSingleSteinerSolution = function(maxItr = 20){

                                               private$fixedTerminalIndices <- super$getNodeDT()[isTerminal == TRUE, .nodeID]
                                               private$solutionGraph <- graph.empty()

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


                                             solutionPool = list()

                                           )
)
