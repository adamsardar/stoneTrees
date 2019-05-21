

#' @importFrom sets set set_union
nodeCentricSteinerForestProblem <- R6Class("nodeCentricSteinerForestProblem",
                                           inherit = subOptimalSteinerProblem,
                                           public = list(

                                             #Overide
                                             initialize = function(network, solverChoice = chooseSolver(), verbose = TRUE, solverTimeLimit = 300, solutionTolerance = 0, solverTrace = as.integer(verbose)){

                                               super$initialize(network,
                                                                solverChoice = solverChoice,
                                                                verbose = verbose,
                                                                presolveGraph = FALSE,
                                                                solutionTolerance = solutionTolerance,
                                                                solverTrace = solverTrace,
                                                                solverTimeLimit = solverTimeLimit)

                                               if(super$getNodeDT()[isTerminal == TRUE, nrow(.SD)] <= 2) stop("Steiner Forest routines require at least 3 fixed terminals (preferably many more!)")

                                               return(invisible(self))
                                             },

                                             sampleMultipleBootstrapSteinerSolutions = function(nBootstraps = 5, maxItr = 0, resamplingProbability= 0.5){

                                               validateSingleInteger(nBootstraps)
                                               validateSingleInteger(maxItr)
                                               validateSinglePositiveSemiDefiniteNumeric(resamplingProbability)

                                               # solve normal steiner tree - this produces a bunch of connectivity constraints and
                                               # will also ensure that the solution is connected
                                               self$findSingleSteinerSolution()
                                               private$metaSolutionIndciesPool <- set_union(self$getBootstrapSolutionPool(), sets::set(private$currentSolutionIndices))

                                               bootItr <- 1

                                               while(bootItr <= nBootstraps){

                                                 private$resampleFixedTerminals(resamplingProbability)

                                                 #Find up to ten degenerate solutions as you can
                                                 super$identifyMultipleSteinerSolutions(maxItr)

                                                 private$metaSolutionIndciesPool <- set_union(self$getBootstrapSolutionPool(), super$getSolutionPool())

                                                 #Flush the parent solution pool as we're about to research for solutions
                                                 private$solutionIndciesPool <- sets::set()

                                                 bootItr %<>% add(1)
                                               }

                                               return(invisible(self))
                                             },

                                             getBootstrapSolutionPool = function(){

                                               return(private$metaSolutionIndciesPool)
                                             },

                                             #Overide
                                             getSolutionPool = function(){

                                               if(identical(parent.frame(), globalenv())) warning("During the Steiner forest process the solution pool constantly being flushed - it is likely to be empty. Use $getBootstrapSolutionPool() for Steiner forest problems.")
                                               return(super$getSolutionPool())
                                             },

                                             getBootstrapSolutionPoolGraphs = function(collapseSols = TRUE){

                                               if(collapseSols){

                                                 #Ensure that the solution pool is up to date when we induce the subgraph. Since we are using a set, there is no cost to this
                                                 return( induced.subgraph(private$searchGraph, V(private$searchGraph)[unique(unlist( self$getBootstrapSolutionPool()))]))
                                               }else{

                                                 return( self$getBootstrapSolutionPool() %>%
                                                           as.list %>%
                                                           lapply( function(indices){ induced.subgraph(private$searchGraph, V(private$searchGraph)[indices])}) )
                                               }
                                             },

                                             #Overide
                                             getSolutionPoolGraphs = function(){

                                               if(identical(parent.frame(), globalenv())) warning("During the Steiner forest process the solution pool constantly being flushed - it is likely to be empty. Use $getBootstrapSolutionPoolGraphs() for Steiner forest problems.")
                                               return(super$getSolutionPoolGraphs())
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

                                             #This will be a set of integer sets - the parent class has a solution pool - here we aggrgate it!
                                             metaSolutionIndciesPool = sets::set()
                                           )
)
