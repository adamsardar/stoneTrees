
#' @importFrom sets set set_union
subOptimalSteinerProblem <- R6Class("subOptimalSteinerProblem",

                                           inherit = nodeCentricSteinerTreeProblem,

                                           public = list(

                                             initialize = function(network, solverChoice = chooseSolver(), verbose = TRUE, presolveGraph = TRUE, solverTimeLimit = 300, solutionTolerance = 0){

                                                super$initialize(network = network,
                                                          solverChoice = solverChoice,
                                                          verbose = verbose,
                                                          presolveGraph = presolveGraph,
                                                          solverTimeLimit = solverTimeLimit)

                                               self$setSolutionTolerance(solutionTolerance + 1E-8)  # Add epsilon

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

                                               for(i in 1:maxItr){

                                                 private$setNoveltyConstraints()

                                                 super$solve()

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
