
#' @importFrom sets set set_union
subOptimalSteinerProblem <- R6Class("subOptimalSteinerProblem",
                                           inherit = nodeCentricSteinerTreeProblem,
                                           public = list(


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

                                           ),

                                           private = list(

                                             # Overide the superclass
                                             generateConstraintMatrix = function(){

                                               return(rbind(super$fixedTerminalConstraints$variables,
                                                            super$nodeDegreeConstraints$variables,
                                                            super$twoCycleConstraints$variables,
                                                            super$connectivityConstraints$variables,
                                                            private$novelSolutionsConstraint$variables) )
                                             },

                                             # Overide the superclass
                                             generateConstraintRHS = function(){

                                               return(c(super$fixedTerminalConstraints$rhs,
                                                        super$nodeDegreeConstraints$rhs,
                                                        super$twoCycleConstraints$rhs,
                                                        super$connectivityConstraints$rhs,
                                                        private$novelSolutionsConstraint$rhs))
                                             },

                                             # Overide the superclass
                                             generateConstraintDirections = function(){

                                               return(c(super$fixedTerminalConstraints$directions,
                                                        super$nodeDegreeConstraints$directions,
                                                        super$twoCycleConstraints$directions,
                                                        super$connectivityConstraints$directions,
                                                        private$novelSolutionsConstraint$directions))
                                             },


                                             # Add a constraint that we cannot have a soluton that we have already seen
                                             # This constraint is not from the original paper, but it is quite simple
                                             # For each solution, sum_i y_i > 0 for i !in a solution
                                             setNoveltyConstraints = function(){

                                               if(length(private$solutionPool) == 0){

                                                 return(Matrix(nrow = 0,
                                                               ncol = vcount(super$searchGraph),
                                                               dimnames = list(NULL, V(super$searchGraph)$name)) ) }


                                               onesMat <- Matrix(1, nrow = 1, ncol = vcount(super$searchGraph))

                                               noveltyConstraintsList <- private$solutionPool %>%
                                                                     as.list %>%
                                                 lapply(function(solIndices){
                                                   noveltyConstraint <- onesMat
                                                   noveltyConstraint[solIndices] <- 0
                                                   return(noveltyConstraint)})

                                               private$novelSolutionsConstraint <- list(
                                                  variables = Reduce(rbind, noveltyConstraintsList),
                                                  direction = ">",
                                                  rhs = 0)
                                             },

                                             solutionPool = sets::set(),
                                             novelSolutionsConstraint = list(),
                                           )

)
