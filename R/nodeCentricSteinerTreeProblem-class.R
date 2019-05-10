
nodeCentricSteinerTreeProblem <- R6Class("nodeCentricSteinerTreeProblem",

                              public = list(




                                initialize = function(){

                                  #TODO parameter checking


                                  private$addFixedTerminalConstraints()

                                },



                                solve = function(){}

                                getFixedTerminalConstraints <- function(){private$fixedTerminalConstraints}


                              ),

                              private = list(

                                # Constraint 3.) Nodes designated as fixed terminals MUST be in the solutions
                                addFixedTerminalConstraints = function(){

                                  fixedTerminals_variables <-
                                    sparseMatrix(i = private$fixedTerminalsIndicies,
                                                 j = private$fixedTerminalsIndicies,
                                                 x = 1,
                                                 dims = c( vcount(private$searchGraph), vcount(private$searchGraph)),
                                                 dimnames = list( paste0("fixedTerminalConstraintFor", 1:vcount(private$searchGraph)),
                                                                  V(private$searchGraph)$name))

                                  private$fixedTerminalConstraints <- list(
                                                          variables = fixedTerminals_variables[private$fixedTerminalsIndicies, ],
                                                          directions = rep("==", nrow(private$fixedTerminalsIndicies)),
                                                          rhs = rep(1, nrow(private$fixedTerminalsIndicies)) )
                                  invisible(self)
                                },

                                # Constraint 5.) Only potential or fixed terminals can have a degree of 1, all other nodes must have degree >= 2
                                addNodeDegreeInequalities = function(){},

                                # Constraint 6.) If a potential node i sits adjacent to a potential terminal j, then j must be included if i is
                                addTwoCycleInequalities = function(){},

                                # Constraint 2.) Aims to endorce connections (via minimal node seperators) for pairs of terminals that are disconnected
                                # Applied lazily over terminal pairs in seperate components
                                addConnectivityConstraints = function(){}

                                searchGraph = graph.empty(),

                                # Work with integers rather than names
                                terminalIndicies = integer(),
                                fixedTerminalIndicies = integer(),
                                potentialTerminalIndicies = integer(),

                                nodeScores = numeric(),

                                fixedTerminalConstraints = list()
                              )
)
