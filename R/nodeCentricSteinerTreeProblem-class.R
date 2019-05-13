
nodeCentricSteinerTreeProblem <- R6Class("nodeCentricSteinerTreeProblem",

                              public = list(

                                initialize = function(){

                                  # TODO parameter checking

                                  # search graph - undirected and simple. NAMED

                                  # edgeDT - both directions
                                  # nodeDT with node indicies

                                  # nodeScores - all -1 if absent

                                  private$addFixedTerminalConstraints()
                                  private$addNodeDegreeInequalities()
                                  private$addTwoCycleInequalities()
                                },

                                solve = function(){}

                                #Access to the constaint matricies will make for easier testing
                                getFixedTerminalConstraints <- function(){private$fixedTerminalConstraints},
                                getNodeDegreeConstraints <- function(){private$nodeDegreeConstraints},
                                getTwoCycleConstraints <- function(){private$twoCycleConstraints}
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
                                addNodeDegreeInequalities = function(){

                                  nodeDegreeInequalities_variables <- get.adjacency(private$searchGraph, sparse = TRUE)
                                  diag(nodeDegreeInequalities_variables) <- -2
                                  diag(nodeDegreeInequalities_variables)[ private$terminalIndicies ] <- -1

                                  private$nodeDegreeConstraints <- list(
                                    variables = nodeDegreeInequalities_variables,
                                    directions = rep(">=", nrow(nodeDegreeInequalities_variables)),
                                    rhs = rep(0, nrow(nodeDegreeInequalities_variables)) )

                                  invisible(self)
                                },

                                # Constraint 6.) If a potential node j sits adjacent to a potential terminal i, then i must be included if j is
                                addTwoCycleInequalities = function(){

                                  # Two cycle constraints are simple - for each edge i->j around a potential terminal i, have a setup which enforces its inclusion if j is present in solutoion
                                  twoCycle_variables <- private$edgeDT[fromNodeID %in% private$potentialTerminalIndicies,
                                                               sparseMatrix(i = c(.edgeID, .edgeID) ,
                                                                            j = c(fromNodeID, toNodeID),
                                                                            x = rep(c(1, -1) , each = length(.edgeID) ),
                                                                            dims = c( max(.edgeID), vcount(private$searchGraph)),
                                                                            dimnames = list( paste("twoCycleOnEdge", 1:max(.edgeID)) , V(private$searchGraph)$name))]

                                  # We only care about edges *from* a potential terminal node
                                  twoCycle_variables  %<>% .[ private$edgeDT[fromNodeID %in% private$potentialTerminalIndicies, .edgeID], ]

                                  private$twoCycleConstraints <- list(
                                    variables = twoCycle_variables,
                                    directions = rep(">=", nrow(twoCycle_variables)),
                                    rhs = rep(0, nrow(twoCycle_variables)) )

                                  invisible(self)
                                },

                                # Constraint 2.) Aims to endorce connections (via minimal node seperators) for pairs of terminals that are disconnected
                                # Applied lazily over terminal pairs in seperate components
                                addConnectivityConstraints = function(){}

                                searchGraph = graph.empty(),

                                # Work with integers rather than names
                                terminalIndicies = integer(),
                                fixedTerminalIndicies = integer(),
                                potentialTerminalIndicies = integer(),

                                nodeScores = numeric(),

                                fixedTerminalConstraints = list(),
                                nodeDegreeConstraints = list(),
                                twoCycleConstraints = list(),

                                edgeDT = data.table(),
                                nodeDT = data.table()
                              )
)
