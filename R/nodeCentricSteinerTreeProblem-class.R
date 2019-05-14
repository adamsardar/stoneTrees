
nodeCentricSteinerTreeProblem <- R6Class("nodeCentricSteinerTreeProblem",

                              public = list(

                                initialize = function(network, solverChoice = "GLPK", verbose = TRUE){

                                  # TODO validate solver (and that solver is available)
                                  private$solver <- solverChoice

                                  #TODO validate verbose
                                  private$verbosity <- verbose

                                  # TODO check network (particularly that it is named)
                                  if(is.directed(network)){warning("Input network is directed and only undirected networks are supported - casting to a simple undirected network.")}
                                  private$searchGraph <- network %>% as.undirected %>% simplify
                                  V(private$searchGraph)$.nodeID <- 1:vcount(private$searchGraph) #Assign a unique node integer to each node

                                  if(length(decompose(private$searchGraph)) != 1) stop("Search network must only have a single connected component.")

                                  # nodeDT with node indicies
                                  private$nodeDT <- get.data.frame(private$searchGraph, what = "vertices") %>% data.table

                                  # Solution status will effectively be provided by the 'inComponent' attribute
                                  private$nodeDT[, inComponent := NA_integer_]

                                  # check for isTerminal vertex attribute: all FALSE if absent, validate if present
                                  if(! "isTerminal" %in% colnames(private$nodeDT)){

                                    private$nodeDT[, isTerminal := FALSE]
                                  }else{

                                    if( ! (is.logical(private$nodeDT$isTerminal) & all(!is.na(private$nodeDT$isTerminal))) ) stop("isTerminal node attributes *must* all be boolean, with no NA's")
                                  }

                                  # check for nodeScore: all -1 if absent, validate if present
                                  if(! "nodeScore" %in% colnames(private$nodeDT)){

                                    private$nodeDT[, nodeScore := -1]
                                  }else{

                                    if( ! (is.numeric(private$nodeDT$nodeScore) & all(!is.na(private$nodeDT$nodeScore))) ) stop("nodeScore node attributes *must* all be numeric values, with no NA's")
                                  }


                                  private$fixedTerminalIndicies <- private$nodeDT[isTerminal == TRUE, .nodeID] # Fixed terminals must be included in a solution
                                  private$potentialTerminalIndicies <- private$nodeDT[nodeScore > 0, .nodeID] # potential terminals are those with nodeScore greater than 0

                                  private$terminalIndicies <- unique(c(private$fixedTerminalIndicies, private$potentialTerminalIndicies))

                                  # Check that there are *some* potential terminals, otherwise error
                                  if(length(private$terminalIndicies) == 0) stop("No potential terminals (fixedTermals or potentialTerminals) presents. Review nodeScore and/or isTerminal vertex attributes!")

                                  # edgeDT - both directions for each arc
                                  private$edgeDT <- get.data.frame(as.directed(private$searchGraph, mode = "mutual"), what = "edges") %>% data.table %>% unique
                                  private$edgeDT[, .edgeID := .I]

                                  private$edgeDT[private$nodeDT, fromNodeID := .nodeID, on = .(from = name)]
                                  private$edgeDT[private$nodeDT, toNodeID := .nodeID, on = .(to = name)]

                                  # Set terminal constraints
                                  # What if no fixed terminals? Skip
                                  private$addFixedTerminalConstraints()

                                  private$addNodeDegreeInequalities()
                                  private$addTwoCycleInequalities()

                                  # populate initial connectivity constraint matrix - initially empty as we only add these constraints 'lazily' (i.e. as needed to repair disconnected solutions)
                                  private$twoCycleConstraints <- list( variables = Matrix(0, sparse = TRUE, nrow = 0,
                                                                                          ncol = vcount(private$searchGraph),
                                                                                          dimnames = list(NULL, V(private$searchGraph)$name)),
                                                                        directions = character(),
                                                                        rhs = numeric() )

                                  invisible(self)
                                },

                                isSolutionConnected = function(){

                                  if(private$nodeDT[,all(is.na(inComponent))]){ return(FALSE) }

                                  #If all nodes in solution (not NA) are in the same comonent, then the solution is connected
                                  return(nrow(private$nodeDT[!is.na(inComponent), .N, by = inComponent]) == 1)
                                  },

                                getSingleSteinerSolution = function(maxItr = 20){

                                  itrCount <- 1

                                  while( (! self$isSolutionConnected() ) & (itrCount <= maxItr) ){

                                    private$solve()

                                    #Add connectivity constraints
                                    private$addConnectivityConstraints()

                                    itrCount %<>% add(1)
                                  }

                                  return(sol)
                                },

                                getNodeDT = function(){ private$nodeDT },
                                getEdgeDT = function(){ private$edgeDT },

                                #Access to the constaint matricies will make for easier testing
                                getFixedTerminalConstraints = function(){ private$fixedTerminalConstraints },
                                getNodeDegreeConstraints = function(){ private$nodeDegreeConstraints },
                                getTwoCycleConstraints = function(){ private$twoCycleConstraints },
                                getConnectivityConstraints = function(){ private$connectivityConstraints }
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
                                                          variables = fixedTerminals_variables[private$fixedTerminalIndicies, ],
                                                          directions = rep("==", length(private$fixedTerminalIndicies)),
                                                          rhs = rep(1, length(private$fixedTerminalIndicies)) )

                                  if(private$verbosity) message("Adding ", length(private$fixedTerminalIndicies)," fixed terminal constraints ...")

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

                                  if(private$verbosity) message("Adding ",nrow(nodeDegreeInequalities_variables)," node degree inequality constraints ...")

                                  invisible(self)
                                },

                                # Constraint 6.) If a potential node j sits adjacent to a potential terminal i, then i must be included if j is
                                addTwoCycleInequalities = function(){

                                  # 2-cycle inequalities (6)
                                  # Observe the following: if node i ∈ V is adjacent to a node j ∈ T_p, so that c_ij< p_j, then if i is part of the optimal solution, j has to be included as well
                                  # y_i ≤ y_j ∀ i ∈ V, j ∈ T_p

                                  # Two cycle constraints are simple - for each edge i->j around a potential terminal i, have a setup which enforces its inclusion if j is present in solutoion
                                  twoCycle_variables <- private$edgeDT[,
                                                               sparseMatrix(i = c(.edgeID, .edgeID) ,
                                                                            j = c(fromNodeID, toNodeID),
                                                                            x = rep(c(1, -1) , each = length(.edgeID) ),
                                                                            dims = c( length(.edgeID), vcount(private$searchGraph)),
                                                                            dimnames = list( paste("twoCycleOnEdge", .edgeID) , V(private$searchGraph)$name))]

                                  # We only care about edges *from* a potential terminal node
                                  twoCycle_variables  %<>% .[ private$edgeDT[fromNodeID %in% private$potentialTerminalIndicies, .edgeID], ]

                                  private$twoCycleConstraints <- list(
                                    variables = twoCycle_variables,
                                    directions = rep(">=", nrow(twoCycle_variables)),
                                    rhs = rep(0, nrow(twoCycle_variables)) )

                                  if(private$verbosity) message("Adding ",nrow(twoCycle_variables)," two-cycle constraints ...")

                                  invisible(self)
                                },

                                # Constraint 2.) Aims to endorce connections (via minimal node seperators) for pairs of terminals that are disconnected
                                # This is the only really complicated constraint set - study of it should be in conjunction with the paper Fischetti et al.
                                # Only applied over terminal pairs in seperate components
                                addConnectivityConstraints = function(){

                                  # Using a dedicated column in private$nodeDT for cluster membership
                                  if(private$nodeDT[,all(is.na(inComponent))]){ stop("No nodes in components for solution! Call solver first?") }

                                  # Break solution into connected components. If only a single component, then we're done
                                  componentsSummaryDT <- private$nodeDT[!is.na(inComponent), .N, by = inComponent]
                                  terminalsInComponentsDT <-  private$nodeDT[!is.na(inComponent) & isTerminal]

                                  solutionIsConnected <- nrow(componentsSummaryDT) == 1

                                  if(private$verbosity){

                                    message(ifelse(solutionIsConnected, "##Solution is connected!##", "##Solution is disconnected!##"))
                                    message("With ", componentsSummaryDT$N," nodes and ", nrow(terminalsInComponentsDT)," terminals.")
                                  }

                                  if(solutionIsConnected){ invisible(self) }

                                  # Create terminal pairs
                                  # Filter for terminal pairs in different components
                                  allTerminalPairs <- combn( names(private$terminalIndicies), 2) %>% t %>% as.data.table
                                  setnames(allTerminalPairs, c("T1nodeID","T2nodeID")) # Note that we are working by node indicies here to avoid using node names

                                  allTerminalPairs[private$nodeDT, c("T1inComponent","T1name") := list(inComponent, name), on = .(T1nodeID = .nodeID)]
                                  allTerminalPairs[private$nodeDT, c("T2inComponent","T2name") := list(inComponent, name), on = .(T2nodeID = .nodeID)]

                                  # Only need keep node pairs where both nodes are in the solution but not in the same component
                                  terminalPairsInDifferentComponents <-
                                       allTerminalPairs[!is.na(T1inComponent) & !is.na(T1inComponent)][T1inComponent != T2inComponent]

                                  # Pre-compute the component surfaces (A(C_i) in paper) - i.e. nodes that are adjacent to C_i but not in C_i
                                  clusterSurfacesDT <- private$nodeDT[!is.na(inComponent)] %>%
                                    merge(private$edgeDT, by.x = ".nodeID", by.y = "fromNodeID") %>%
                                    .[,.(componentSurfaceNodeID = .SD[!toNodeID %in% .nodeID, unique(toNodeID)]), by = inComponent]

                                  # Matrix to store connectivity constraints
                                  allConnectivityConstraints <- Matrix(0, sparse = TRUE,
                                                                       nrow = nrow(terminalPairsInDifferentComponents),
                                                                       ncol = vcount(private$searchGraph),
                                                                       dimnames = list(NULL, V(private$searchGraph)$name) )

                                  graphDiameter <- diameter(private$searchGraph)

                                  if(private$verbosity) message("Adding ",nrow(terminalPairsInDifferentComponents)," connectivity constraints based upon node seperators")

                                  # For each terminal pair i,j; compute R_j (reachable set from R_j in graph ommiting C_i) and compute the minimum node-seperator set N_ij = A(C_i) intersect R_j
                                  for(r in 1:nrow(terminalPairsInDifferentComponents)){

                                    terminalPair <- terminalPairsInDifferentComponents[r,]

                                    C_i_componentSurfaceNodeIDs <- clusterSurfacesDT[ inComponent ==  terminalPair$T1inComponent, componentSurfaceNodeID]

                                    graph_omitCi <- delete_vertices(fullSearchGraph,
                                                                    V(fullSearchGraph)[ private$nodeDT[inComponent ==  terminalPair$T1inComponent, .nodeID] ])

                                    # Note that this is the graph LESS the Ci componenet - hence we can't precompute upfront.
                                    # Notice the use of terminal name, not ID
                                    # Reffered to as Rj in the paper
                                    nodeNamesReachableFromJ <- V(make_ego_graph(graph_omitCi,
                                                                                nodes = V(graph_omitCi)[name == terminalPair$T2name],
                                                                                order = graphDiameter)[[1]])$name

                                    # Nodes adjacent to C_i and in R_j consititues a minimal node seperator 𝒩(i,j)
                                    # Identify minimal node seperator set 𝒩(i,j) = A(C_i) intersect R_j
                                    minNodeSep_ij_nodeIDs <- intersect(C_i_componentSurfaceNodeIDs,
                                                                       private$nodeDT[name %in% nodeNamesReachableFromJ, .nodeID])

                                    # Add connectivity constraint
                                    # y(N) ≥ y_i + y_j -1 ∀ i,j ∈ T, i ≠ j  ∀N ∈ 𝒩(i,j)
                                    # i.e. if you have i then you must have a node seperator to have j included in the result
                                    # This acts as a way of repairing a disconnected solution in a lazy fashion (there are exponentially many node seperators, hence enumerating them all up front is not feasible)
                                    allConnectivityConstraints[r, terminalPair$T1nodeID] <- -1
                                    allConnectivityConstraints[r, terminalPair$T2nodeID] <- -1

                                    allConnectivityConstraints[r, minNodeSep_ij_nodeIDs] <- 1
                                  }

                                  # Append connectivity constraints matrix to existing
                                  private$twoCycleConstraints <- list( variables = rbind(private$twoCycleConstraints$variables, allConnectivityConstraints) )

                                  private$twoCycleConstraints$directions <- rep(">=", nrow(private$twoCycleConstraints$variables) )
                                  private$twoCycleConstraints$rhs <- rep(-1, nrow(private$twoCycleConstraints$variables) )

                                  invisible(self)
                                },

                                # Optimise under current constraints using a solver agnostic interface
                                solve = function(){

                                  if(private$verbosity) message("SOLVING ...")

                                  if(length(private$potentialTerminalIndicies) == 1){

                                    message("Only a single potential termianl - a trivial solution")
                                    solutionIndicies <- private$potentialTerminalIndicies
                                  }else{

                                    GLPKsolution <- Rglpk_solve_LP(

                                      obj = private$nodeDT[order(.nodeID), nodeScore],

                                      mat = rbind(private$fixedTerminalConstraints$variables,
                                                  private$nodeDegreeConstraints$variables,
                                                  private$twoCycleConstraints$variables,
                                                  private$connectivityConstraints$variables),

                                      dir = c(private$fixedTerminalConstraints$directions,
                                              private$nodeDegreeConstraints$directions,
                                              private$twoCycleConstraints$directions,
                                              private$connectivityConstraints$directions),

                                      rhs = c(private$fixedTerminalConstraints$rhs,
                                              private$nodeDegreeConstraints$rhs,
                                              private$twoCycleConstraints$rhs,
                                              private$connectivityConstraints$rhs),

                                      max = TRUE,

                                      control = list(verbose = private$verbosity),

                                      types = "B")

                                    solutionIndicies <- which(GLPKsolution$solution > 0)
                                  }


                                  solIndex <<- solutionIndicies


                                  #
                                  graphOfSolution <- induced_subgraph(private$searchGraph,
                                                                      V(private$searchGraph)[solutionIndicies])

                                  disconnectedComponentList <- decompose(graphOfSolution)

                                  private$nodeDT[,inComponent := NA_integer_]

                                  for(i in 1:length(disconnectedComponentList)){

                                    private$nodeDT[.nodeID %in% V(disconnectedComponentList[[i]])$.nodeID, inComponent := i]
                                  }

                                  invisible()
                                },

                                searchGraph = graph.empty(),

                                # Work with integers rather than names
                                terminalIndicies = integer(),
                                fixedTerminalIndicies = integer(),
                                potentialTerminalIndicies = integer(),

                                fixedTerminalConstraints = list(),
                                nodeDegreeConstraints = list(),
                                twoCycleConstraints = list(),

                                connectivityConstraints = list(),

                                edgeDT = data.table(),
                                nodeDT = data.table(),

                                solver = character(),
                                verbosity = logical()
                              )

)
