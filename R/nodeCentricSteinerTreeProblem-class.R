
nodeCentricSteinerTreeProblem <- R6Class("nodeCentricSteinerTreeProblem",

                                         public = list(

                                           initialize = function(network, solverChoice = chooseSolver(), verbose = TRUE, presolveGraph = TRUE, solverTimeLimit = 600){

                                             private$solver <- validateSolverChoice(solverChoice)
                                             private$solverTimeLimit <- validateSingleInteger(solverTimeLimit)
                                             private$verbosity <- validateFlag(verbose)

                                             interactomeName <- deparse(substitute(network))

                                             validateIsNetwork(network)
                                             if(is.directed(network)){warning("Input network is directed and only undirected networks are supported - casting to a simple undirected network.")}


                                             private$graphPresolved <- validateFlag(presolveGraph)
                                             inputGraph <- network %>% as.undirected %>% simplify
                                             if(presolveGraph){inputGraph <- condenseSearchGraph(inputGraph)} #graph condensation is a presolve step

                                             private$searchGraph <-  set_graph_attr(inputGraph, "SearchNetwork", interactomeName)
                                             V(private$searchGraph)$.nodeID <- 1:vcount(private$searchGraph) #Assign a unique node integer to each node
                                             E(private$searchGraph)$.edgeID <- 1:ecount(private$searchGraph)

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

                                             private$fixedTerminalIndices <- private$nodeDT[isTerminal == TRUE, .nodeID] # Fixed terminals must be included in a solution
                                             private$potentialTerminalIndices <- private$nodeDT[nodeScore > 0, .nodeID] # potential terminals are those with nodeScore greater than 0

                                             # Check that there are *some* terminals, otherwise error
                                             if(length( unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices))) == 0) stop("No potential terminals (fixedTermals or potentialTerminals) presents. Review nodeScore and/or isTerminal vertex attributes!")

                                             # edgeDT - both directions for each arc
                                             private$edgeDT <- get.data.frame(as.directed(private$searchGraph, mode = "mutual"), what = "edges") %>% data.table %>% unique

                                             private$edgeDT[private$nodeDT, fromNodeID := .nodeID, on = .(from = name)]
                                             private$edgeDT[private$nodeDT, toNodeID := .nodeID, on = .(to = name)]

                                             # Set terminal constraints
                                             # What if no fixed terminals? Skip
                                             private$addFixedTerminalConstraints()

                                             private$addNodeDegreeInequalities()
                                             private$addTwoCycleInequalities()
                                             private$addConnectivityConstraints() #These will be empty at this stage

                                             return(invisible(self))
                                           },

                                           isSolutionConnected = function(){

                                             if(length(private$currentSolutionIndices) == 0){ return(FALSE) }

                                             return(is.connected( self$getCurrentSolutionGraph() ))
                                           },

                                           getCurrentSolutionGraph = function(){ return( induced.subgraph(private$searchGraph, V(private$searchGraph)[ private$currentSolutionIndices ]) ) },

                                           getTerminals = function(){

                                             return(list(fixedTerminals = private$fixedTerminalIndices,
                                                         potentialTerminals = private$potentialTerminalIndices,
                                                         terminals = unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices)) ))
                                           },

                                           findSingleSteinerSolution = function(maxItr = 20){

                                             itrCount <- 1

                                             while( (! self$isSolutionConnected() ) & (itrCount <= maxItr) ){

                                               private$solve()

                                               #Add connectivity constraints
                                               private$addConnectivityConstraints()

                                               itrCount %<>% add(1)
                                             }

                                             if(itrCount == maxItr) warning("Maximum number of solver iterations reached. In all likelihood the solution has not converged and may well be disconnected! Check!")

                                             return( uncondenseGraph( self$getCurrentSolutionGraph() ) ) #Uncondense graph undoes the graph presolve (or does nothing if the presolve step is omitted)
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
                                               sparseMatrix(i = private$fixedTerminalIndices,
                                                            j = private$fixedTerminalIndices,
                                                            x = 1,
                                                            dims = c( vcount(private$searchGraph), vcount(private$searchGraph)),
                                                            dimnames = list( paste0("fixedTerminalConstraintFor", 1:vcount(private$searchGraph)),
                                                                             V(private$searchGraph)$name))

                                             private$fixedTerminalConstraints <- list(
                                               variables = fixedTerminals_variables[private$fixedTerminalIndices, ],
                                               directions = rep("==", length(private$fixedTerminalIndices)),
                                               rhs = rep(1, length(private$fixedTerminalIndices)) )

                                             if(private$verbosity) message("Adding ", length(private$fixedTerminalIndices)," fixed terminal constraints ...")

                                             invisible(self)
                                           },

                                           # Constraint 5.) Only potential or fixed terminals can have a degree of 1, all other nodes must have degree >= 2
                                           addNodeDegreeInequalities = function(){

                                             nodeDegreeInequalities_variables <- get.adjacency(private$searchGraph, sparse = TRUE)
                                             diag(nodeDegreeInequalities_variables) <- -2
                                             diag(nodeDegreeInequalities_variables)[unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices))] <- -1

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
                                                                                               dimnames = list( paste0("twoCycleOnEdge", .edgeID) , V(private$searchGraph)$name))]

                                             # We only care about edges *from* a potential terminal node
                                             twoCycle_variables  %<>% .[ private$edgeDT[fromNodeID %in% private$potentialTerminalIndices, .edgeID], ]

                                             private$twoCycleConstraints <- list(
                                               variables = twoCycle_variables,
                                               directions = rep(">=", nrow(twoCycle_variables)),
                                               rhs = rep(0, nrow(twoCycle_variables)) )

                                             if(private$verbosity) message("Adding ",nrow(twoCycle_variables)," two-cycle constraints ...")

                                             invisible(self)
                                           },

                                           # Constraint 2.) Aims to endorce connections (via minimal node seperators) for pairs of terminals that are disconnected
                                           # This is the only really complicated constraint set in the function - study of it should be in conjunction with the paper Fischetti et al.
                                           # Only applied over terminal pairs in seperate components
                                           # TODO break out some of this functionality into smaller functions - it screams technical debt!
                                           addConnectivityConstraints = function(){

                                             # Using a dedicated column in private$nodeDT for cluster membership
                                             if(private$nodeDT[,all(is.na(inComponent))]){
                                               #i.e no nodes in any clusters - solver has to be called.

                                               private$connectivityConstraints <- list( variables = Matrix(0, sparse = TRUE, nrow = 0,
                                                                                                           ncol = vcount(private$searchGraph),
                                                                                                           dimnames = list(NULL, V(private$searchGraph)$name)),
                                                                                        directions = character(),
                                                                                        rhs = numeric() )

                                               if(private$verbosity) message("No connectivity constraints to add at this stage ...")

                                               return(invisible(self))
                                             }

                                             # Break solution into connected components. If only a single component, then we're done
                                             componentsSummaryDT <- private$nodeDT[!is.na(inComponent), .N, by = inComponent]
                                             terminalsInComponentsDT <-  private$nodeDT[!is.na(inComponent)][.nodeID %in% unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices))]

                                             if(private$verbosity){

                                               message(ifelse( self$isSolutionConnected() , "##Solution is connected!##", "##Solution is disconnected!##"))
                                               message("With ", componentsSummaryDT[,sum(N)]," nodes, ", nrow(terminalsInComponentsDT)," potential terminal(s) and ",
                                                       terminalsInComponentsDT[isTerminal == TRUE, nrow(.SD)], " fixed terminal(s) across ", nrow(componentsSummaryDT), " component(s).")
                                             }

                                             # Only add additional constraints if required
                                             if(! self$isSolutionConnected() ){

                                               # Create terminal pairs
                                               # Filter for terminal pairs in different components
                                               allTerminalPairs <- combn( unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices)), 2) %>% t %>% as.data.table
                                               setnames(allTerminalPairs, c("T1nodeID","T2nodeID")) # Note that we are working by node indicies here to avoid using node names

                                               allTerminalPairs[private$nodeDT, "T1inComponent" := inComponent, on = .(T1nodeID = .nodeID)]
                                               allTerminalPairs[private$nodeDT, "T2inComponent" := inComponent, on = .(T2nodeID = .nodeID)]

                                               # Only need keep node pairs where both nodes are in the solution but not in the same component
                                               terminalPairsInDifferentComponents <-
                                                 allTerminalPairs[!is.na(T1inComponent) & !is.na(T1inComponent)][T1inComponent != T2inComponent]

                                               # Pre-compute the component surfaces (A(C_i) in paper) - i.e. nodes that are adjacent to C_i but not in C_i
                                               clusterSurfacesDT <- private$nodeDT[!is.na(inComponent)] %>%
                                                 merge(private$edgeDT, by.x = ".nodeID", by.y = "fromNodeID") %>%
                                                 .[,.(componentSurfaceNodeID = .SD[!toNodeID %in% .nodeID, unique(toNodeID)]), by = inComponent]

                                               # Any edge that starts from a member of a cluster is stored in this table
                                               clusterEdgesDT <- private$nodeDT[!is.na(inComponent)] %>%
                                                 merge(private$edgeDT, by.x = ".nodeID", by.y = "fromNodeID")

                                               #A pot to keep the conConstraints in
                                               conConstraints <- list()

                                               # The computationally expensive piece
                                               for(T1component in unique(terminalPairsInDifferentComponents$T1inComponent) ){

                                                 termainPairsWithT1incomponent <- terminalPairsInDifferentComponents[T1inComponent == T1component]
                                                 termainPairsWithT1incomponent[, constraintIndex := .I]

                                                 clusterEdgeIDs <- clusterEdgesDT[inComponent == T1component, unique(.edgeID )]

                                                 #Crucially the nodes are still present from the original component
                                                 searchGraphWithoutComponent <- delete.edges(private$searchGraph,
                                                                                             E(private$searchGraph)[.edgeID %in% clusterEdgeIDs ] )

                                                 # For each terminal pair i,j; compute R_j (reachable set from R_j in graph ommiting C_i) and compute the minimum node-seperator set 𝒩(i,j) = A(C_i) intersect R_j

                                                 # R_j in the network wihtout C_i
                                                 # Note that this is the graph LESS the Ci componenet - hence we can't precompute upfront.
                                                 reachabilityMatrix <- distances(searchGraphWithoutComponent,
                                                                                 v = V(searchGraphWithoutComponent)[termainPairsWithT1incomponent$T2nodeID])
                                                 reachabilityMatrix[is.infinite(reachabilityMatrix)] <- 0

                                                 # A(C_i)
                                                 clusterSurfaceNodes <- clusterSurfacesDT[inComponent == T1component, componentSurfaceNodeID]
                                                 componentSurfaceMatrix <- sparseMatrix( i = rep(1:nrow(reachabilityMatrix), each = length(clusterSurfaceNodes) ) ,
                                                                                         j =  rep(clusterSurfaceNodes, nrow(reachabilityMatrix)),
                                                                                         x = TRUE,
                                                                                         dims = dim(reachabilityMatrix),
                                                                                         dimnames = dimnames(reachabilityMatrix))

                                                 # Add connectivity constraint
                                                 # y(N) ≥ y_i + y_j -1 ∀ i,j ∈ T, i ≠ j  ∀N ∈ 𝒩(i,j)
                                                 # i.e. if you have i then you must have a node seperator to have j included in the result
                                                 # This acts as a way of repairing a disconnected solution in a lazy fashion (there are exponentially many node seperators, hence enumerating them all up front is not feasible)

                                                 # 𝒩(i,j) = A(C_i) intersect R_j: y(N) above
                                                 minimumNodeSeperatorsMatrix <- ( (reachabilityMatrix) & componentSurfaceMatrix)

                                                 # These are the -1's: y_i + y_j
                                                 termianlPairValues <- termainPairsWithT1incomponent[, sparseMatrix( i =  rep(constraintIndex, times = 2) ,
                                                                                                                     j =  c(T1nodeID,T2nodeID) ,
                                                                                                                     x =  -1 ,
                                                                                                                     dims = dim(minimumNodeSeperatorsMatrix),
                                                                                                                     dimnames = dimnames(minimumNodeSeperatorsMatrix))]

                                                 minimumNodeSeperatorsMatrix <- minimumNodeSeperatorsMatrix + termianlPairValues

                                                 conConstraints %<>% c(minimumNodeSeperatorsMatrix)
                                               }

                                               allConnectivityConstraints <- Reduce(rbind, conConstraints)

                                               if(private$verbosity) message("Adding ",nrow(allConnectivityConstraints)," connectivity constraints based on node-seperators ...")

                                               # Append connectivity constraints matrix to existing variables, building up a pool of constraints that dictate connectivity
                                               private$connectivityConstraints <- list( variables = rbind(private$connectivityConstraints$variables, allConnectivityConstraints) )

                                               private$connectivityConstraints$directions <- rep(">=", nrow(private$connectivityConstraints$variables) )
                                               private$connectivityConstraints$rhs <- rep(-1, nrow(private$connectivityConstraints$variables) )
                                             }

                                             return(invisible(self))
                                           },

                                           # Optimise under current constraints using a solver agnostic interface
                                           solve = function(){

                                             if(private$verbosity) message("SOLVING ...")

                                             if(length(unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices)) ) == 1){

                                               if(private$verbosity) message("Only a single terminal (in the possibly presolved graph) - a trivial solution")
                                               private$currentSolutionIndices <- unique(c(private$fixedTerminalIndices, private$potentialTerminalIndices))
                                             }else{


                                               functionArgs <- list(

                                                 cVec = private$nodeDT[order(.nodeID), nodeScore],

                                                 Amat = rbind(private$fixedTerminalConstraints$variables,
                                                              private$nodeDegreeConstraints$variables,
                                                              private$twoCycleConstraints$variables,
                                                              private$connectivityConstraints$variables),

                                                 bVec = c(private$fixedTerminalConstraints$rhs,
                                                          private$nodeDegreeConstraints$rhs,
                                                          private$twoCycleConstraints$rhs,
                                                          private$connectivityConstraints$rhs),

                                                 vtypeVec = "B",

                                                 senseVec = c(private$fixedTerminalConstraints$directions,
                                                              private$nodeDegreeConstraints$directions,
                                                              private$twoCycleConstraints$directions,
                                                              private$connectivityConstraints$directions),

                                                 cplexParamList = list(trace = as.integer(private$verbosity),
                                                                       tilim = private$solverTimeLimit)
                                               )

                                               MILPsolve <- switch(private$solver ,
                                                                   RCPLEX = do.call("solver_CPLEX", functionArgs),
                                                                   CPLEXAPI = do.call("solver_CPLEXapi", functionArgs),
                                                                   LPSOLVE = do.call("solver_LPSOLVE", functionArgs),
                                                                   RGLPK = do.call("solver_GLPK", functionArgs),
                                                                   LPSYMPHONY = do.call("solver_SYMPHONY", functionArgs))

                                               solVec <- round(MILPsolve$solution)

                                               private$currentSolutionIndices <- which(solVec > 0)
                                             }

                                             disconnectedComponentList <- decompose( self$getCurrentSolutionGraph() )
                                            # The nodeDT table and inComponent variable keeps track of which node is where
                                             private$nodeDT[,inComponent := NA_integer_]

                                             for(i in 1:length(disconnectedComponentList)){

                                               # Using nodeIDs to track membership
                                               private$nodeDT[.nodeID %in% V(disconnectedComponentList[[i]])$.nodeID, inComponent := i]
                                             }

                                             return(invisible(self))
                                           },

                                           searchGraph = graph.empty(),

                                           currentSolutionIndices = integer(),

                                           # Work with integers rather than names
                                           fixedTerminalIndices = integer(),
                                           potentialTerminalIndices = integer(),

                                           fixedTerminalConstraints = list(),
                                           nodeDegreeConstraints = list(),
                                           twoCycleConstraints = list(),

                                           connectivityConstraints = list(),

                                           edgeDT = data.table(),
                                           nodeDT = data.table(),

                                           solver = character(),
                                           solverTimeLimit = integer(),

                                           graphPresolved = logical(),
                                           verbosity = logical()
                                         )
)
