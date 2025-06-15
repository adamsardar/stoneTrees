#' Construct an object representation of a multiple-solution Steiner tree/maximum weight connected subgraph (MWCS) problem
#' 
#' @description
#' Collect degenerate and sub-optimal solutions to Steiner problems (MStTP or MWCS) with uniform or no edge weights.
#'
#' @details 
#' Rather than find just a single solution to a MStTP/MWCS, one can populate a solution pool with multiple degenerate/tolerable solutions.
#'
#' This class is derived from *nodeCentricSteinerTreeProblem*: all methods available in the superclass are available here. The difference is that after
#' each acceptable solution is found, the solution is a.) stored in a solution pool and b.) used to generate a 'novelty' constraint on future solutions.
#' 
#' @examples
#' library(igraph)
#'
#'  # Maximum-Weight Connected Subgraph (MWCS) - find sub-optimal solutions
#'
#'  ## Vertex attribute details node costs/prizes
#'  head(V(lymphomaGraph)$nodeScore)
#'
#'  lymphoma_multiMWCS = subOptimalSteinerProblem$new(lymphomaGraph, solutionTolerance = 0.5)
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
#' @importFrom sets set_union
#' @export
subOptimalSteinerProblem = R6Class(
  "subOptimalSteinerProblem",
  inherit = nodeCentricSteinerTreeProblem,
  public = list(
    #' @description
    #' Constructor for the subOptimalSteinerProblem class. Alongside the arguments for the super-class constructor, there is also 'solutionTolerance', which instructs the object as to the gap between optimal and observed solution that is acceptable.
    #' @param network Search network, with either boolean isTerminal and/or continuous nodeScores recorded for each node in the search network
    #' @param solverChoice (optional) Select your preffered solver, or rely on default
    #' @param verbose Controls print verbosity of routine
    #' @param presolveGraph Whether to include the speed optimisation routine to coalesce adjacent nodes, decreasing the search space (default:TRUE - strongly recommended)
    #' @param solverTimeLimit Constrain how long, in seconds, each invocation of the MILP solver can take
    #' @param solverTrace Control how much detail to request from the solver  
    #' @param solutionTolerance If sub-optimal solutions are to be collected, what tolerance to use? (default: 0)
    #' @return A `nodeCentricSteinerForestProblem` object, ready to collect solutions
    initialize = function(
      network,
      solverChoice = chooseSolver(),
      verbose = TRUE,
      presolveGraph = TRUE,
      solverTimeLimit = 300,
      solutionTolerance = 0,
      solverTrace = as.integer(verbose)
    ) {
      super$initialize(
        network = network,
        solverChoice = solverChoice,
        verbose = verbose,
        presolveGraph = presolveGraph,
        solverTrace = solverTrace,
        solverTimeLimit = solverTimeLimit
      )

      self$setSolutionTolerance(solutionTolerance + 1E-10) # Add epsilon to handle small fluctuations

      private$setNoveltyConstraints()

      return(invisible(self))
    },

    #' @description
    #' Internal function aggegating solutions
    #' @return A set of integers that map to nodes
    getSolutionPool = function() {
      return(private$solutionIndicesPool)
    },

    #' @description
    #' Internal function aggegating solutions
    #' @param collapseSols  Return all graphs collapsed, or a list of graphs
    #' @return A set of integers that map to nodes
    getSolutionPoolGraphs = function(collapseSols = TRUE) {
      if (collapseSols) {
        #Ensure that the solution pool is up to date when we induce the subgraph. Since we are using a set, there is no cost to this
        return(uncondenseGraph(induced.subgraph(
          private$searchGraph,
          V(private$searchGraph)[unique(unlist(self$getSolutionPool()))]
        )))
      } else {
        return(
          self$getSolutionPool() %>%
            as.list %>%
            lapply(function(indices) {
              induced.subgraph(
                private$searchGraph,
                V(private$searchGraph)[indices]
              )
            })
        ) %>%
          lapply(uncondenseGraph)
      }
    },

    #' @description
    #' Compute the scores of the solutions in the solution pool. These are in the same order as the list of graphs returned by $getSolutionPoolGraphs(FALSE)
    #' @return A set of continuous scores of solutions
    getSolutionPoolScores = function() {
      return(
        self$getSolutionPool() %>%
          as.list %>%
          sapply(function(indices) {
            super$getNodeDT()[.nodeID %in% indices, sum(nodeScore)]
          }) %>%
          unlist
      )
    },


    #' @description 
    #' Returns the optimum score from solutions in the solution pool
    #' @return Continuous value of the top score
    getOptimumScore = function() {
      return(max(self$getSolutionPoolScores(), na.rm = TRUE))
    },

    #' @description 
    #' Internal utility function that forces new solutions to be found
    #' @return Matrix of existing solutions
    getNoveltyConstraints = function() {
      return(private$novelSolutionsConstraint)
    },

    #' @description 
    #' setter for the solver solution tolerance
    #' @return Continuous value of the top score
    getSolutionTolerance = function() {
      return(private$tolerance)
    },

    #' @description 
    #' setter for the solver solution tolerance
    #' @param x Single continuous number to use for search 
    #' @return self
    setSolutionTolerance = function(x) {
      check_number_decimal(x, min = 0, allow_infinite = FALSE)
      private$tolerance = x
      return(invisible(self))
    },

    
    #' @description
    #' Provide a list of all connectivity constraint calls
    #' @return A list of constraints
    getNconnectivityConstraintsCalls = function() {
      private$nConnectivityConstraintsCalls
    },


    
    #' @description
    #' The multi-solution version of `$findSingleSteinerSolution`
    #' @param maxItr As we sample the sub-optimal solutions, what is the maximum number of within-tolerance solutions to aggregate?
    #' @return self
    identifyMultipleSteinerSolutions = function(maxItr = 10) {
      check_number_whole(maxItr, min = 0)

      self$findSingleSteinerSolution()

      private$nConnectivityConstraintsCalls = self$getNconnectivityConstraintsCalls()

      private$solutionIndicesPool = set_union(
        self$getSolutionPool(),
        sets::set(private$currentSolutionIndices)
      )

      multiSteinerItr = 1

      super$nConnectivityConstraintsCalls = 0

      while (multiSteinerItr <= maxItr) {
        private$setNoveltyConstraints()

        super$solve()
        multiSteinerItr %<>% add(1)

        if (vcount(super$getCurrentSolutionGraph()) == 0) {
          message(
            "STOP iteration, solution not found. No more novelty constraint added"
          )

          break()
        } else {
          #add solution graph if connected, else add connectivity constraints and resolve
          if (super$isSolutionConnected()) {
            #If the absolute difference between scores is within tolerance, then add to pool
            if (
              abs(super$getCurrentSolutionScore() - self$getOptimumScore()) <=
                private$tolerance
            ) {
              private$solutionIndicesPool = set_union(
                self$getSolutionPool(),
                sets::set(private$currentSolutionIndices)
              )

              private$nConnectivityConstraintsCalls = c(
                self$getNconnectivityConstraintsCalls(),
                super$nConnectivityConstraintsCalls
              )

              super$nConnectivityConstraintsCalls = 0
            } else {
              message(
                "Next feasible solution is outside of solution tolerance! Consider increasing it with $setSolutionTolerance(x) method?"
              )
              break()
            }
          } else {
            super$addConnectivityConstraints()

            super$nConnectivityConstraintsCalls = super$nConnectivityConstraintsCalls %<>%
              add(1)
          }
        }
      }

      return(invisible(self))
    }
  ),

  private = list(
    # Overide the superclass
    gatherConstraintObjects = function() {
      return(list(
        private$fixedTerminalConstraints,
        private$nodeDegreeConstraints,
        private$twoCycleConstraints,
        private$connectivityConstraints,
        private$novelSolutionsConstraint
      ))
    },

    # Add a constraint that we cannot have a solution that we have already seen
    # This constraint is not from the original paper, but it is quite simple
    # For each solution, sum_i y_i > 0 for i !in a solution
    setNoveltyConstraints = function() {
      noveltyConstraintsList = private$solutionIndicesPool %>%
        as.list %>%
        lapply(function(solIndices) {
          noveltyConstraint = Matrix(
            1,
            nrow = 1,
            ncol = vcount(private$searchGraph),
            sparse = TRUE
          )
          noveltyConstraint[solIndices] = 0
          return(noveltyConstraint)
        })

      if (private$verbosity) {
        message(
          "Adding ",
          length(noveltyConstraintsList),
          " novelty constraint(s) ..."
        )
      }

      #Deal with empty solution pools - add a matrix with no rows but the correct columns
      noveltyConstraintsList %<>%
        c(list(Matrix(nrow = 0, ncol = vcount(private$searchGraph))))

      noveltyConstraintsMatrix = Reduce(rbind, noveltyConstraintsList)

      private$novelSolutionsConstraint = list(
        variables = noveltyConstraintsMatrix,
        directions = rep(">=", nrow(noveltyConstraintsMatrix)),
        rhs = rep(1, nrow(noveltyConstraintsMatrix))
      )
    },

    solutionIndicesPool = sets::set(),

    novelSolutionsConstraint = list(),

    tolerance = numeric()
  )
)
