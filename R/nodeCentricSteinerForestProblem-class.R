#' Construct an object representation of the bootstraped Steiner Tree process (aka Steiner Forest routine)
#' 
#' @description
#' Solve multiple bootstrap Minimum Steiner Tree problems (aka Steiner Forest procedure)
#'
#' @details
#' Sub-sample seed sets and gather robust solutions to seed subsets.
#' 
#' Given a set of seeds/fixed terminals a Minimum Steiner Tree can be found. One might well be interested in studying the common nodes that
#' would be included with, say, just 50% of the seed set. This process is known as a 'bootstrap' in statistics and this class looks to
#' repeatedly sample seeds to produce a consensus set of MStTP solutions. Sub-solutions can also be collected, albeit at an increased burden on the
#' solver (and therefore dramatically increasing the time).
#'
#' This class is derived from *subOptimalSteinerProblem* and in turn *nodeCentricSteinerTreeProblem*: all methods available in the superclass are available here. The difference is that after
#' each acceptable solution is found, the solution is a.) stored in a bootstrap solution pool and b.) used to generate a 'novelty' constraint on future solutions. For each bootstrap run, the solution pool
#' is flushed and the process re-rerun. In the end, all of the boostrap solutions are in the bootstrap solution pool.
#'
#' @examples
#' library(igraph)
#'
#' #Prepare a simple seed-based Steiner sampling in a reasonable sized network
#' fixedTerminalLymphomaGraph = lymphomaGraph
#' V(fixedTerminalLymphomaGraph)$isTerminal = FALSE
#' V(fixedTerminalLymphomaGraph)[nodeScore > 0]$isTerminal = TRUE
#' fixedTerminalLymphomaGraph = delete_vertex_attr(fixedTerminalLymphomaGraph, "nodeScore")
#'
#'
#' # Example of solving *just* the single-solution Minimum Steiner Tree Problem
#' MStTPsingle = nodeCentricSteinerTreeProblem$new(fixedTerminalLymphomaGraph)
#' MStTPsingle$findSingleSteinerSolution()
#'
#' #Solve multiple bootstrap Steiner Trees (Steiner Forest)
#' SteinFor = nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph)
#'
#' #Run two bootstrap routines (resample fixed terminals and solve) and
#' #ALSO run the sub-optimal solution searcher thrice
#' SteinFor$sampleMultipleBootstrapSteinerSolutions(nBootstraps = 2, maxItr = 3)
#' ## Takes around a minute using RGLPK as the solver
#'
#' @references Fischetti M, Leitner M, Ljubić I, Luipersbeck M, Monaci M, Resch M, et al. Thinning out Steiner trees: a node-based model for uniform edge costs. Math Program Comput. dimacs11.cs.princeton.edu; 2017;9: 203–229.
#' @references Beisser D, Klau GW, Dandekar T, Müller T, Dittrich MT. BioNet: An R-Package for the functional analysis of biological networks. Bioinformatics. 2010;26: 1129–1130.
#' @references \url{https://en.wikipedia.org/wiki/Steiner_tree_problem}
#' @family SteinerProblemSolver
#' @seealso nodeCentricSteinerTreeProblem
#' @seealso subOptimalSteinerProblem
#' @importFrom sets set_union
#' @importFrom purrr map
#' @export
nodeCentricSteinerForestProblem = R6Class(
  "nodeCentricSteinerForestProblem",
  inherit = subOptimalSteinerProblem,
  public = list(
    
    #' @description
    #' Constructor for the nodeCentricSteinerForestProblem class. Note the loss of 'presolveGraph'; the repeated resampling of fixed terminal nodes prevents this.
    #' @param network Search network, with either boolean isTerminal and/or continuous nodeScores recorded for each node in the search network
    #' @param solverChoice (optional) Select your preffered solver, or rely on default
    #' @param verbose Controls print verbosity of routine
    #' @param presolveGraph Whether to include the speed optimisation routine to coalesce adjacent nodes, decreasing the search space (default:TRUE - strongly recommended)
    #' @param solverTimeLimit Constrain how long, in seconds, each invocation of the MILP solver can take
    #' @param solverTrace Control how much detail to request from the solver  
    #' @param solutionTolerance If sub-optimal solutions are to be collected, what tolerance to use? (default: 0)
    #' @param RNGseed Integer seed to control ootstrap sampling process
    #' @return A `nodeCentricSteinerForestProblem` object, ready to collect solutions
    initialize = function(
      network,
      solverChoice = chooseSolver(),
      verbose = TRUE,
      solverTimeLimit = 300,
      solutionTolerance = 0,
      solverTrace = as.integer(verbose),
      RNGseed = sample.int(size = 1, .Machine$integer.max)
    ) {
      private$seedPool = RNGseed

      super$initialize(
        network,
        solverChoice = solverChoice,
        verbose = verbose,
        presolveGraph = FALSE,
        solutionTolerance = solutionTolerance,
        solverTrace = solverTrace,
        solverTimeLimit = solverTimeLimit
      )

      if (super$getNodeDT()[isTerminal == TRUE, nrow(.SD)] <= 2) {
        stop(
          "Steiner Forest routines require at least 3 fixed terminals (preferably many more!)"
        )
      }

      return(invisible(self))
    },

    #' @description
    #' Run the bootstrap procedure nBootstraps times, each time resampling seeds with pSuccess = resamplingProbability, collecting degenerate or suboptimal solutions for maxItr times.
    #' @param nBootstraps Number of bootstrap samplings
    #' @param maxItr The maximum number of constraint/solve cycles that the process should attempt
    #' @param resamplingProbability chance of node being chosen in a resampled subset
    #' @return self
    sampleMultipleBootstrapSteinerSolutions = function(
      nBootstraps = 5,
      maxItr = 0,
      resamplingProbability = 0.5
    ) {
      check_number_whole(nBootstraps, min = 1)
      check_number_whole(maxItr, min = 0)
      check_number_decimal(resamplingProbability, min = 0, max = 1)

      # solve normal steiner tree - this produces a bunch of connectivity constraints and
      # will also ensure that the solution is connected
      self$findSingleSteinerSolution()
      private$metasolutionIndicesPool = set_union(
        self$getBootstrapSolutionPool(),
        sets::set(private$currentSolutionIndices)
      )

      bootstrapItr = 1

      while (bootstrapItr <= nBootstraps) {
        if (private$verbosity) {
          message("Bootstrap ", bootstrapItr)
        }

        private$resampleFixedTerminals(resamplingProbability)

        #Find up to ten degenerate solutions as you can
        super$identifyMultipleSteinerSolutions(maxItr)

        private$nConnectivityConstraintsCallsPool = c(
          self$getNconnectivityConstraintsCallsPool(),
          super$getNconnectivityConstraintsCalls()
        )

        private$metasolutionIndicesPool = set_union(
          self$getBootstrapSolutionPool(),
          super$getSolutionPool()
        )

        #Flush the parent solution pool as we're about to research for solutions
        private$solutionIndicesPool = sets::set()

        bootstrapItr = bootstrapItr + 1
      }

      return(invisible(self))
    },

    #' @description
    #' Extract the set of node IDs appearing in solutions
    #' @return A set of integers that map to nodes
    getBootstrapSolutionPool = function() {
      return(private$metasolutionIndicesPool)
    },

    #' @description
    #' Getting for the 
    #' @return Matrix for use in MILP problem
    getNconnectivityConstraintsCallsPool = function() {
      return(private$nConnectivityConstraintsCallsPool)
    },

    #' @description
    #' Internal function aggegating solutions
    #' @return A set of integers that map to nodes
    getSolutionPool = function() {
      if (identical(parent.frame(), globalenv())) {
        warning(
          "During the Steiner forest process the solution pool constantly being flushed - it is likely to be empty. Use $getBootstrapSolutionPool() for Steiner forest problems."
        )
      }
      return(super$getSolutionPool())
    },


    #' @description
    #' Either return a list of solutions within tolerance (collapseSols = FALSE) or pool all solutions together and return a single graph (collapseSols = TRUE, defaults)
    #' @param collapseSols Return all graphs collapsed, or a list of graphs
    #' @return Either a list of graph solutions, or a collapsed version of all answers
    getBootstrapSolutionPoolGraphs = function(collapseSols = TRUE) {
      if (collapseSols) {
        #Ensure that the solution pool is up to date when we induce the subgraph. Since we are using a set, there is no cost to this
        return(induced.subgraph(
          private$searchGraph,
          V(
            private$searchGraph
          )[unique(unlist(self$getBootstrapSolutionPool()))]
        ))
      } else {
        return(
          self$getBootstrapSolutionPool() |>
            as.list() |>
            map(
              ~ induced.subgraph(
                private$searchGraph,
                V(private$searchGraph)[.x]
              )
            )
        )
      }
    },

    #' @description
    #' Internal function aggegating solutions
    #' @return A set of integers that map to nodes
    getSolutionPoolGraphs = function() {
      if (identical(parent.frame(), globalenv())) {
        warning(
          "During the Steiner forest process the solution pool constantly being flushed - it is likely to be empty. Use $getBootstrapSolutionPoolGraphs() for Steiner forest problems."
        )
      }
      return(super$getSolutionPoolGraphs())
    },

    #' @description
    #' overides the parent classes method and freshly regenerates the Steiner solution afresh each time. This is because we resample the seeds repeatedly.
    #' @param maxItr The maximum number of constraint/solve cycles that the process should attempt
    #' @return An induced subgraph of the best solution found given the constraints and interation cap
    findSingleSteinerSolution = function(maxItr = 20) {
      private$fixedTerminalIndices = super$getNodeDT()[
        isTerminal == TRUE,
        .nodeID
      ]
      private$currentSolutionIndices = integer()

      return(super$findSingleSteinerSolution(maxItr = maxItr))
    },

    #' @description
    #' Utility method to extract first seed for random sampling
    #' @return Integer seed
    getInitialSeed = function() {
      head(private$seedPool, n = 1)
    },

    #' @description
    #' Utility method to extract most recent seed for random sampling
    #' @return Integer seed
    getLatestSeed = function() {
      tail(private$seedPool, n = 1)
    },

    #' @description
    #' Utility method to return all seeds
    #' @return Integer seeds
    getAllSeeds = function() {
      return(private$seedPool)
    }
  ),
  private = list(
    addSeedToPool = function(seed) {
      private$seedPool = c(private$seedPool, seed)
      invisible(self)
    },

    sampleNewSeed = function() {
      return(sample.int(n = .Machine$integer.max, size = 1))
    },

    generateNextRNGseed = function() {
      latestSeed = self$getLatestSeed()
      set.seed(latestSeed)
      newSeed = private$sampleNewSeed()

      private$addSeedToPool(newSeed)

      invisible(self)
    },

    resampleFixedTerminals = function(pSuccess = 0.5) {
      if (private$verbosity) {
        message("Bootstrap sampling seeds/fixed terminals ...")
      }

      #Flush the set of existing terminals
      private$fixedTerminalIndices = integer()

      private$generateNextRNGseed()

      while (length(private$fixedTerminalIndices) == 0) {
        #Looks complicated, but really it just resamples the isTerminal/fixedTerminal nodeIDs
        private$fixedTerminalIndices = super$getNodeDT()[
          isTerminal == TRUE,
          .SD[runif(n = .N) <= pSuccess, .nodeID]
        ]

        if (any(duplicated(private$fixedTerminalIndices))) {
          warning(
            "Duplicated fixed terminals - this shouldn't be possible. Please contact package maintainer!: ",
            private$fixedTerminalIndices
          )
        }

        # I am unsure as to why there needs to be a unique call here - there should never be duplicated .nodeID values.
      }

      super$addFixedTerminalConstraints() #Regenerate the fixed terminal constraints now that we have a different fixed terminal set via bootstrap

      #The connectivity constraints relate to a different set of terminals (and therefore a different problem) - flush them
      if (private$verbosity) {
        message(
          "Flush existing connectivity constraints now that we have new fixed terminals"
        )
      }
      super$flushConnectivityConstraints()

      return(invisible(self))
    },

    metasolutionIndicesPool = sets::set(), #This will be a aggregated set of integer sets - the parent class has a solution pool - here we aggregate it!
    seedPool = integer(),
    nConnectivityConstraintsCallsPool = numeric()
  )
)
