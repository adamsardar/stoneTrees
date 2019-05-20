context("Inspect steiner forest routines")

library(igraph)
library(sets)

#Prepare a simple seed-based Steiner sampling in a reasonable sized network
fixedTerminalLymphomaGraph <- lymphomaGraph
V(fixedTerminalLymphomaGraph)$isTerminal <- FALSE
V(fixedTerminalLymphomaGraph)[nodeScore > 0]$isTerminal <- TRUE
V(fixedTerminalLymphomaGraph)$nodeScore <- -1

test_that("Inspect object construction",{

  expect_error(nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, presolveGraph = TRUE),
               regexp = "presolveGraph", label = "Graphs can't be presolved for Steiner forest")

})

test_that("Run a small Steiner forest routine and check the results", {

  steinFor <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE)

  expect_false(steinFor$.__enclos_env__$private$graphPresolved, label = "No presolve on the graph when permuting constraints")
  expect_true(set_is_empty(steinFor$getSolutionPool()), label = "Before calling the solver, solution pool must be empty")

  simpleTreeSol <- steinFor$findSingleSteinerSolution()

  expect_silent( steinFor$sampleMultipleBootstrapSteinerSolutions(10))
  expect_false( set_is_empty(steinFor$getSolutionPool()), label = "Once populated, pool must be empty")

  expect_gte(vcount(steinFor$getSolutionPoolGraphs(collapseSols = TRUE)), vcount(simpleTreeSol), label = "Forest must be larger than tree by definition")

  expect_true( is.igraph(steinFor$getSolutionPoolGraphs(collapseSols = TRUE)), label = "Collapsed solution should be an igraph")

  expect_equal( length(steinFor$getSolutionPoolGraphs(collapseSols = FALSE)), length(steinFor$getSolutionPool()), label = "Solution pool must be equal to graph list in size")
  expect_true( all(sapply(steinFor$getSolutionPoolGraphs(collapseSols = FALSE), is.igraph)), label = "All solutions should be igraphs")

  expect_true(is.connected(steinFor$getSolutionPoolGraphs(collapseSols = TRUE)), label = "Solution must be connected")
  expect_true( all( sapply(steinFor$getSolutionPoolGraphs(collapseSols = FALSE),is.connected) ), label = "Solution must be connected")
})

test_that("Checking ability to find many k-steiner trees in the karate graph",{

  karateGraph_MSTP_kStein <- nodeCentricSteinerForestProblem$new(karateGraph, verbose = FALSE)$sampleMultipleBootstrapSteinerSolutions(10)

  expect_true( is.igraph(karateGraph_MSTP_kStein$getSolutionPoolGraphs(collapseSols = TRUE)), label = "Collapsed solution should be an igraph")

  expect_equal( length(karateGraph_MSTP_kStein$getSolutionPoolGraphs(collapseSols = FALSE)), length(karateGraph_MSTP_kStein$getSolutionPool()), label = "Solution pool must be equal to graph list in size")
  expect_true( all(sapply(karateGraph_MSTP_kStein$getSolutionPoolGraphs(collapseSols = FALSE), is.igraph)), label = "All solutions should be igraphs")

  expect_true(is.connected(karateGraph_MSTP_kStein$getSolutionPoolGraphs(collapseSols = TRUE)), label = "Solution must be connected")
  expect_true( all( sapply(karateGraph_MSTP_kStein$getSolutionPoolGraphs(collapseSols = FALSE),is.connected) ), label = "Solution must be connected")
})
