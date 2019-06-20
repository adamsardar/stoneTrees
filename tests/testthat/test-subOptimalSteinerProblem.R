context("Inspecting node-centric steiner tree object creation")

library(igraph)
library(data.table)
library(sets)

test_that("Inspect sub-optimal solution searcher construction and answers found for MWCS",{

  expect_error(subOptimalSteinerProblem$new(lymphomaGraph, solutionTolerance = -1, verbose = FALSE), regexp = "positive", label = "negative tolerences are not allowed")

  expect_silent(testLymphoma <- subOptimalSteinerProblem$new(lymphomaGraph, solutionTolerance = 0.5, verbose = FALSE))

  expect_equal(testLymphoma$getSolutionPoolGraphs(collapseSols = FALSE), list(), label = "Before calling identifyMultipleSteinerSolutions() there should be no solutions in pool")
  expect_true( is.null(testLymphoma$getSolutionPoolScores()), label = "Before populating the solution pool, the solution scores must be null")

  singleSolGraph <- testLymphoma$findSingleSteinerSolution()

  expect_equal(testLymphoma$getSolutionPoolGraphs(collapseSols = FALSE), list(), label = "Before calling identifyMultipleSteinerSolutions() there should be no solutions in pool, even if findSingleSteinerSolution() has been called")

  expect_message(testLymphoma$identifyMultipleSteinerSolutions(), "solution tolerance")

  expect_gt(vcount(testLymphoma$getSolutionPoolGraphs(collapseSols = TRUE)), vcount(singleSolGraph), label = "More solutions must create a larger graph")

  expect_true(all(sapply(testLymphoma$getSolutionPoolGraphs(collapseSols = FALSE), is.igraph)))
  expect_true(all(sapply(testLymphoma$getSolutionPoolGraphs(collapseSols = FALSE), is.connected)))

  expect_true( is.igraph(testLymphoma$getSolutionPoolGraphs(collapseSols = TRUE)) )
  expect_true( is.connected(testLymphoma$getSolutionPoolGraphs(collapseSols = TRUE)) )

  expect_lte( diff(range(testLymphoma$getSolutionPoolScores())), testLymphoma$getSolutionTolerance(), label = "Range in solution scores must be within tolerence")

  solutionOutsideTol <- testLymphoma$findSingleSteinerSolution()

  expect_gt( diff(range( c(testLymphoma$getSolutionPoolScores(), testLymphoma$getCurrentSolutionScore())) ), testLymphoma$getSolutionTolerance(), label = "The next solution should be outside solution tolerence")

  expect_error(testLymphoma$setSolutionTolerance( "1" ), regexp = "positive")

  testLymphoma$setSolutionTolerance(1)

  testLymphoma$identifyMultipleSteinerSolutions()

  expect_lte( diff(range(testLymphoma$getSolutionPoolScores())), testLymphoma$getSolutionTolerance(), label = "Range in solution scores must be within tolerence")
})


test_that("Inspect sub-optimal solution searcher construction and answers found for MStP (without nodeScores)",{

<<<<<<< HEAD
  expect_silent(testKarate <- subOptimalSteinerProblem$new(karateGraph, solutionTolerance = 0, verbose = TRUE))

  expect_message(testKarate$identifyMultipleSteinerSolutions(), "solution tolerance")
=======
  expect_silent(testKarate <- subOptimalSteinerProblem$new(karateGraph, solutionTolerance = 0, verbose = FALSE))
>>>>>>> 7a0d6b3... Fix broken test logic

  expect_equal(vcount(testKarate$getSolutionPoolGraphs()), 5, label = "visual inspection of the graph shows that there are two 4 node solutions")

  expect_equal( nrow(testKarate$getNoveltyConstraints()$variables), length(testKarate$getSolutionPoolScores()), label = "There should be an equal number of novelty constraints as solutions in pool")

  expect_equal( nrow(testKarate$getNoveltyConstraints()$variables), length(testKarate$getNoveltyConstraints()$directions), label = "There should be an equal number of novelty constraints as directions")
  expect_equal( nrow(testKarate$getNoveltyConstraints()$variables), length(testKarate$getNoveltyConstraints()$rhs), label = "There should be an equal number of novelty constraints as rhs")

  expect_equal( vcount(testKarate$findSingleSteinerSolution()), 5)

  expect_equal( diff(range(testKarate$getSolutionPoolScores())), 0, label = "All solutions should be degenerate as solutionTolerance was set to 0")
})

