context("Compare results of the default solver, lpSolve, with all other supported solvers available")

library(igraph)

steinObject = nodeCentricSteinerTreeProblem$new(lymphomaGraph, verbose = FALSE, presolve = FALSE, solverChoice = "LPSOLVE")
lpSolveSolution = steinObject$findSingleSteinerSolution()

test_that("Check that lpSolve is producing to correct values (lymphoma graph is a DIMACS problem with a known solution)", {

  expect_equal(vcount(lpSolveSolution),46)
  expect_gt(sum(V(lpSolveSolution)$nodeScore),70)
  expect_true(is.connected(lpSolveSolution))
})


test_that("Compare cplexAPI result", {

  if(!"cplexAPI" %in% .packages(all.available = TRUE)) skip("cplexAPI not installed")

    steinObject = nodeCentricSteinerTreeProblem$new(lymphomaGraph, verbose = FALSE, presolve = FALSE, solverChoice = "CPLEXAPI")
    cplexapiSolution = steinObject$findSingleSteinerSolution()

    expect_equal(vcount(lpSolveSolution), vcount(cplexapiSolution))
    expect_equal(sum(V(lpSolveSolution)$nodeScore), sum(V(cplexapiSolution)$nodeScore))
    expect_true(is.connected(cplexapiSolution))

})

test_that("Compare Rcbc result", {
  
  
  if(!"rcbc" %in% .packages(all.available = TRUE)) skip("rcbc not installed")
  
  steinObject = nodeCentricSteinerTreeProblem$new(lymphomaGraph, presolve = FALSE, verbose = FALSE, solverChoice = "RCBC")
  cbcSolution = steinObject$findSingleSteinerSolution()
  
  expect_equal(vcount(lpSolveSolution), vcount(cbcSolution))
  expect_equal(sum(V(lpSolveSolution)$nodeScore), sum(V(cbcSolution)$nodeScore))
  expect_true(is.connected(cbcSolution))
  
})



test_that("Compare Rglpk result", {


  if(!"Rglpk" %in% .packages(all.available = TRUE)) skip("Rglpk not installed")

    steinObject = nodeCentricSteinerTreeProblem$new(lymphomaGraph, presolve = FALSE, verbose = FALSE, solverChoice = "RGLPK")
    glpkSolution = steinObject$findSingleSteinerSolution()

    expect_equal(vcount(lpSolveSolution), vcount(glpkSolution))
    expect_equal(sum(V(lpSolveSolution)$nodeScore), sum(V(glpkSolution)$nodeScore))
    expect_true(is.connected(glpkSolution))

})

test_that("Compare lpsymphony result", {

  if(!"lpsymphony" %in% .packages(all.available = TRUE)) skip("lpsymphony not installed")

    steinObject = nodeCentricSteinerTreeProblem$new(lymphomaGraph, verbose = FALSE, presolve = FALSE, solverChoice = "LPSYMPHONY")
    symphonySolution = steinObject$findSingleSteinerSolution()

    expect_equal(vcount(lpSolveSolution), vcount(symphonySolution))
    expect_equal(sum(V(lpSolveSolution)$nodeScore), sum(V(symphonySolution)$nodeScore))
    expect_true(is.connected(symphonySolution))

})
