context('Testing graph condensation procedures')

test_that("Check presolve on lymphoma",{

  lymphoma_MWCS = nodeCentricSteinerTreeProblem$new( lymphomaGraph, verbose = FALSE)$findSingleSteinerSolution()

  condensedGraph = condenseSearchGraph(lymphomaGraph)
  condensed_MWCS = nodeCentricSteinerTreeProblem$new( condensedGraph, verbose = FALSE, presolveGraph = FALSE)$findSingleSteinerSolution()

  expect_equal(sum(V(condensed_MWCS)$nodeScore), sum(V(lymphoma_MWCS)$nodeScore))

  expect_equal(sort(V(uncondenseGraph(condensed_MWCS))$name), sort(V(lymphoma_MWCS)$name))
})

test_that("Check presolve on karate graph",{

  karateGraph_MStP = nodeCentricSteinerTreeProblem$new( karateGraph, verbose = FALSE)$findSingleSteinerSolution()

  condensedGraph = condenseSearchGraph(karateGraph)
  condensed_MStTP = nodeCentricSteinerTreeProblem$new( condensedGraph, verbose = FALSE, presolveGraph = FALSE)$findSingleSteinerSolution()

  expect_equal(sort(V(condensed_MStTP)[isTerminal]$name),sort(V(karateGraph_MStP)[isTerminal]$name))

  expect_equal(sort(V(uncondenseGraph(condensed_MStTP))$name),sort(V(karateGraph_MStP)$name))
})


test_that("Check that a fully connected input graph is dealt with properly",{

  connectedSubNetNodes = c('q','g','f')

  V(karateGraph)$isTerminal = FALSE
  V(karateGraph)[name %in% connectedSubNetNodes]$isTerminal = TRUE

  noPresolveRes = nodeCentricSteinerTreeProblem$new( karateGraph, verbose = FALSE, presolveGraph = FALSE)$findSingleSteinerSolution()
  presolvedRes = nodeCentricSteinerTreeProblem$new( karateGraph, verbose = FALSE, presolveGraph = TRUE)$findSingleSteinerSolution()

  presolveStein = nodeCentricSteinerTreeProblem$new( karateGraph, verbose = FALSE, presolveGraph = TRUE)

  expect_equal(vcount(noPresolveRes),vcount(presolvedRes))
  expect_equal(ecount(noPresolveRes),ecount(presolvedRes))
  expect_equal(sort(V(noPresolveRes)$name),sort(V(presolvedRes)$name))
})
