context('Testing graph condensation procedures')

test_that("Check presolve on lymphoma",{

  lymphoma_igraph <- readMWCSgraph('./testData/ACTMOD/lymphoma.stp')
  lymphoma_MWCS <- solveNodeCentricSteinerTreeProblem(lymphoma_igraph,verbosity=0)

  condensedGraph <- condenseSearchGraph(lymphoma_igraph)
  condensed_MWCS <- solveNodeCentricSteinerTreeProblem(condensedGraph,verbosity=0,presolveGraph=FALSE)

  expect_equal(sum(V(condensed_MWCS)$nodeScore),sum(V(lymphoma_MWCS)$nodeScore))
  expect_equal(sort(V(condensed_MWCS)$name),sort(V(lymphoma_MWCS)$name))
})

test_that("Check presolve on karate graph",{

  data("karateGraph")

  karateGraph_MStP <- solveNodeCentricSteinerTreeProblem(karateGraph,verbosity=0)

  condensedGraph <- condenseSearchGraph(karateGraph)
  condensed_MStTP <- solveNodeCentricSteinerTreeProblem(condensedGraph,verbosity=0,presolveGraph = FALSE)

  expect_equal(sort(V(condensed_MStTP)[isTerminal]$name),sort(V(karateGraph_MStP)[isTerminal]$name))
  expect_equal(sort(V(condensed_MStTP)$name),sort(V(karateGraph_MStP)$name))
})


test_that("Check that a fully connected input graph is dealt with properly",{

  data("karateGraph")

  connectedSubNetNodes <- c('q','g','f')

  V(karateGraph)$isTerminal <- FALSE
  V(karateGraph)[name %in% connectedSubNetNodes]$isTerminal <- TRUE

  noPresolveRes <- solveNodeCentricSteinerTreeProblem(karateGraph, presolveGraph = FALSE)
  presolvedRes <- solveNodeCentricSteinerTreeProblem(karateGraph, presolveGraph = TRUE)

  expect_equal(vcount(noPresolveRes),vcount(presolvedRes))
  expect_equal(ecount(noPresolveRes),ecount(presolvedRes))
  expect_equal(sort(V(noPresolveRes)$name),sort(V(presolvedRes)$name))
})
