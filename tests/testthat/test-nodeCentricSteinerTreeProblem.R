context("Inspecting node-centric steiner tree object creation")

library(igraph)
library(data.table)

#Define a utility function - the allows us to perform a load of checks across networks of different kinds
inspectNodeCentricSteinerTreeObjectCreation <- function(testGraph){

  test_that("Inspect object creation", {

    if(is.directed(testGraph)){

      expect_warning(testSteinObject <- nodeCentricSteinerTreeProblem$new(testGraph, verbose = FALSE),
                     regexp = "directed",
                     info = "Directed graphs are cast to undirected with a warning")
    }

    testGraph <- as.undirected(testGraph)

    expect_message(testSteinObject <- nodeCentricSteinerTreeProblem$new(testGraph, verbose = TRUE),
                   regexp = "constraints",
                   info = "verbose flag should control diagnostic printing")

    expect_silent(nodeCentricSteinerTreeProblem$new(testGraph, verbose = FALSE))

    # Inspect fixed terminal constraints
    expect_equal( dim(testSteinObject$getFixedTerminalConstraints()$variables),
                  c( nrow(testSteinObject$getNodeDT()[isTerminal == TRUE]), vcount(testGraph)),
                  info = "There should be an equal number of fixed termianls and respective constraints (with nVertices varaibles)")

    expect_equal( nrow(testSteinObject$getFixedTerminalConstraints()$variables),
                  length(testSteinObject$getFixedTerminalConstraints()$directions),
                  info = "The number of variables and inequality directions must match")

    expect_equal(testSteinObject$getFixedTerminalConstraints()$directions,
                 rep("==", nrow(testSteinObject$getNodeDT()[isTerminal == TRUE])),
                 info = "Fixed terminals are equality constraints")

    expect_equal( nrow(testSteinObject$getFixedTerminalConstraints()$variables),
                  length(testSteinObject$getFixedTerminalConstraints()$rhs),
                  info = "The number of variables and inequality values must match")

    expect_equal(testSteinObject$getFixedTerminalConstraints()$rhs,
                 rep(1, nrow(testSteinObject$getNodeDT()[isTerminal == TRUE])),
                 info = "Fixed terminals are simple 1=1 equalities")

    expect_equivalent( rowSums(testSteinObject$getFixedTerminalConstraints()$variables),
                       rep(1, nrow(testSteinObject$getNodeDT()[isTerminal == TRUE]) ),
                       info = "Fixed terminal constraints are simple - a single 1 in the terminal variables")

    # Inspect node degree inequalities
    expect_equal( dim(testSteinObject$getNodeDegreeConstraints()$variables),
                  c(vcount(testGraph), vcount(testGraph)),
                  info = "The should be a constraint for every node and a varabile for every node, hence |V|x|V|")

    expect_equal( nrow(testSteinObject$getNodeDegreeConstraints()$variables),
                  length(testSteinObject$getNodeDegreeConstraints()$directions),
                  info = "The number of variables and inequality directions must match")

    expect_equal( nrow(testSteinObject$getNodeDegreeConstraints()$variables),
                  length(testSteinObject$getNodeDegreeConstraints()$rhs),
                  info = "The number of variables and inequality values must match")

    expect_true( all(rowSums(testSteinObject$getNodeDegreeConstraints()$variables) < degree(testGraph)) )

    expect_true( all(diag(testSteinObject$getNodeDegreeConstraints()$variables) %in% c(-1,-2)),
                 info = "All nodes must be -1 or -2 on diagonal")

    expect_true( all(diag(testSteinObject$getNodeDegreeConstraints()$variables)[testSteinObject$getTerminals()$terminals] == -1),
                 info = "All terminal nodes must be -1 on the diagonal (tree can stop there)")

    # Inspect two-cycle inequalities
    expect_equal( dim(testSteinObject$getTwoCycleConstraints()$variables),
                  c( nrow(merge(testSteinObject$getNodeDT()[nodeScore > 0], testSteinObject$getEdgeDT(), by.x = ".nodeID", by.y = "fromNodeID")),
                     vcount(testGraph) ),
                  info = "Two cycle inequalities refer to potential terminals, hence there should be as many constraints as there are positive node scores")

    expect_equal( nrow(testSteinObject$getTwoCycleConstraints()$variables),
                  length(testSteinObject$getTwoCycleConstraints()$directions),
                  info = "The number of variables and inequality directions must match")

    expect_equal( nrow(testSteinObject$getTwoCycleConstraints()$variables),
                  length(testSteinObject$getTwoCycleConstraints()$rhs),
                  info = "The number of variables and inequality values must match")

    expect_equivalent(rowSums(testSteinObject$getTwoCycleConstraints()$variables),
                      rep(0, nrow(merge(testSteinObject$getNodeDT()[nodeScore > 0], testSteinObject$getEdgeDT(), by.x = ".nodeID", by.y = "fromNodeID"))),
                      info = "Two-cycle inqualities are a seies of 1,-1 pairs: one for each edge out from a potential terminal")

    expect_equivalent(testSteinObject$getTwoCycleConstraints()$rhs,
                      rep(0, nrow(merge(testSteinObject$getNodeDT()[nodeScore > 0], testSteinObject$getEdgeDT(), by.x = ".nodeID", by.y = "fromNodeID"))),
                      info = "Two-cycle inqualities are a seies of 1,-1 pairs: one for each edge out from a potential terminal")

    expect_equivalent(testSteinObject$getTwoCycleConstraints()$directions,
                      rep(">=", nrow(merge(testSteinObject$getNodeDT()[nodeScore > 0], testSteinObject$getEdgeDT(), by.x = ".nodeID", by.y = "fromNodeID"))),
                      info = "Two-cycle inqualities (>=) are a seies of 1,-1 pairs: one for each edge out from a potential terminal")

    # Check that connectivity constraints are initialised correctly
    expect_true(is.numeric(testSteinObject$getConnectivityConstraints()$rhs))
    expect_true(length(testSteinObject$getConnectivityConstraints()$rhs) == 0)

    expect_true(is.character(testSteinObject$getConnectivityConstraints()$directions))
    expect_true(length(testSteinObject$getConnectivityConstraints()$directions) == 0)

    expect_true( all(dim(testSteinObject$getConnectivityConstraints()$variables) == c(0,vcount(testGraph))) )

    expect_equal( nrow(testSteinObject$getConnectivityConstraints()$variables),
                  length(testSteinObject$getConnectivityConstraints()$directions),
                  info = "Connectivity constraints are added lazily - there should be none at the start!")

    expect_equal( nrow(testSteinObject$getConnectivityConstraints()$variables),
                  length(testSteinObject$getConnectivityConstraints()$rhs),
                  info = "Connectivity constraints are added lazily - there should be none at the start!")

  })

  return(invisible())
}

test_that("Ensure that incorrect constructor inputs are not tolerated", {

  expect_error(nodeCentricSteinerTreeProblem$new(letters),
               regexp = "igraph",
               info = "Class should fail when a non-igraph object is given")

  expect_error(nodeCentricSteinerTreeProblem$new(lymphomaGraph, solverChoice = "noSuchSolver"),
               regexp = "solver",
               info = "Class should fail when a nonsense solver choice is given")

  expect_error(nodeCentricSteinerTreeProblem$new(lymphomaGraph, verbose = "notAboolean"),
               regexp = "verbose",
               info = "Class should fail when a nonsense verbosity choice is given")
})

# Seed based - small
inspectNodeCentricSteinerTreeObjectCreation(karateGraph)

# nodeScore based - medium sized
inspectNodeCentricSteinerTreeObjectCreation(lymphomaGraph)

# seed based - medium sized
inspectNodeCentricSteinerTreeObjectCreation(gene42_igraph)

# TODO - include more Stein Lib examples
