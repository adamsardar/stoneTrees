context('Testing the Heinz solver against some SteinLib instances')

library(data.table)

test_that("Trying Heinz MStTP against the gene42.stp instance in the GENE SteinLib set",{

  if(!Sys.which("heinz") %like% "heinz"){skip("Can't run this test without heinz installed")}

  gene42_igraph <- readMStTPgraph('./testData/GENE/gene42.stp')

  suppressWarnings(gene42_MSTP <- findHeinzMStTP(gene42_igraph))
  expect_true(is.connected(gene42_MSTP))
  expect_true(vcount(gene42_MSTP) <= 1.05*126) # Within 5% of the known optimum

  expect_false(is.null(graph_attr(gene42_MSTP)$SearchNetwork))
  expect_equal(graph_attr(gene42_MSTP)$SearchNetwork,"gene42_igraph")

  expect_false(is.null(graph_attr(gene42_MSTP)$SamplingMethod))
  expect_equal(graph_attr(gene42_MSTP)$SamplingMethod,"Minimum Spanning Tree")
})

test_that("Trying Heinz MStTP against the gene61a.stp instance in the GENE SteinLib set",{

  if(!Sys.which("heinz") %like% "heinz"){skip("Can't run this test without heinz installed")}

  gene61a_igraph <- readMStTPgraph('./testData/GENE/gene61a.stp')

  suppressWarnings(gene61a_MSTP <- findHeinzMStTP(gene61a_igraph))
  expect_true(is.connected(gene61a_MSTP))
  expect_true(vcount(gene61a_MSTP) <= 1.05*205) # Within 5% of the known optimum

  expect_false(is.null(graph_attr(gene61a_MSTP)$SearchNetwork))
  expect_equal(graph_attr(gene61a_MSTP)$SearchNetwork,"gene61a_igraph")

  expect_false(is.null(graph_attr(gene61a_MSTP)$SamplingMethod))
  expect_equal(graph_attr(gene61a_MSTP)$SamplingMethod,"Minimum Spanning Tree")
})

test_that("Trying Heinz MStTP against the gene61f.stp instance in the GENE SteinLib set",{

  if(!Sys.which("heinz") %like% "heinz"){skip("Can't run this test without heinz installed")}

  gene61f_igraph <- readMStTPgraph('./testData/GENE/gene61f.stp')

  suppressWarnings(gene61f_MSTP <- findHeinzMStTP(gene61f_igraph))
  expect_true(is.connected(gene61f_MSTP))
  expect_true(vcount(gene61f_MSTP) <= 1.05*198) # Within 5% of the known optimum
})

test_that("Trying Heinz MWCS against the HCMV.stp instance in the ACTMOD SteinLib set",{

  if(! Sys.which("heinz") %like% "heinz"){skip("Can't run this test without heinz installed")}

  HCMV_igraph <- readMWCSgraph('./testData/ACTMOD/HCMV.stp')

  suppressWarnings(HCMV_MWCS <- findHeinzMWCS(HCMV_igraph))

  expect_equal(vcount(HCMV_MWCS),17)
  expect_true(is.connected(HCMV_MWCS))

  expect_false(is.null(graph_attr(HCMV_MWCS)$SearchNetwork))
  expect_equal(graph_attr(HCMV_MWCS)$SearchNetwork,"HCMV_igraph")

  expect_false(is.null(graph_attr(HCMV_MWCS)$SamplingMethod))
  expect_equal(graph_attr(HCMV_MWCS)$SamplingMethod,"Minimum Spanning Tree Of Maximum-Weight Connected Subgraph")

  suppressWarnings(HCMV_MWCS <- findHeinzMWCS(HCMV_igraph,induceSubgraph = TRUE))

  expect_false(is.null(graph_attr(HCMV_MWCS)$SamplingMethod))
  expect_equal(graph_attr(HCMV_MWCS)$SamplingMethod,"Induced Subgraph Of Maximum-Weight Connected Subgraph")
})

test_that("Trying Heinz MWCS against the lymphoma.stp instance in the ACTMOD SteinLib set",{

  if(!Sys.which("heinz") %like% "heinz"){skip("Can't run this test without heinz installed")}

  lymphoma_igraph <- readMWCSgraph('./testData/ACTMOD/lymphoma.stp')

  suppressWarnings(lymphoma_MWCS <- findHeinzMWCS(lymphoma_igraph))

  expect_equal(vcount(lymphoma_MWCS),46)
  expect_true(is.connected(lymphoma_MWCS))
})
