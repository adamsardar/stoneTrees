library(einSteiner)
context("Checking that the APSP function works as expected")

test_that("Look at APSP for gene42",{ 
  
  seedNodes <- V(gene42_igraph)[isTerminal]$name
  
  gene42APSPmod <- calculateAPSPmodule(seedNodes,gene42_igraph)
  
  expect_true(is.igraph(gene42APSPmod))
  expect_true(is.connected(gene42APSPmod))
  expect_error(calculateAPSPmodule(c("not a","set","of","seeds"),gene42_igraph))
  expect_warning(calculateAPSPmodule(c(seedNodes,"not a","set","of","seeds"),gene42_igraph))
  
  gene42_igraph %<>% delete_vertex_attr("name")
  expect_error(calculateAPSPmodule(seedNodes,gene42_igraph))
})

NAnamefGraph <- readRDS("./testData/NAnamedNetword.RDS")

test_that("Check Behaviour On A Graph With NA Values - with a seed set containing an NA node",{
  
  seedNodes <- c("Y","J","C",NA)
  
  apspMod <- calculateAPSPmodule(seedNodes,NAnamefGraph)

  expect_warning(calculateAPSPmodule(seedNodes,NAnamefGraph,omitNA = FALSE))
  
  expect_true(is.igraph(apspMod))
  expect_true(is.connected(apspMod))  

  expect_false( any( is.na(V(apspMod)$name)),info = "There shouldn't be any NA nodes in the answer here (looking at the input graph)")
})

test_that("Check Behaviour On A Graph With NA Values - with a path known to include an NA-named node",{
  
  seedNodes <- c("B","H","Z","J")
  
  apspMod <- calculateAPSPmodule(seedNodes,NAnamefGraph)
  
  expect_true(is.igraph(apspMod))
  expect_true(is.connected(apspMod))  
  expect_equal(vcount(apspMod),6, info = "Correct module size is 6, but if too many NA's are included in the subgraph induction step we shall have more nodes that we want")
  expect_equal(vcount( induced_subgraph(apspMod,V(apspMod)[is.na(name)]) ),1, info="There should only be one NA named node in this graph") 
})
