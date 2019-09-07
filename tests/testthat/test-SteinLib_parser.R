context('Testing the SteinLib parser')

test_that("Checking the MStTP STP  file format reader on SteinLib for GENE problem gene 42",{

  gene42_igraph <- readMStTPgraph('./testData/GENE/gene42.stp')
  expect_equal(vcount(gene42_igraph),335)
  expect_equal(ecount( as.undirected(gene42_igraph)),456)
  expect_equal( length(which( V(gene42_igraph)$isTerminal)),43)
  expect_equal( length( unique(E(gene42_igraph)$edgeWeight )),1)
})

test_that("Checking the MStTP STP  file reader on SteinLib for GENE problem gene gene61a",{

  gene61a_igraph <- readMStTPgraph('./testData/GENE/gene61a.stp')
  expect_equal(vcount(gene61a_igraph),395)
  expect_equal(ecount( as.undirected(gene61a_igraph)),512)
  expect_equal( length(which( V(gene61a_igraph)$isTerminal)),82)
  expect_equal( length( unique(E(gene61a_igraph)$edgeWeight )),1)
})

test_that("Checking the MStTP STP file format reader on SteinLib for GENE problem gene gene61b",{

  gene61b_igraph <- readMStTPgraph('./testData/GENE/gene61b.stp')
  expect_equal(vcount(gene61b_igraph),570)
  expect_equal(ecount( as.undirected(gene61b_igraph)),808)
  expect_equal( length(which( V(gene61b_igraph)$isTerminal)),82)
  expect_equal( length( unique(E(gene61b_igraph)$edgeWeight )),1)
})

test_that("Checking the MWCS STP file format reader on SteinLib for ACTMOD problem gene lymphoma.stp",{

  lymphoma_igraph <- readMWCSgraph('./testData/ACTMOD/lymphoma.stp')
  expect_equal(vcount(lymphoma_igraph),2034)
  expect_equal(ecount( lymphoma_igraph),7756)
  expect_false(is.directed(lymphoma_igraph))

  #Check that nodes have be attributed the correct scores
  expect_equal(V(lymphoma_igraph)[name == '575']$nodeScore,-6.14,tolerance = 0.01)
  expect_equal(V(lymphoma_igraph)[name == '1319']$nodeScore,-7.23,tolerance = 0.01)
  expect_equal(V(lymphoma_igraph)[name == '554']$nodeScore,-6.94,tolerance = 0.01)
})

test_that("Checking the MWCS STP file format reader on SteinLib for ACTMOD problem gene HCMV.stp",{

  HCMV_igraph <- readMWCSgraph('./testData/ACTMOD/HCMV.stp')
  expect_equal(vcount(HCMV_igraph),3863)
  expect_equal(ecount( HCMV_igraph),29293)
  expect_false(is.directed(HCMV_igraph))

  #Check that nodes have be attributed the correct scores
  expect_equal(V(HCMV_igraph)[name == '1359']$nodeScore,-3.462,tolerance = 0.01)
  expect_equal(V(HCMV_igraph)[name == '3107']$nodeScore,-1.296,tolerance = 0.01)
  expect_equal(V(HCMV_igraph)[name == '551']$nodeScore,-1.648,tolerance = 0.01)
})

test_that("Checking the MWCS STP file format from SteinLib for GENE problem gene gene42.stp",{

  gene42_in <- readMStTPgraph('./testData/GENE/gene42.stp')
  
  temp_file_gene42 <- tempfile(pattern = "gene42")
  writeMStTPfile_heinzFormat(gene42_in, temp_file_gene42)
  
  gene42_igraph_alt <- readMWCSgraph(temp_file_gene42)

  expect_identical(vcount(gene42_in),vcount(gene42_igraph_alt))
  expect_identical(ecount(gene42_in),ecount(gene42_igraph_alt))
  
  expect_identical(length( V(gene42_in)[isTerminal == TRUE]$name ), 
                   length(V(gene42_igraph_alt)[nodeScore > 0]$name))
})
