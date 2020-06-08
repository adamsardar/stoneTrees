library(stoneTrees)
library(igraph)
library(microbenchmark)

fixedTerminalLymphomaGraph <- lymphomaGraph
V(fixedTerminalLymphomaGraph)$isTerminal <- FALSE
V(fixedTerminalLymphomaGraph)[nodeScore > 0]$isTerminal <- TRUE
V(fixedTerminalLymphomaGraph)$nodeScore <- -1

# I notice that the Steiner forest is actually much harder than MWCS to solve - probably the symmetry of having no node-weights

testSteinFor <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rcbc")

set.seed(2345)
solDetails <- list()

for(i in 1:100){
  
  print(i)

  #Ensure that our solvers are both solving the same type of problem
  # We can't just run two bootstrap Stiener forest routines 
  testSteinFor$.__enclos_env__$private$resampleFixedTerminals()

  testSteinFor$.__enclos_env__$private$solver <- "RGLPK"
  rglpkTime <- system.time(RGLPKsol <- testSteinFor$findSingleSteinerSolution())
  
  testSteinFor$.__enclos_env__$private$solver <- "RCBC"
  rcpcTime <- system.time(RCPCsol <- testSteinFor$findSingleSteinerSolution()) 

  solDetails %<>%
    c(list(data.table(RCPCtime = rcpcTime["elapsed"], 
               RCPCvcount = vcount(RCPCsol),
               RGLPKtime = rglpkTime["elapsed"],
               RGLPKvcount = vcount(RGLPKsol),
               nFixedTerminals = length(V(RCPCsol)[isTerminal]))))
}


rbindlis

system.time( nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rcbc")$findSingleSteinerSolution() )
system.time( nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rglpk")$findSingleSteinerSolution() )

system.time()
system.time( )

bench <- microbenchmark(
  {set.seed(1234); nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rcbc")$sampleMultipleBootstrapSteinerSolutions()$getBootstrapSolutionPoolGraphs()},
  {set.seed(1234);  nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rglpk")$sampleMultipleBootstrapSteinerSolutions()$getBootstrapSolutionPoolGraphs()}
)



#Fixed seeds? Averaged over set