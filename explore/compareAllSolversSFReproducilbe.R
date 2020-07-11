devtools::load_all()
library(igraph)
library(purrr)

library(data.table)
library(ggplot2)




# reformat node score for MStT (not MWCS) and add Terminal attributes
fixedTerminalLymphomaGraph <- lymphomaGraph
V(fixedTerminalLymphomaGraph)$isTerminal <- FALSE
V(fixedTerminalLymphomaGraph)[nodeScore > 0]$isTerminal <- TRUE
V(fixedTerminalLymphomaGraph)$nodeScore <- -1

#check that all solvers are working
##nodeCentricSteinerTreeProblem
#allSolversSolution <- map(stoneTrees_solvers, ~{ nodeCentricSteinerTreeProblem$new(fixedTerminalLymphomaGraph,
#                                                                                   solverChoice = .x)$findSingleSteinerSolution() } )
#allSolversSolution

##nodeCentricSteinerForestProblem
#allSolversSolution <- map(stoneTrees_solvers, ~{ nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph,
#                                                                                     solverChoice = .x)$sampleMultipleBootstrapSteinerSolutions(nBootstraps = 1)  } )
#map(allSolversSolution, ~{.x$getBootstrapSolutionPoolGraphs()})

set.seed(2345)

#set the number of trials
numberOfTrial <- 30
bootstrapIteraction <- 100
nsuboptimalSol <- 5

Solver <- "rcbc"

steinForestSuboptimalBenchDT <- map(1:numberOfTrial, function(i){


  print(i)

  # Note, The two process are not solving the exact same problem. Stochastic algorithms are a pain

time <- system.time( SteinForSolverX <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph,
                                                                               solverChoice = Solver,
                                                                               verbose = TRUE,
                                                                               solverTrace = 0)$sampleMultipleBootstrapSteinerSolutions(nBootstraps = bootstrapIteraction,
                                                                                                                                        maxItr = nsuboptimalSol) )

    SolutionPool <- SteinForSolverX$getBootstrapSolutionPoolGraphs(collapseSols = FALSE)


    save(SolutionPool, file = paste0("./results/Reproducible_steinForestSuboptimalSolver_",Solver,"_Trial", i,".RData"))



    Sizes <- map_int(SolutionPool, vcount)

  return(data.table(trial = paste(i),
                   solver = paste(Solver),
                    time = time["elapsed"],
                    vcount = vcount(SteinForSolverX$getBootstrapSolutionPoolGraphs()),
                    modulesVcount = list(Sizes[-length(Sizes)]),
                    Niteration = length(Sizes[-length(Sizes)])))



}) %>% rbindlist()

#test1 <- rbindlist(steinForestSuboptimalBenchDT)


fwrite(steinForestSuboptimalBenchDT, paste0("./results/Reproducible_steinForestSuboptimalBench_",Solver, ".tsv"))
