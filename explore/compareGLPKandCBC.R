devtools::load_all()
library(igraph)
#library(microbenchmark)


fixedTerminalLymphomaGraph <- lymphomaGraph
V(fixedTerminalLymphomaGraph)$isTerminal <- FALSE
V(fixedTerminalLymphomaGraph)[nodeScore > 0]$isTerminal <- TRUE
V(fixedTerminalLymphomaGraph)$nodeScore <- -1

# I notice that the Steiner forest is actually much harder than MWCS to solve - probably the symmetry of having no node-weights

# 10 bootstraps, 5x suboptimal
glpktim <- system.time( SteinForGLPK <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, solverChoice = "RGLPK")$sampleMultipleBootstrapSteinerSolutions(nBootstraps = 100, maxItr = 5) )
## Takes around a minute using RGLPK as the solver
cbctime <- system.time( SteinForCBC <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, solverChoice = "RCBC")$sampleMultipleBootstrapSteinerSolutions(nBootstraps = 100, maxItr = 5) )



testSteinFor <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = TRUE, solverChoice = "rcbc")

set.seed(2345)

fwrite(data.table(RCBCtime = character(), 
                  RCBCvcount = integer(),
                  RGLPKtime = character(),
                  RGLPKvcount = integer(),
                  nFixedTerminals = integer()),
  "/tmp/solDetails2.tsv", sep = "\t")

for(i in 1:100){

  print(i)

  #Ensure that our solvers are both solving the same type of problem
  # We can't just run two bootstrap Steiner forest routines 
  testSteinFor$.__enclos_env__$private$resampleFixedTerminals()

  testSteinFor$.__enclos_env__$private$flushConnectivityConstraints()
  testSteinFor$.__enclos_env__$private$solver <- "RGLPK"
  rglpkTime <- system.time(RGLPKsol <- testSteinFor$findSingleSteinerSolution())

  testSteinFor$.__enclos_env__$private$flushConnectivityConstraints()  
  testSteinFor$.__enclos_env__$private$solver <- "RCBC"
  rcbcTime <- system.time(RCBCsol <- testSteinFor$findSingleSteinerSolution()) 
  
  fwrite(data.table(RCBCtime = rcbcTime["elapsed"], 
               RCBCvcount = vcount(RCBCsol),
               RGLPKtime = rglpkTime["elapsed"],
               RGLPKvcount = vcount(RGLPKsol),
               nFixedTerminals = length(V(RCBCsol)[isTerminal])),
         file = "/tmp/solDetails2.tsv", sep = "\t", quote = FALSE, append = TRUE)
}


solDetailsDT <- fread("/tmp/solDetails.tsv")

solDetailsDT  %>%
  .[,trial := .I] %>%
  melt(id.vars = "trial") %>%
  .[variable %like% "time"] %>%
  ggplot(aes(x = trial, y = value)) +
  geom_point(aes(colour = variable)) +
  theme_bw() +
  scale_colour_manual(values = c("dodgerblue","black"))


solDetailsDT %>%
  .[,trial := .I] %>%
  melt(id.vars = "trial") %>%
  .[variable %like% "time"] %>%
  ggplot(aes(x = variable, y = value/60)) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(aes(colour = variable)) +
  scale_y_log10() +
  theme_bw() + theme(legend.position = "none") +
  labs(x = "Solver", y = "Time (minues) [log scale]", title = "Comparison of CBC and GLPK solving single iterations ofthe Steiner Forest problem")

solDetailsDT[,.N,by = .(RCPCvcount <= RGLPKvcount)]
solDetailsDT[,.N,by = .(RCPCvcount >= RGLPKvcount)]

system.time( nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rcbc")$findSingleSteinerSolution() )
system.time( nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rglpk")$findSingleSteinerSolution() )


bench <- microbenchmark(
  {set.seed(1234); nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rcbc")$sampleMultipleBootstrapSteinerSolutions()$getBootstrapSolutionPoolGraphs()},
  {set.seed(1234);  nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rglpk")$sampleMultipleBootstrapSteinerSolutions()$getBootstrapSolutionPoolGraphs()}
)


## Sample from bootstrap

set.seed(2345)

fwrite(data.table(RCBCtime = character(), 
                  RCBCvcount = integer(),
                  RGLPKtime = character(),
                  RGLPKvcount = integer(),
                  nFixedTerminals = integer()),
       "~/steinForestBench.tsv", sep = "\t")

for(i in 1:100){
  
  print(i)
  
  #Ensure that our solvers are both solving the same type of problem
  # We can't just run two bootstrap Steiner forest routines 
  rglpkTime <- system.time(RGLPKsol <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rglpk")$sampleMultipleBootstrapSteinerSolutions())
  
  rcbcTime <- system.time(RCBCsol <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, verbose = FALSE, solverChoice = "rcbc")$sampleMultipleBootstrapSteinerSolutions())
            
  
                
  fwrite(data.table(RCBCtime = rcbcTime["elapsed"], 
                    RCBCvcount = vcount(RCBCsol$getBootstrapSolutionPoolGraphs()),
                    RGLPKtime = rglpkTime["elapsed"],
                    RGLPKvcount = vcount(RGLPKsol),
                    
                    
         file = "~/steinForestBench.tsv", sep = "\t", quote = FALSE, append = TRUE)
}


steinForBench <- fread("~/steinForestBench.tsv")

steinForBench[, trial := .I]

pTime <- steinForBench[,.(RCBCtime, RGLPKtime, trial)] %>% 
  melt(id.vars = "trial") %>%
  ggplot(aes(x = variable, y = as.numeric(value)/60)) +
  geom_boxplot(aes(colour = variable)) +
  geom_jitter(aes(colour = variable)) +
  theme_bw() +
  labs(x = "Solver", y = "Time (min)")

pTime + scale_y_log10()


steinForBench[,.(RCBCvcount, RGLPKvcount, trial)] %>%
  melt(id.vars = "trial") %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot(aes(colour = variable)) +
  geom_jitter(aes(colour = variable))

# Networks are bigger. Is this because RCBC is creating worse solutions?


set.seed(2345)

fwrite(data.table(RCBCtime = character(), 
                  RCBCvcount = integer(),
                  RGLPKtime = character(),
                  RGLPKvcount = integer(),
                  nFixedTerminals = integer()),
       "~/steinForestSuboptimalBench.tsv", sep = "\t")


for(i in 1:15){
  
  print(i)

  # 10 bootstraps, 5x suboptimal
  # Note, The two process are not solving the exact same problem. Stochastic algorithms are a pain
  
message("Starting GLPK")
  glpktim <- system.time( SteinForGLPK <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, solverChoice = "RGLPK", verbose = FALSE)$sampleMultipleBootstrapSteinerSolutions(nBootstraps = 100, maxItr = 5) )
message("Completed GLPK")

message("Starting CBC")
  cbctime <- system.time( SteinForCBC <- nodeCentricSteinerForestProblem$new(fixedTerminalLymphomaGraph, solverChoice = "RCBC", verbose = FALSE)$sampleMultipleBootstrapSteinerSolutions(nBootstraps = 100, maxItr = 5) )
message("Completed CBC")  

  cbcSizes <- sapply(SteinForCBC$getBootstrapSolutionPoolGraphs(collapseSols = FALSE), vcount)
  glpkSizes <- sapply(SteinForGLPK$getBootstrapSolutionPoolGraphs(collapseSols = FALSE), vcount)

  # The final item in the set is the standard Steiner Tree solution, omit
  glpkSizes[-length(glpkSizes)] %>% length  
  cbcSizes[-length(cbcSizes)] %>% length
  
  
  fwrite(data.table(CBCtime = cbctime["elapsed"], 
                    CBCvcount = vcount(SteinForCBC$getBootstrapSolutionPoolGraphs()),
                    CBCmodulesVcount = list( glpkSizes[-length(glpkSizes)]),
                    GLPKtime = glpktim["elapsed"],
                    GLPKvcount = vcount(SteinForGLPK$getBootstrapSolutionPoolGraphs()),
                    GLPKmodulesVcount = list( cbcSizes[-length(cbcSizes)])
                    ),
                    
                    
                    file = "~/steinForestSuboptimalBench.tsv", sep = "\t", quote = FALSE, append = TRUE)
}


steinForestSuboptimalBenchDT <- fread("~/steinForestSuboptimalBench.tsv")
steinForestSuboptimalBenchDT[, trial := .I]

# Much faster
steinForestSuboptimalBenchDT[,.(CBCtime, GLPKtime, trial)] %>%
  melt(id.vars = "trial") %>%
  ggplot(aes(x = variable, y = value/60, colour = variable)) +
  geom_boxplot() +
  geom_jitter() +
  scale_y_continuous( breaks = seq(0,100,5)) +
  theme_bw(base_size = 20) +
  labs(x = "solver", y = "Time (min)",
       title = "Time taken to solve Steiner forest problem on lymphoma graph",
       subtitle = "100 bootstraps, collecting up to 5 degenerate solutions on each run.")

# But much larger
steinForestSuboptimalBenchDT[,.(CBCvcount, GLPKvcount, trial)]  %>%
  melt(id.vars = "trial") %>%
  ggplot(aes(x = variable, y = value, colour = variable)) +
  geom_jitter() +
  geom_boxplot() + 
  theme_bw(base_size = 20) +
  labs(x = "solver", y = "Aggregated module size",
       title = "Size of Steiner forest solution (i.e. all of the subproblem graphs combined)",
       subtitle = "CBC solutions are much larger than GLPK")



GLPKmodules <- steinForestSuboptimalBenchDT[,.(solver = "GLPK", size = unlist(strsplit(GLPKmodulesVcount, "\\|"))), by = trial]
CBCmodules <- steinForestSuboptimalBenchDT[,.(solver = "CBC", size = unlist(strsplit(CBCmodulesVcount, "\\|"))), by = trial]


solutionModules <- rbind(CBCmodules, GLPKmodules)
solutionModules[,size := as.integer(size)]

solutionModules %>%
  ggplot(aes(x = solver, y = size, colour = solver)) +
  geom_jitter() +
  geom_boxplot() + 
  theme_bw(base_size = 20) +
  labs(x = "solver", y = "Module size (nodes) of a single Steiner forest sub-problem",
       title = "Individual tree solutions to resampled terminals",
       subtitle = "Collected across 15 runs of 100x bootstraps\nCBC solutions are, on average, the same size as GLPK solutions.")


# But there are more CBC solutions collected

solutionModules[,.N,by = .(trial, solver)] %>%
  ggplot(aes(x = solver, y = N, colour = solver)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 20) +
  labs(title = "Size of solution pool (i.e. the number of distinct solutions to Steiner forest subproblems) collected",
       subtitle = "Far more solutions are collected by CBC than GLPK.\n(15 runs of Steiner forest routines with 100x bootstraps)",
       x = "Solver", y = "Number of solutions")

