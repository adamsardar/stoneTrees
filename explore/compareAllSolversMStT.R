devtools::load_all()
library(igraph)
library(purrr)
library(dplyr)
library(data.table)
library(ggplot2)




# reformat node score for MStT (not MWCS) and add Terminal attributes
fixedTerminalLymphomaGraph <- lymphomaGraph
V(fixedTerminalLymphomaGraph)$isTerminal <- FALSE
V(fixedTerminalLymphomaGraph)[nodeScore > 0]$isTerminal <- TRUE
V(fixedTerminalLymphomaGraph)$nodeScore <- -1


set.seed(2345)

#set the number of trials
numberOfTrial <- 3
nsuboptimalSol <- 500


MStTSuboptimalBenchDT <- map(1:numberOfTrial, function(i){

  print(i)

  # Note, The two process are not solving the exact same problem. Stochastic algorithms are a pain

Results <- map(stoneTrees_solvers[1:3], safely(~{ #lpsymphony is creating issues

  message(paste("Starting", .x))

    time <- system.time( SteinForSolverX <- subOptimalSteinerProblem$new(fixedTerminalLymphomaGraph,
                                                                          solverChoice = .x,
                                                                          verbose = TRUE,
                                                                         solverTrace = 0)$identifyMultipleSteinerSolutions(maxItr = nsuboptimalSol) )

    SolutionPool <- SteinForSolverX$getSolutionPoolGraphs(collapseSols = FALSE)

    ScorePool <- SteinForSolverX$getSolutionPoolScores()

    save(SolutionPool, file = paste0("./results/MStTSolver_",.x,"_Trial", i,".RData"))

     Sizes <- map_int(SolutionPool, vcount)

  return(data.table(solver = paste(.x),
                    time = time["elapsed"],
                    vcount = vcount(SteinForSolverX$getSolutionPoolGraphs()),
                    modulesVcount = list(Sizes[-length(Sizes)]),
                    Niteration = length(Sizes[-length(Sizes)]),
                    modulesScore = list(ScorePool)))

   message(paste("Completed", .x))

}, otherwise = NA_real_))

ResultsDT <- unlist(Results, recursive = FALSE) %>%
             .[map_lgl(., ~{is.data.table(.x)})]%>%
             rbindlist()

return(ResultsDT[,trial := i])

}) %>% rbindlist()

fwrite(MStTSuboptimalBenchDT, "./results/MStTSuboptimalBench.tsv")






MStTSuboptimalBenchDT <- fread("./results/MStTSuboptimalBench.tsv")%>%
  .[, `:=`(modulesVcount = strsplit(modulesVcount, "\\|"),
           modulesScore = strsplit(modulesScore, "\\|"))]


MStTSuboptimalBenchDT[trial == 1, .(solver, time, vcount, NSolutions = Niteration)]

MStTSuboptimalBenchDT[, .(`TotalTime(hours)` = sum(time)/60/60), by = .(solver)]

MStTSuboptimalBenchDT[, .(TotalTime = sum(time)/60/60), by = .(solver)][,TotalTime] %>% sum

MStTSuboptimalBenchDT[, .(AvSise = mean(vcount)), by = .(solver)]

#
p1 <- ggplot(MStTSuboptimalBenchDT, aes(x = solver, y = time/60)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, aes(color = as.factor(trial))) +
  #scale_y_continuous( breaks = seq(0,100,5)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")+
  labs(y = "Time (min)",
       title = "Time taken to solve Steiner forest problem")

p2 <- ggplot(MStTSuboptimalBenchDT, aes(x = solver, y = time/60)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, aes(color = as.factor(trial))) +
  #scale_y_continuous( breaks = seq(0,100,5)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")+
  scale_y_log10()+
  labs(y = "Time (logScale)",
       title = "Time taken to solve Steiner forest problem")


#
p3 <- ggplot(MStTSuboptimalBenchDT, aes(x = solver, y = vcount))+
  geom_boxplot() +
  geom_jitter(width = 0.1, aes(color = as.factor(trial))) +
  theme_bw(base_size = 20) +
  #scale_color_hue()+
  theme(legend.position = "none")+
  labs( y = "Aggregated module size",
        title = "Size of Steiner forest solution")


p4 <- MStTSuboptimalBenchDT[,.(modulesVcount = unlist(modulesVcount)), by = .(solver, time, vcount, Niteration, trial)] %>%

  ggplot(., aes(x = solver, y = as.numeric(modulesVcount))) +
  geom_jitter(width = 0.1, aes(color = as.factor(trial)), alpha = 0.2) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")+
  labs(y = "vCount of all modules in the pool",
       title = "Individual tree solutions to resampled terminals")


p5 <- ggplot(MStTSuboptimalBenchDT, aes(x = solver, y = Niteration)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, aes(color = as.factor(trial))) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")+
  labs(y = "Number of solutions",
       title = "Size of solution pool")


plot <- egg::ggarrange(p1,p2,p3,p4,p5, top =  paste(numberOfTrial , "collecting up to",  nsuboptimalSol, " solutions on each run."))



ggsave(plot = plot, filename = "./results/CBCvsCPLEXvcGPLKPerformanceMStT.png", width = 16, height = 15)
