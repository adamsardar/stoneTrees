  <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/adamsardar/stoneTrees.svg?branch=master)](https://travis-ci.org/adamsardar/stoneTrees)
  [![Codecov test coverage](https://codecov.io/gh/adamsardar/stoneTrees/branch/master/graph/badge.svg)](https://codecov.io/gh/adamsardar/stoneTrees?branch=master)
  [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/adamsardar/stoneTrees?branch=master&svg=true)](https://ci.appveyor.com/project/adamsardar/stoneTrees)
  <!-- badges: end -->
  
# stoneTrees - an R package for solving Steiner tree problems

A package dedicated to finding solutions to [Steiner tree problems](https://en.wikipedia.org/wiki/Steiner_tree_problem) in graphs using Integer Linear Programming (ILP). Motivation stems from a need for solutions to the Minimum Steiner Tree (MStT) and Maximum-Weight Connected Sub-graph (MWCS) problems in computational biology. For example:

* (Dittrich et al.)[https://www.ncbi.nlm.nih.gov/pubmed/18586718] (2008) demonstrate how MWCS is a means by which to combine per-gene expression data with mechanistic protein-protein interaction networks and extract functional modules in a data-driven way. 
* (Liang et al.)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5723949/] (2017) use MStT to infer phylogenies of mutated cancer cells from a dataset of copy-number variation (CNV).
* (Lam, Alexandersson and Pachter)[https://www.liebertpub.com/doi/abs/10.1089/10665270360688156] show that sequence alignment can be mapped to a Steiner tree problem.

This package serves as a faithful implementation of "Thinning out Steiner Trees" (with a few bells and whistles added on the sides) by Fischetti et al. (2017). This algorithm that more-or-less won the [DIMACS11 competition on algorithms for solving Steiner Tree problems](http://dimacs11.zib.de/). A major advancement of Fischetti et al. is that ILP variables are only intriduced for nodes, rather than edges *as well*, which dramatically decreases the run-time of the process.

## Installation

```
> devtools::install_github("adamsardar/stoneTrees") 
> library(stoneTrees)
```

### ILP Solvers

`stoneTrees` aims to be compatible with several ILP solvers; at current a user can chose from [`lpsymphony`](https://www.bioconductor.org/packages/release/bioc/html/lpsymphony.html), [`Rglpk`](https://cran.r-project.org/web/packages/Rglpk/index.html), [`lpSolve`](https://cran.r-project.org/web/packages/lpSolve/index.html) and [`Rcplex`](https://cran.r-project.org/web/packages/Rcplex/index.html) .  The default solver is `lpSolve`. However, there is a strong recommendation of the open-source `Rglpk` solver or, better yet, the proprietary `Rcplex` solver. Installation of these last two packages, whilst relatively straightforward, is complex enough to affect the choice of default.

`Rglpk` can be easily installed. On Linux install the glpk-dev package (`apt install libglpk-dev`); on Mac OSX you can use brew (`brew install glpk`) and on Windows you can [follow the community install guild](http://winglpk.sourceforge.net/). Following that `install.packages("Rglpk")` should work.

## Usage

Problems are constructed using the $new method of the appropriate problem class and then a collector method is called. For the base MStT or MWCS problem, this involves calling `nodeCentricSteinerTreeProblem$new()` followed by `$findSingleSteinerSolution()`.

By way of example, look at the `lymphoma` test dataset that comes with the package. It is a MWCS problem; a graph with node score attributes detailing prizes/costs of node inclusion:

```
library(igraph)
library(ggplot2)
qplot(V(lymphomaGraph)$nodeScore)
```

Most nodes have negative weights (costs) - the MWCS looks to group as many positive weights (prizes) together in a connected sub-graph.

```
>  lymphomaMWCS <- nodeCentricSteinerTreeProblem$new(lymphomaGraph)
>  lymphomaMWCS$findSingleSteinerSolution()
IGRAPH 6674598 UN-- 46 50 -- 
+ attr: nodeNameSep (g/c), SearchNetwork (g/c), name (v/c), nodeScore (v/n)
+ edges from 6674598 (vertex names):
 [1] 57 --58   58 --59   58 --96   96 --98   96 --106  90 --143  137--143  143--155  123--180  59 --380  28 --429 
[12] 143--473  180--491  490--501  96 --543  143--543  180--543  542--551  551--556  490--599  98 --608  62 --615 
[23] 143--650  551--650  28 --675  143--675  528--696  4  --808  490--808  528--808  543--808  696--808  4  --814 
[34] 63 --814  528--814  543--814  696--814  28 --927  143--931  98 --1059 615--1059 28 --1155 927--1155 556--896 
[45] 96 --1797 573--1797 96 --1720 543--962  501--1619 675--1879
```

If the user is interested in collecting sub-optimal solutions, then a different constructor is used that allows one to specify the solution tolerance parameter.

```
> lymphoma_multiMWCS <- subOptimalSteinerProblem$new(lymphomaGraph, solutionTolerance = 1)
> lymphoma_multiMWCS$identifyMultipleSteinerSolutions()
> lymphoma_multiMWCS$getSolutionPoolGraphs()
IGRAPH f0365e4 UN-- 57 84 -- 
+ attr: nodeNameSep (g/c), SearchNetwork (g/c), name (v/c), nodeScore (v/n)
+ edges from f0365e4 (vertex names):
 [1] 57 --58   58 --59   58 --96   96 --98   96 --106  90 --143  137--143  143--155  123--180  63 --373  90 --373 
[12] 320--373  59 --380  28 --429  143--473  180--491  373--491  314--494  314--501  490--501  62 --512  380--512 
[23] 429--512  473--512  314--528  96 --543  143--543  180--543  314--543  512--551  542--551  512--556  551--556 
[34] 320--599  490--599  98 --608  62 --615  143--650  551--650  28 --675  143--675  320--696  373--696  528--696 
[45] 528--697  696--697  494--757  491--759  494--759  528--759  543--759  697--759  512--766  4  --808  373--808 
[56] 490--808  494--808  528--808  543--808  696--808  4  --814  63 --814  314--814  373--814  528--814  543--814 
[67] 696--814  759--814  28 --927  143--931  98 --1059 615--1059 28 --1155 927--1155 556--896  96 --1797 573--1797
[78] 96 --1720 543--962  494--1619 501--1619 512--1769 675--1879 757--1987
```

Notice here that there are three stages: build the object (`subOptimalSteinerProblem$new()`), identify solutions and add them to the pool of distinct solutions (`lymphoma_multiMWCS$identifyMultipleSteinerSolutions()`) and finally extract solutions from the pool as graphs (`lymphoma_multiMWCS$getSolutionPoolGraphs()`).

If the user is interested in the bootstrapped Steiner Tree problem (aka the Steiner Forest problem), whereby initial seeds are repeatedly sub-sampled and the resultant Steiner problem solved, then a third class is used. This is a seed/terminal-based routine, so a graph with terminals must be used: one included in the package is the `karateGraph` problem. The methods used are: `nodeCentricSteinerForestProblem$new()`, `$sampleMultipleBootstrapSteinerSolutions()` which populates a solution pool with random draws from the Steiner forest problem and `$getBootstrapSolutionPoolGraphs()` which returns the aggregated result.

```
> nodeCentricSteinerForestProblem$new(karateGraph)$sampleMultipleBootstrapSteinerSolutions()$getBootstrapSolutionPoolGraphs()
IGRAPH 53c01e2 UN-- 5 7 -- Zachary
+ attr: name (g/c), SearchNetwork (g/c), name (v/c), isTerminal (v/l), .nodeID (v/n)
+ edges from 53c01e2 (vertex names):
[1] o--G w--G E--G o--H w--H E--H G--H
```

Notice the chaining of methods together.

# References:

D. Beisser, G. W. Klau, T. Dandekar, T. Mueller and M. Dittrich (2010) BioNet: an R-package for the Functional Analysis of Biological Networks. Bioinformatics.

Fischetti M, Leitner M, Ljubić I, Luipersbeck M, Monaci M, Resch M, et al. Thinning out Steiner trees: a node-based model for uniform edge costs. Math Program Comput. dimacs11.cs.princeton.edu; 2017

Dittrich MT, Klau GW, Rosenwald A, Dandekar T, Müller T. Identifying functional modules in protein-protein interaction networks: An integrated exact approach. Bioinformatics. 2008
