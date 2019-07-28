# stoneTrees - an R package for solving Stiener-tree problems

A package dedicated to finding minimum Steiner trees in networks. Particularly biological networks; which tend to be very sparse and on the easier end of the spectrum of problems. This package is particularly interested in Minimum Steiner Tree (MStT) and Maximum-Weight Connected Subgraph (MWCS) problems.

This package serves as a faithful implementation of "Thinning out Steiner Trees" (with a few bells and whistles added on the sides), which is the algorithm that more-or-less won the [DIMACS11 competition on algorithms for solving Steiner Tree problems](http://dimacs11.zib.de/).

## Installation

```
> devtools::install_github("adamsardar/stoneTrees")
> library(stoneTrees)
```

The default solver in the package is `lpSolve`. However, empirical observation leads to the suggestion of using [Rglpk](https://cran.r-project.org/web/packages/Rglpk/index.html) or, better yet, [Rcplex](https://cran.r-project.org/web/packages/Rcplex/index.html). The installation of these last two packages, whilst relatively simple, are sufficiently complex on Windows to affect the solver reccomendation.

`Rglpk` can be easily installed. On linux search for the the glpk-dev package (`apt install libglpk-dev` for Debian flavoured distros); on mac you can use port (`port install glpk`) and on Windows you can [follow the community install guild](http://winglpk.sourceforge.net/). Following that `> install.packages("Rglpk")` should proceed smoothly.

## Usage

Problems are constructed using the $new method of the appropriate instatialisation and then a collector method is called. For the base Steiner tree or MWCS problem (or a blend of the two), this is nodeCentricSteinerTreeProblem$new() followed by $findSingleSteinerSolution():

Using a test dataset that comes with the package: `lymphomaGraph`, which is a MWCS problem. It has a node score attribute detailing prizes/costs of node inclusion:

```
library(igraph)
library(ggplot2)
qplot(V(lymphomaGraph)$nodeScore)
```

Most nodes have negative weights (costs) - the MWCS looks to group as many positive weights (prizes) together.

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

If the user is interested in collecting sub-optimal solutions, then a difference constructor is used that allows one to specicy the solution tolerence that is acceptable. Notice here that there are three stages: build the object (`subOptimalSteinerProblem$new()`), identify solutions (`lymphoma_multiMWCS$identifyMultipleSteinerSolutions()`) and finally extract solutions from the pool as graphs (`lymphoma_multiMWCS$getSolutionPoolGraphs()`).

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

Finally, if the user is interested in the bootstrapped Steiner Tree problem (aka the Steiner Forest problem), then a third class is used. This is a seed/terminal-based routine, so a graph with terminals must be used: one included in the package is the `karateGraph` problem. The methods used are: `nodeCentricSteinerForestProblem$new()`, `$sampleMultipleBootstrapSteinerSolutions()` which populated a solution pool with random draws from the Steiner forest problem and `$getBootstrapSolutionPoolGraphs()` which returns the aggregated result.

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
