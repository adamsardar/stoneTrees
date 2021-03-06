% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subOptimalNodeCentricSteinerTreeProblem-class.R
\docType{class}
\name{subOptimalSteinerProblem}
\alias{subOptimalSteinerProblem}
\title{Collect degenerate and sub-optimal solutions to Steiner problems (MStTP or MWCS) with uniform or no edge weights.}
\format{
R6Class \code{subOptimalSteinerProblem} Construct an object representation of a multiple-solution Steiner tree/maximum weight connected subgraph (MWCS) problem
}
\description{
Rather than find just a single solution to a MStTP/MWCS, one can populate a solution pool with multiple degenerate/tolerable solutions.
}
\details{
This class is derived from *nodeCentricSteinerTreeProblem*: all methods available in the superclass are available here. The difference is that after
each acceptable solution is found, the solution is a.) stored in a solution pool and b.) used to generate a 'novelty' constraint on future solutions.
}
\section{methods}{


Alongisde those for *nodeCentricSteinerTreeProblem*
\describe{
   \item{\code{new(network, solverChoice = chooseSolver(), verbose = TRUE, presolveGraph = TRUE, solverTimeLimit = 300, solverTrace = as.integer(verbose), solutionTolerance = 0)}}{Constructor for the subOptimalSteinerProblem class. Alongside the arguments for the super-class constructor, there is also 'solutionTolerance', which instructs the object as to the gap between optimal and observed solution that is acceptable.}
   \item{\code{identifyMultipleSteinerSolutions(maxItr = 10)}}{Add solutions to the solution pool. maxItr is an argument dictating the number of runs through the optimsation procedure.}
   \item{\code{getSolutionPoolGraphs(collapseSols = TRUE)}}{Either return a list of solutions within tolerance (collapseSols = FALSE) or pool all solutions together and return a single graph (collapseSols = TRUE, defaults)}
   \item{\code{getSolutionPoolScores()}}{Compute the scores of the solutions in the solution pool. These are in the same order as the list of graphs returned by $getSolutionPoolGraphs(FALSE)}
   \item{\code{getOptimumScore()}}{Returns the optimum score from solutions in the solution pool}
   \item{\code{getSolutionTolerance()}}{Retreive the tolerance that permits a solution to be added to the solution pool in future calls to $identifyMultipleSteinerSolutions()}
   \item{\code{setSolutionTolerance(x)}}{Alter the tolerance that permits a solution to be added to the solution pool in future calls to $identifyMultipleSteinerSolutions()}
}
}

\examples{
library(igraph)

 # Maximum-Weight Connected Subgraph (MWCS) - find sub-optimal solutions

 ## Vertex attribute details node costs/prizes
 head(V(lymphomaGraph)$nodeScore)

 lymphoma_multiMWCS = subOptimalSteinerProblem$new(lymphomaGraph, solutionTolerance = 0.5)

 #Populate the solution pool with multiple solutions - notice the
 lymphoma_multiMWCS$identifyMultipleSteinerSolutions()

 lymphoma_multiMWCS$getSolutionPoolGraphs(collapseSols = FALSE)

 lymphoma_multiMWCS$getSolutionPoolScores()

 #All solution scores are within tolerance
 diff(range(lymphoma_multiMWCS$getSolutionPoolScores()))
}
\references{
Fischetti M, Leitner M, Ljubić I, Luipersbeck M, Monaci M, Resch M, et al. Thinning out Steiner trees: a node-based model for uniform edge costs. Math Program Comput. dimacs11.cs.princeton.edu; 2017;9: 203–229.

Beisser D, Klau GW, Dandekar T, Müller T, Dittrich MT. BioNet: An R-Package for the functional analysis of biological networks. Bioinformatics. 2010;26: 1129–1130.

\url{https://en.wikipedia.org/wiki/Steiner_tree_problem}
}
\seealso{
nodeCentricSteinerTreeProblem

Other SteinerProblemSolver: 
\code{\link{nodeCentricSteinerForestProblem}},
\code{\link{nodeCentricSteinerTreeProblem}}
}
\concept{SteinerProblemSolver}
\section{Super class}{
\code{\link[stoneTrees:nodeCentricSteinerTreeProblem]{stoneTrees::nodeCentricSteinerTreeProblem}} -> \code{subOptimalSteinerProblem}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{subOptimalSteinerProblem$new()}}
\item \href{#method-getSolutionPool}{\code{subOptimalSteinerProblem$getSolutionPool()}}
\item \href{#method-getSolutionPoolGraphs}{\code{subOptimalSteinerProblem$getSolutionPoolGraphs()}}
\item \href{#method-getSolutionPoolScores}{\code{subOptimalSteinerProblem$getSolutionPoolScores()}}
\item \href{#method-getOptimumScore}{\code{subOptimalSteinerProblem$getOptimumScore()}}
\item \href{#method-getNoveltyConstraints}{\code{subOptimalSteinerProblem$getNoveltyConstraints()}}
\item \href{#method-getSolutionTolerance}{\code{subOptimalSteinerProblem$getSolutionTolerance()}}
\item \href{#method-setSolutionTolerance}{\code{subOptimalSteinerProblem$setSolutionTolerance()}}
\item \href{#method-getNconnectivityConstraintsCalls}{\code{subOptimalSteinerProblem$getNconnectivityConstraintsCalls()}}
\item \href{#method-identifyMultipleSteinerSolutions}{\code{subOptimalSteinerProblem$identifyMultipleSteinerSolutions()}}
\item \href{#method-clone}{\code{subOptimalSteinerProblem$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="findSingleSteinerSolution">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-findSingleSteinerSolution}{\code{stoneTrees::nodeCentricSteinerTreeProblem$findSingleSteinerSolution()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getConnectivityConstraints">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getConnectivityConstraints}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getConnectivityConstraints()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getCurrentSolutionGraph">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getCurrentSolutionGraph}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getCurrentSolutionGraph()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getCurrentSolutionScore">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getCurrentSolutionScore}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getCurrentSolutionScore()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getEdgeDT">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getEdgeDT}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getEdgeDT()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getFixedTerminalConstraints">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getFixedTerminalConstraints}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getFixedTerminalConstraints()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getNodeDT">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getNodeDT}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getNodeDT()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getNodeDegreeConstraints">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getNodeDegreeConstraints}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getNodeDegreeConstraints()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getTerminals">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getTerminals}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getTerminals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getTwoCycleConstraints">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getTwoCycleConstraints}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getTwoCycleConstraints()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="getnConnectivityConstraintsCalls">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-getnConnectivityConstraintsCalls}{\code{stoneTrees::nodeCentricSteinerTreeProblem$getnConnectivityConstraintsCalls()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="stoneTrees" data-topic="nodeCentricSteinerTreeProblem" data-id="isSolutionConnected">}\href{../../stoneTrees/html/nodeCentricSteinerTreeProblem.html#method-isSolutionConnected}{\code{stoneTrees::nodeCentricSteinerTreeProblem$isSolutionConnected()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$new(
  network,
  solverChoice = chooseSolver(),
  verbose = TRUE,
  presolveGraph = TRUE,
  solverTimeLimit = 300,
  solutionTolerance = 0,
  solverTrace = as.integer(verbose)
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolutionPool"></a>}}
\if{latex}{\out{\hypertarget{method-getSolutionPool}{}}}
\subsection{Method \code{getSolutionPool()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getSolutionPool()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolutionPoolGraphs"></a>}}
\if{latex}{\out{\hypertarget{method-getSolutionPoolGraphs}{}}}
\subsection{Method \code{getSolutionPoolGraphs()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getSolutionPoolGraphs(collapseSols = TRUE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolutionPoolScores"></a>}}
\if{latex}{\out{\hypertarget{method-getSolutionPoolScores}{}}}
\subsection{Method \code{getSolutionPoolScores()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getSolutionPoolScores()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getOptimumScore"></a>}}
\if{latex}{\out{\hypertarget{method-getOptimumScore}{}}}
\subsection{Method \code{getOptimumScore()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getOptimumScore()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getNoveltyConstraints"></a>}}
\if{latex}{\out{\hypertarget{method-getNoveltyConstraints}{}}}
\subsection{Method \code{getNoveltyConstraints()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getNoveltyConstraints()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolutionTolerance"></a>}}
\if{latex}{\out{\hypertarget{method-getSolutionTolerance}{}}}
\subsection{Method \code{getSolutionTolerance()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getSolutionTolerance()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setSolutionTolerance"></a>}}
\if{latex}{\out{\hypertarget{method-setSolutionTolerance}{}}}
\subsection{Method \code{setSolutionTolerance()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$setSolutionTolerance(x)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getNconnectivityConstraintsCalls"></a>}}
\if{latex}{\out{\hypertarget{method-getNconnectivityConstraintsCalls}{}}}
\subsection{Method \code{getNconnectivityConstraintsCalls()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$getNconnectivityConstraintsCalls()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-identifyMultipleSteinerSolutions"></a>}}
\if{latex}{\out{\hypertarget{method-identifyMultipleSteinerSolutions}{}}}
\subsection{Method \code{identifyMultipleSteinerSolutions()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$identifyMultipleSteinerSolutions(maxItr = 10)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{subOptimalSteinerProblem$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
