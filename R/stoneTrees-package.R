#' Solve Minimum Steiner Tree (MStT) and Maximum-Weight Connected Subgraph Problems Using Integer Linear Programming
#'
#' A package dedicated to finding minimum Steiner trees in networks. Particularly biological ones,
#' which tend to be very sparse and on the easier end of the spectrum. This package is especially
#' interested in Minimum Steiner Tree (MStT) and Maximum-Weight Connected Subgraph (MWCS) problems.
#'
#' This package serves as a faithful implementation of "Thinning out Steiner Trees" (with a few bells and whistles added on the sides).
#'
#' @docType package
#' @name stoneTrees
#' @author Adam Sardar
#' @references Fischetti M, Leitner M, Ljubić I, Luipersbeck M, Monaci M, Resch M, et al. Thinning out Steiner trees: a node-based model for uniform edge costs. Math Program Comput. dimacs11.cs.princeton.edu; 2017
#' @import R6
#' @importFrom data.table ".__T__$:base" ".__T__$<-:base" ".__T__[[<-:base" ".__T__[<-:base" ".__T__[:base"  rbindlist as.data.table data.table ":=" %like% fread setkey setnames
#' @import igraph
#' @import Matrix
#' @importFrom magrittr "%<>%"
NULL
