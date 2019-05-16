#' R implementation of Steiner Tree solver
#'
#' A package dedicated to finding minimum Steiner trees in networks. Particularly biological ones,
#' which tend to be very sparse and on the easier end of the spectrum. This package is particularly
#' interested in Minimum Steiner Tree (MStT), Maximum-Weight Connected Subgraph (MWCS) and
#' Prize-Collecting Steiner Tree (PCST) problems.
#'
#' This package serves as a faithful implementation of "Thinning out Steiner Trees".
#'
#' @docType package
#' @name stoneTrees
#' @author Adam Sardar
#' @import R6
#' @import data.table
#' @import igraph
#' @import Matrix
#' @import Rglpk
NULL
