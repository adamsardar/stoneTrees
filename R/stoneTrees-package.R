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
#' @import igraph
#' @import Matrix
#' @import data.table
#' @importFrom magrittr "%<>%"
#' @importFrom rlang is_vector
#' @importFrom rlang is_list
#' @importFrom rlang is_na
#' @importFrom rlang abort
#' @importFrom rlang ffi_standalone_is_bool_1.0.7
#' @importFrom rlang is_string
#' @importFrom rlang ffi_standalone_check_number_1.0.7
#' @importFrom rlang  env_get_list
#' @importFrom rlang `%||%`
#' @importFrom rlang is_missing
#' @importFrom rlang is_logical
#' @importFrom rlang caller_arg
#' @importFrom rlang caller_env
NULL
