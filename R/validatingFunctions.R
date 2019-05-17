#Useful validation routines
globalVariables(c("."))

#' @importFrom  ensurer ensure
validateIsNetwork <- function(network2validate){

  network2validate %<>% ensure(is.igraph,
                               err_desc = "Input network must be an igraph object")

  network2validate %<>% ensure(length(decompose(., mode = "weak")) == 1,
                               err_desc = "Input network must a single connected component (consider using igraph::decompose?)")

  invisible(network2validate)
}

#' @importFrom  ensurer ensure
validateSquareMatrix <- function(Matrix2validate){

  Matrix2validate %<>% ensure(toupper(class(.)) %like% "MATRIX",
                              nrow(.) > 0,
                              nrow(.) == ncol(Matrix2validate),
                              is.numeric(.[1,1]),
                              err_desc = "Candidate matrix must be a numeric square matrix or Matrix-derived class")

  invisible(Matrix2validate)
}

#' @importFrom  ensurer ensure
validateSingleInteger <- function(int2evaluate){

  int2evaluate %<>% ensure(is.numeric,
                      is.whole.number,
                      length(.) == 1,
                      err_desc = "Expecting a single numeric integer!")

  invisible(int2evaluate)
}


#' @importFrom ensurer ensure
#' @importFrom stringr str_c
validateFlag <- function(flag2validate){

  flagName <- deparse(substitute(flag2validate))

  flag2validate %<>% ensure(is.logical(.),
                          length(.) == 1,
                          err_desc = str_c(flagName," must be a single boolean"))

  invisible(flag2validate)
}

#' @importFrom  ensurer ensure
solverChoiceValidator <- function(candidateSolverChoice){

  candidateSolverChoice %<>%  ensure(is.character,
                                     toupper(.) %in% c("RCPLEX","CPLEXAPI","RGLPK","LPSYMPHONY","LPSOLVE"),
                                     err_desc = "At current only RCPLEX, CPLEXAPI, RGLPK, LPSYMPHONY and LPSOLVE solvers are supported")

  invisible(toupper(candidateSolverChoice))
}

is.whole.number <- function(x){ x == as.integer(x)}
