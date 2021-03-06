#Useful validation routines
globalVariables(c("."))

#' @importFrom  ensurer ensure
#' @import igraph
validateIsNetwork = function(network2validate, singleWeakComponent = TRUE, isDirected = NULL){

  network2validate %<>% ensure(is.igraph,
                               err_desc = "Input network must be an igraph object")

  if(singleWeakComponent){

    network2validate %<>% ensure(length(decompose(., mode = "weak")) == 1,
                                 err_desc = "Input network must a single connected component (consider using igraph::decompose?)")
  }

  if(!is.null(isDirected)){

    validateFlag(isDirected)

    network2validate %>% ensure(is.directed(.) == isDirected,
                                err_desc = paste0("Expecting is.directed() output of graph to be ",isDirected) )
  }

  if("nodeScore" %in% vertex_attr_names(network2validate)){

    V(network2validate)$nodeScore %>% ensure(is.numeric,
                                             !any(is.na(.)),
                                             err_desc = "nodeScore values for input graph must be numeric vectors with no NA entries")
  }

  if("isTerminal" %in% vertex_attr_names(network2validate)){

    V(network2validate)$isTerminal %>% ensure(is.logical,
                                             all(!is.na(.)),
                                             err_desc = "isTerminal values for input graph must be logical vectors without NA entries")
  }

  invisible(network2validate)
}

#' @importFrom  ensurer ensure
validateSquareMatrix = function(Matrix2validate){

  Matrix2validate %<>% ensure(toupper(class(.)) %like% "MATRIX",
                              nrow(.) > 0,
                              nrow(.) == ncol(Matrix2validate),
                              is.numeric(.[1,1]),
                              err_desc = "Candidate matrix must be a numeric square matrix or Matrix-derived class")

  invisible(Matrix2validate)
}

#' @importFrom  ensurer ensure
validateSingleInteger = function(int2evaluate){

  int2evaluate %<>% ensure(is.numeric,
                      is.whole.number,
                      length(.) == 1,
                      err_desc = "Expecting a single numeric integer!")

  invisible(int2evaluate)
}

#' @importFrom  ensurer ensure
validateSinglePositiveSemiDefiniteNumeric = function(num2evaluate){

  num2evaluate %<>% ensure(is.numeric,
                           length(.) == 1,
                           . >= 0,
                           err_desc = "Expecting a single positive semi-definite numeric value!")

  invisible(num2evaluate)
}

#' @importFrom ensurer ensure
#' @importFrom stringr str_c
validateFlag = function(flag2validate){

  flagName = deparse(substitute(flag2validate))

  flag2validate %<>% ensure(is.logical(.),
                          length(.) == 1,
                          err_desc = str_c(flagName," must be a single boolean"))

  invisible(flag2validate)
}

#' @importFrom  ensurer ensure
validateSolverChoice = function(candidateSolverChoice){

  candidateSolverChoice %<>%  ensure(is.character,
                                     toupper(.) %in% toupper(stoneTrees_solvers),
                                     err_desc = paste0( c("At current the supported solvers are:", stoneTrees_solvers), collapse = " "))

  invisible(toupper(candidateSolverChoice))
}

is.whole.number = function(x){ x == as.integer(x)}
