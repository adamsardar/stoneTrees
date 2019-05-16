#' Function to map a constraint matrix to a cplex compatible row,column, value set of triples
#'
#' This function is copied wholesale from Rcplex: See Rcplex:::toCPXMatrix. However, Rcplex is a difficult package
#' to include as a dependency as it requires explicit compilation with links to proprietry header/library files.
#'
#' @param Amat Constraint matrix. Can be of type sparseMatrix or simple_triplet_matrix
#'
#' @return A list compatible with the copyLpCPLEX in cplexAPI
#' @importFrom methods as is
#' @seealso cplexAPI::copyLpCPLEX
toCPXMatrix <- function (Amat) {
  
  if (is.null(Amat)) {  return(NULL)}
  
  else if (is(Amat, "sparseMatrix")) {
    
    Amat <- as(Amat, "dgCMatrix") # Must be column major
    matbeg <- Amat@p
    matcnt <- diff(c(Amat@p, length(Amat@x)))
    matind <- Amat@i
    matval <- Amat@x
  }else if (is(Amat, "simple_triplet_matrix")) {
    
    matbeg <- c(0L, cumsum(tabulate(Amat$j, Amat$ncol)))
    matcnt <- tabulate(Amat$j, Amat$ncol)
    matind <- Amat$i - 1L
    matval <- Amat$v
  }  else {

    matbeg <- (0L:(ncol(Amat) - 1L)) * nrow(Amat)
    matcnt <- rep(nrow(Amat), ncol(Amat))
    matind <- rep(0L:(nrow(Amat) - 1L), ncol(Amat))
    matval <- as.vector(Amat)
  }
  
  return(list(matbeg = as.integer(matbeg),
              matcnt = as.integer(matcnt), 
              matind = as.integer(matind),
              matval = as.double(matval)))
}