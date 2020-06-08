globalVariables(c("."))

stoneTrees_solvers <- c("cplexAPI", "Rglpk", "rcbc", "lpSolve" ,"lpsymphony")

# Choose the best solver from those available
chooseSolver <- function(){

  return( toupper(stoneTrees_solvers[stoneTrees_solvers %in% .packages(all.available = "TRUE")][1]) )
}

# This function uses the lower-level cplexAPI package, which does not suffer from the segfault bug in Rcplex for larger matricies
# The use of cplexAPI:: so as to acess the namespace is horrid, but necessary for a clean build process without CPLEX being available
#' @importFrom parallel detectCores
solver_CPLEXapi <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0), ...){

  if(!"cplexAPI" %in% .packages(all.available = TRUE) ) stop("cplexAPI must be installed in order to use the CPLEX solver")

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}

  if(!is.null(cplexParamList$trace)) cplexParamList$trace %<>% as.integer
  if(!is.null(cplexParamList$tilim)) cplexParamList$tilim %<>% as.numeric
  if(!is.null(cplexParamList$nThreads)) cplexParamList$nThreads %<>% as.integer

  Acpx <- toCPXMatrix(Amat)

  env <- cplexAPI::openEnvCPLEX()

  prob <- cplexAPI::initProbCPLEX(env, pname = "RankMatrixDecompostion")

  if(cplexParamList$trace > 0){ cplexAPI::setIntParmCPLEX(env, cplexAPI::CPX_PARAM_SCRIND, cplexAPI::CPX_ON) }
  if(!is.null(cplexParamList$tilim)){ cplexAPI::setDblParmCPLEX(env, cplexAPI::CPX_PARAM_TILIM, cplexParamList$tilim) }

  if(!is.null(cplexParamList$nThreads)){

    cplexAPI::setIntParmCPLEX(env, cplexAPI::CPX_PARAM_THREADS, cplexParamList$nThreads)
  }else{

    cplexAPI::setIntParmCPLEX(env, cplexAPI::CPX_PARAM_THREADS, as.integer( max(c(detectCores() - 4,1, na.rm = TRUE))) )
  }

  cplexAPI::copyLpCPLEX(env = env,
                        lp = prob,
                        nCols = ncol(Amat),
                        nRows = nrow(Amat),
                        lpdir = cplexAPI::CPX_MAX,
                        objf = unname(cVec),
                        rhs = bVec,
                        sense = unlist(list("<" = "l", "<=" = "L", "==" = "E", ">=" = "G", ">" = "g")[senseVec]),
                        matbeg = Acpx$matbeg,
                        matcnt = Acpx$matcnt,
                        matind = Acpx$matind,
                        matval = Acpx$matval,
                        lb = rep(0,length(cVec)) ,
                        ub = rep(1,length(cVec)),
                        rngval = NULL)

  cplexAPI::copyColTypeCPLEX(env = env,
                             lp = prob,
                             xctype = vtypeVec)

  cplexAPI::mipoptCPLEX(env, prob)

  MILPsolve <- cplexAPI::solutionCPLEX(env, prob)

  cplexAPI::delProbCPLEX(env, prob)

  cplexAPI::closeEnvCPLEX(env)

  MILPsolve$solution <- MILPsolve$x

  return(MILPsolve)
}

solver_GLPK <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0), ...){

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}

  if(!"Rglpk" %in% .packages(all.available = TRUE) ) stop("Rglpk must be installed in order to use the GLPK solver")

  glpkParamList <- list(verbose = FALSE, presolve = TRUE, tm_limit = 0, canonicalize_status = TRUE)
  glpkParamList$tm_limit <-  ifelse(!is.null(cplexParamList$tilim), as.integer(cplexParamList$tilim)*1000, as.integer(300*1000))
  glpkParamList$verbose <-  ifelse(cplexParamList$trace == 0, FALSE, TRUE)

  MILPsolve <- Rglpk::Rglpk_solve_LP(obj = cVec,
                                     mat = Amat,
                                     rhs = bVec,
                                     dir = senseVec,
                                     max = TRUE,
                                     types = vtypeVec,
                                     control = glpkParamList)

  return(MILPsolve)
}



#' @import lpSolve
solver_LPSOLVE <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0), nSols = 1){

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}

  if(!"lpSolve" %in% .packages(all.available = TRUE) ) stop("lpSolve must be installed in order to use the GLPK solver")

  validateSingleInteger(nSols)

  MILPsolve <- lpSolve::lp(objective.in = cVec,
                           const.mat = as.matrix(Amat),
                           const.rhs = bVec,
                           const.dir = senseVec,
                           all.bin = TRUE,
                           num.bin.solns = nSols,
                           direction = "max")

  return(MILPsolve)
}

solver_SYMPHONY <-  function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0), ...){

  if(! all( c("lpsymphony","slam") %in% .packages(all.available = TRUE)) ) stop("lpsymphony and slam must both be installed in order to use the GLPK solver")

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}

  MILPsolve <- lpsymphony::lpsymphony_solve_LP(obj = cVec,
                                               mat = slam::as.simple_triplet_matrix(Amat),
                                               rhs = bVec,
                                               dir = senseVec,
                                               max = TRUE,
                                               types = vtypeVec,
                                               time_limit = ifelse( is.numeric(cplexParamList$tilim), as.integer(cplexParamList$tilim), -1),
                                               verbosity = ifelse( is.numeric(cplexParamList$trace), as.integer(cplexParamList$trace) -2, -2))
  return(MILPsolve)
}



solver_CBC <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0), nSols = 1){
  
  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}
  
  if(!"rcbc" %in% .packages(all.available = TRUE) ) stop("rcbc must be installed in order to use the CBC solver")
  
  #Prepare the column bounds
  colUB <- rep.int(Inf, ncol(Amat))
  colLB <- rep.int(-Inf, ncol(Amat))
  
  colUB[vtypeVec == "B"] <- 1L
  colLB[vtypeVec == "B"] <- 0L
  
  # Prepare the row bounds
  rowUB <- rep.int(Inf, nrow(Amat))
  rowLB <- rep.int(-Inf, nrow(Amat))
  
  rowUB[senseVec == "=="] <- bVec[senseVec == "=="]
  rowLB[senseVec == "=="] <- bVec[senseVec == "=="]
  
  rowUB[senseVec == ">="] <- Inf
  rowLB[senseVec == ">="] <- bVec[senseVec == ">="]

  rowUB[senseVec == "<="] <- bVec[senseVec == "<="]
  rowLB[senseVec == "<="] <- -Inf
  
  MILPsolve <- rcbc::cbc_solve(
    obj = cVec,
    mat = Amat,
    
    row_ub = rowUB,
    row_lb = rowLB,
    
    col_ub = colUB,
    col_lb = colLB,
    
    is_integer = (vtypeVec == "B"),
    
    max = TRUE,
    cbc_args = list("sec" = ifelse( is.numeric(cplexParamList$tilim), as.integer(cplexParamList$tilim), -1),
                    "logLevel" = ifelse( is.numeric(cplexParamList$trace), as.integer(cplexParamList$trace), 0) )
  )
  
  MILPsolve$solution <- MILPsolve$column_solution
  
  return(MILPsolve)
}
