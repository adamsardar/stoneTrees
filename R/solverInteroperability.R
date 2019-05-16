globalVariables(c("."))

# Note that Rcplex can segfault when solving for large decomposition problems
#' @importFrom utils installed.packages
solver_CPLEX <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0)){

  if(!"Rcplex" %in% rownames(installed.packages())) stop("Rcplex must be installed in order to use the CPLEX solver")

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}

  if(!is.null(cplexParamList$trace)) cplexParamList$trace %<>% as.integer
  if(!is.null(cplexParamList$tilim)) cplexParamList$tilim %<>% as.integer

  MILPsolve <- Rcplex::Rcplex(cvec = cVec,
                        Amat = Amat,
                        bvec = bVec,
                        sense = senseVec,
                        objsense = "max",
                        vtype = vtypeVec,
                        control = cplexParamList,
                        n=1)

  MILPsolve$solution <- MILPsolve$xopt

   return(MILPsolve)
}

# This function uses the lower-level cplexAPI package, which does not suffer from the segfault bug in Rcplex for larger matricies
# The use of cplexAPI:: so as to acess the namespace is horrid, but necessary for a clean build process without CPLEX being available
#' @importFrom utils installed.packages
#' @importFrom parallel detectCores
solver_CPLEXapi <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0)){
  
  if(!"cplexAPI" %in% rownames(installed.packages())) stop("cplexAPI must be installed in order to use the CPLEX solver")
  
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
                        sense = senseVec,
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

#' @importFrom utils installed.packages
solver_GLPK <- function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0)){

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}

  if(!"Rglpk" %in% rownames(installed.packages())) stop("Rglpk must be installed in order to use the GLPK solver")

  glpkParamList <- list(verbose = FALSE, presolve = TRUE, tm_limit = 0, canonicalize_status = TRUE)
  glpkParamList$tm_limit <-  ifelse(!is.null(cplexParamList$tilim), as.integer(cplexParamList$tilim)*1000, as.integer(300*1000))
  glpkParamList$verbose <-  ifelse(cplexParamList$trace == 0, FALSE, TRUE)

  MILPsolve <- Rglpk::Rglpk_solve_LP(obj = cVec,
                               mat = Amat,
                               rhs = bVec,
                               dir = unlist(list(l = "<", L = "<=", E = "==", G = ">=", g = ">")[senseVec]),
                               max = TRUE,
                               types = vtypeVec,
                               control = glpkParamList)
  
  return(MILPsolve)
}

#' @import lpsymphony
#' @importFrom slam as.simple_triplet_matrix
solver_SYMPHONY <-  function(cVec, Amat, senseVec, bVec=0, vtypeVec="B", cplexParamList = list(trace = 0)){

  if(length(bVec) == 1){bVec %<>% rep(nrow(Amat))}
  if(length(senseVec) == 1){senseVec %<>% rep(nrow(Amat))}
  if(length(vtypeVec) == 1){vtypeVec %<>% rep(ncol(Amat))}

  MILPsolve <- lpsymphony_solve_LP(obj = cVec,
                             mat = as.simple_triplet_matrix(Amat),
                             rhs = bVec,
                             dir = unlist(list(l = "<", L = "<=", E = "==", G = ">=", g = ">")[senseVec]),
                             max = TRUE,
                             types = vtypeVec,
                             time_limit = ifelse( is.numeric(cplexParamList$tilim), as.integer(cplexParamList$tilim), -1),
                             verbosity = ifelse( is.numeric(cplexParamList$trace), as.integer(cplexParamList$trace) -2, -2))
  return(MILPsolve)
}
