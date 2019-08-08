#' Wrapper function, the optimization itself is using the \code{cmaes::cma_es()} function.
#' Use this function as a template to create a \code{phaseSolver} with single individual.
#'
#' @title Phase optimization using (1+1)-CMA-ES
#' @param contextVector The current best solution
#' @param fun function to be optimized
#' @param groupMember Vector containing the current group member, i.e., which variables among the \code{contextVector} to be optimized in the current phase.
#' @param ... Further arguments passed to \code{fun}.
#' @param ubound upper bound for the group
#' @param lbound lower bound for the group
#' @param popSize population size. Ignored in cmaes_1
#' @return Returns a list containing x (decision variable) and y (objective value)
#' @export

cmaes_1 <- function(contextVector,func,groupMember,...,ubound=1,lbound=0,popSize=1){
  groupSize <- length(groupMember)
  individual <- contextVector[groupMember]

  bestPop <- cmaes::cma_es(individual,
                           fn=func,
                           contextVector=contextVector,groupMember=groupMember,...,
                           lower=lbound,
                           upper=ubound)
  bestObj <- fun(bestPop,...)

  return(list(x=bestPop,y=bestObj))
}


#' Wrapper function, the optimization itself is using the \code{DEoptim::DEoptim()} function.
#' Use this function as a template to create a population-based \code{phaseSolver}.
#'
#' @title Phase optimization using Differential Evolution
#' @param contextVector The current best solution
#' @param fun function to be optimized
#' @param groupMember Vector containing the current group member, i.e., which variables among the \code{contextVector} to be optimized in the current phase.
#' @param ... Further arguments passed to \code{fun}.
#' @param ubound upper bound for the group
#' @param lbound lower bound for the group
#' @param popSize population size.
#' @return Returns a list containing x (decision variable) and y (objective value)
#' @export

optimDE <- function(contextVector,fun,groupMember,...,ubound=1,lbound=0,popSize=length(groupMember)){
  groupSize <- length(groupMember)
  individual <- contextVector[groupMember]

  ##### generate population. Ignored here #####
  subpop <- t(randtoolbox::sobol(groupSize,popSize)) # vary the current group variables. Not used in this function. Serves as an example.
  subpop <- t(subpop)   # DEoptim is row major, need to transpose
  ##### end ignore #####

  optimized <- DEoptim::DEoptim(fn=fun,
                                lower=lbound,
                                upper=ubound,
                                control=DEoptim.control(NP=populationSize),
                                contextVector=contextVector,groupMember=groupMember,...) # initial pop can be added here if wanted using "initialpop=subpop"

  return(list(x=optimized$bestmem,y=optimized$bestval))
}


#' Wrapper function, the optimization itself is using the \code{cmaes::cma_es()} function.
#'
#' @title Phase optimization using (mu+mu)-CMA-ES
#' @param contextVector The current best solution
#' @param fun function to be optimized
#' @param groupMember Vector containing the current group member, i.e., which variables among the \code{contextVector} to be optimized in the current phase.
#' @param ... Further arguments passed to \code{fun}.
#' @param ubound upper bound for the group
#' @param lbound lower bound for the group
#' @param popSize population size.
#' @return Returns a list containing x (decision variable) and y (objective value)
#' @export

cmaes <- function(contextVector,func,groupMember,...,ubound=rep(1,length(groupMember)),lbound=rep(0,length(groupMember)),popSize=length(groupMember),maxit=10000){
  groupSize <- length(groupMember)
  individual <- contextVector[groupMember]

  best <- cmaes::cma_es(individual,
                           fn=func,
                           contextVector=contextVector,groupMember=groupMember,...,
                           lower=lbound,
                           upper=ubound,
                           control=list(mu=popSize,lambda=popSize,maxit=maxit))
  bestPop <- best$par
  bestObj <- best$value

  return(list(x=bestPop,y=bestObj,counts=best$counts))
}
