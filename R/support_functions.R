#' Evaluate individual with the specified test function. Non-feasible solution are given Inf as objective values.
#' @title Evaluate objective values of a single individual
#' @param individual The individual to be evaluated
#' @param fun A string containing which problem are being solved.
#' @param ... Further parameters used by \code{fun}
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(8,1)
#' EvaluateIndividual(individual,WFG4,3) # the 3 is passed to WFG4 nObj
#' @export
EvaluateIndividual <- function(individual,fun,...){
  nVar <- length(individual)

  objective <- do.call((fun),args = list(individual,...))
  return(objective)
}



#' Evaluate a population with the specified test function. Non-feasible solution are given Inf as objective values.
#' @title Evaluate objective value of a set of individuals
#' @param pop The population to be evaluated
#' @param fun A string containing which problem are being solved.
#' @param ... Further parameters used by \code{fun}
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' pop <- runif(8,50) # 8 variables, 50 individuals
#' EvaluateIndividual(pop,WFG4,3) # the 3 is passed to WFG4 nObj
#' @export
EvaluatePopulation <- function(pop,fun,...){
  popSize <- ncol(pop)
  popObjective <- NULL
#  timest <- Sys.time()
  popObjective <- foreach::foreach(individualIndex=1:popSize,.combine = 'cbind')%dopar%{
    indObjective <- EvaluateIndividual(pop[,individualIndex],fun,...)
  }
#  timediff<-Sys.time()-timest

  #return(list(pop=popObjective,time=timediff))
  return(popObjective)
}

#' Create initial sample using LHS method. The variables will be ranged between 0-1
#' @title Initialize population with Latin Hypercube Sampling
#' @param numberOfIndividuals The number of individual in the population. Integer > 0.
#' @param chromosomeLength The number of variables per individual
#' @param binaryEncoding whether to use binary encoding or real encoding. Default to FALSE
#' @param minVal Minimum value of the resulting sample
#' @param maxVal Maximum value of the resulting sample
#' @param samplingMethod Not used
#' @return A matrix of size chromosomeLength x nIndividual.
#' @examples
#' nVar <- 14
#' nIndividual <- 100
#' InitializePopulationLHS(nIndividual,nVar,FALSE)
#'
InitializePopulationLHS <- function(numberOfIndividuals,chromosomeLength,minVal=0,maxVal=1) {
  #population<-optimumLHS(n=numberOfIndividuals,k=chromosomeLength,maxSweeps=10,eps=.1,verbose=FALSE)
  population<-lhs::randomLHS(n=numberOfIndividuals,k=chromosomeLength)

  population<-t(population)
  population<-population * (maxVal-minVal) + minVal

  return(population)
}

subfunction <- function(population,contextVector,groupMember,mainfun,...){
  if(is.matrix(population))
    popSize <- nrow(population)
  if(is.vector(population))
    popSize <- 1

  nVar <- length(contextVector)
  contextVector <- rep(contextVector,popSize)

  contextVector <- t(matrix(contextVector,nrow=nVar))
  contextVector[,groupMember] <- population

   #print(mainfun)
  objectiveValue <- mainfun(contextVector,...)
}
