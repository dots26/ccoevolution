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
  if(is.matrix(population)){
    popSize <- nrow(population)
  }
  if(is.vector(population))
    popSize <- 1

  nVar <- length(contextVector)
  contextVector <- rep(contextVector,popSize)

  contextVector <- t(matrix(contextVector,nrow=nVar))
  contextVector[,groupMember] <- population


  objectiveValue <- mainfun(contextVector,...)
}

subfunctionCMA <- function(population,contextVector,groupMember,mainfun,...){
  if(is.matrix(population)){
    population <- t(population)
    popSize <- nrow(population)
  }
  if(is.vector(population))
    popSize <- 1

  nVar <- length(contextVector)
  contextVector <- rep(contextVector,popSize)

  contextVector <- t(matrix(contextVector,nrow=nVar))
  contextVector[,groupMember] <- population


  objectiveValue <- mainfun(contextVector,...)
}

tfun <- function(x,mainfun,...){
  res <- mainfun(x,...)
  return(t(res))
}

modelPrediction <- function(newPoint,model){
  # transpose because cmaes ordering is not compatible with spot
  if(is.matrix(newPoint))
    a <- predict(model,t(newPoint))$y
  return(a)
}
