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
InitializePopulationLHS <- function(ncol,nrow,minVal=0,maxVal=1) {
  #population<-optimumLHS(n=numberOfIndividuals,k=chromosomeLength,maxSweeps=10,eps=.1,verbose=FALSE)
  population<-lhs::randomLHS(n=ncol,k=nrow)

  population<-t(population)
  population<-population * (maxVal-minVal) + minVal

  return(population)
}

InitializePopulationSobol <- function(ncol,nrow,minVal=0,maxVal=1) {
  #population<-optimumLHS(n=numberOfIndividuals,k=chromosomeLength,maxSweeps=10,eps=.1,verbose=FALSE)
  population<- randtoolbox::sobol(n=ncol,dim=nrow,scrambling = 3)

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

multilevel <- function(group,mainfun,NP=NULL,
                       ubound,lbound,
                       budget=3000000,eval_interval=100000,
                       fname_prefix="datalog_",control=NULL,...){
  nGroup <- length(group)
  globalBest <- Inf
  if(is.null(NP)){
    NP <- rep(20,nGroup)
  }
  nVar <- length(cv)
  groupLength <- pracma::zeros(1,nGroup)
  for(i in 1:nGroup){
    groupLength[i] <- length(group[[i]])
  }

  #assembledPopulation <- zeros(popsize,nVar)
  #pop <- vector(mode='list',length = nGroup)
  pop <- InitializePopulationSobol(groupLength[1],NP[1],
                                 minVal = lbound[group[[1]]],
                                 maxVal = ubound[group[[1]]])

  currentGroup <- group[[1]]
  currentVector <- cv
  globalBestPop <- cv
  #TO-DO: globalbest? bestval? check!
  print(currentVector)
  print('S')
  bestval <- Inf
  valCollection <- NULL
  for(i in 1:NP[1]){
    currentVector[currentGroup] <- pop[i,]
    funValue <- EvaluateLowerLevel(mainfun,
                                   currentVector,
                                   group[[2]],
                                   group[c(-1,-2)],
                                   NP[2],
                                   NP[c(-1,-2)],
                                   lbound=lbound,ubound=ubound,
                                   prevnEval = nEval,
                                   eval_interval,bestval=bestval,...)
    print(paste0('neval:',nEval))
    print(currentVector)
    valCollection <- append(valCollection,funValue$value)
    print(valCollection)
    if(funValue$value < globalBest){
      bestval <- funValue$value
      cv[currentGroup] <<- funValue$par
      globalBestPop <- cv
      globalBest <- funValue$value
    }
  }
  # build surrogate here
  print('top level surrogate build')
  print(pop)
  print(valCollection)
  # rf_model <- SPOT::buildKriging((pop),matrix(valCollection))
  rf_model <- randomForest::randomForest((pop),(valCollection),ntree=1000)

  pop <- InitializePopulationSobol(ncol=groupLength[1],nrow=NP[1],
                                 minVal = lbound[group[[1]]],
                                 maxVal = ubound[group[[1]]])

  while(nEval < budget){
    # sequential optimization
    print("top level surrogate update")
    predictedBest <- sansde(fname=evalModel,
                            pop=pop,
                            model=rf_model,
                            bestmem = cv[currentGroup],
                            bestval = bestval,
                            Lbound=lbound[currentGroup],
                            Ubound=ubound[currentGroup],
                            control=list(itermax=100))$par

    currentVector <- cv
    currentVector[currentGroup] <- predictedBest
    newVal <- EvaluateLowerLevel(fun= mainfun,
                                 currentVector = currentVector,
                                 nextGroup =  group[[2]],
                                 furtherGroup = group[c(-1,-2)],
                                 NP_next = NP[2],
                                 NP_further =  NP[c(-1,-2)],
                                 lbound=lbound,ubound=ubound,
                                 nEval,eval_interval,bestval=Inf,
                                 ...)
    print(paste0('neval:',nEval,newVal$value))
    if(newVal$value < globalBest){
      cv[currentGroup] <<- newVal$par
      globalBestPop <- cv
      globalBest <- funValue$value
    }
    valCollection <- append(valCollection,newVal$value)
    pop <- rbind(pop,predictedBest)

    print(bestval)
    print(t(pop[,1:4]))
    print(valCollection)

    rf_model <- randomForest::randomForest((pop),(valCollection),ntree=1000)
    # rf_model <- SPOT::buildKriging((pop),matrix(valCollection))
  }

  return(list(par=globalBestPop,value=globalBest))
}



evalModel <- function(x,model){
  pred <- predict(model,x)
  return(pred)
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
