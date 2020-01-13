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

  if(is.vector(population)) contextVector <- population

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

multilevel <- function(group,contextVector,mainfun,NP=NULL,
                       ubound,lbound,
                       budget=3000000,eval_interval=100000,
                       fname_prefix="datalog_",control=NULL,infill=c("ei","best","setmap"),
                       surrogate_type=1,bestval=Inf,...){
  nGroup <- length(group)
  infill <- infill[1]
  globalBest <- Inf
  if(is.null(NP)){
    NP <- rep(20,nGroup)
  }
  nVar <- length(contextVector)
  groupLength <- pracma::zeros(1,nGroup)
  for(i in 1:nGroup){
    groupLength[i] <- length(group[[i]])
  }

  pop <- InitializePopulationSobol(groupLength[1],NP[1],
                                   minVal = lbound[group[[1]]],
                                   maxVal = ubound[group[[1]]])

  currentGroup <- group[[1]]
  currentVector <- contextVector
  globalBestPop <- contextVector
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
                                   eval_interval,bestval=bestval,infill=infill,...)
    print(paste0('neval:',nEval))
    print(currentVector)
    valCollection <- append(valCollection,funValue$value)
    print(valCollection)
    if(funValue$value < globalBest){
      bestval <- funValue$value
      cv[currentGroup] <<- funValue$par
      globalBestPop <- cv
      globalBest <- funValue$value
      gb <<- globalBest
    }
  }
  # build surrogate here
  print('top level surrogate build')
  rf_model <- SPOT::buildKriging((pop),matrix(valCollection))
  # rf_model <- randomForest::randomForest((pop),(valCollection),ntree=1000)

  while(nEval < budget){
    # sequential optimization
    print("top level surrogate update")
    if(infill=="best"){
      fname<-bestPredModel
      # predictedBest <- sansde(fname=fname,
      #                         pop=pop,
      #                         model=rf_model,
      #                         bestmem = cv[currentGroup],
      #                         bestval = bestval,
      #                         Lbound=lbound[currentGroup],
      #                         Ubound=ubound[currentGroup],
      #                         control=list(itermax=100))$par
      predicted <- cma_es(fn=fname,
                              par=contextVector[currentGroup],
                              model=rf_model,
                              lower=lbound[currentGroup],
                              upper=ubound[currentGroup],
                              control=list(maxit=2000,mu=2,lambda=10,
                                           sigma=0.3*max(ubound[currentGroup]-lbound[currentGroup])))
      predictedBest <- predicted$par
    }else if(infill=="ei"){
      fname<-EIModel
      # predictedBest <- sansde(fname=fname,
      #                         pop=pop,
      #                         model=rf_model,
      #                         bestmem = cv[currentGroup],
      #                         bestval = bestval,
      #                         Lbound=lbound[currentGroup],
      #                         Ubound=ubound[currentGroup],
      #                         control=list(itermax=100),minVal=bestval)$par
      predicted <- cma_es(fn=fname,
                              par=contextVector[currentGroup],
                              model=rf_model,
                              lower=lbound[currentGroup],
                              upper=ubound[currentGroup],
                              minVal=gb,
                              control=list(maxit=10000,mu=10,lambda=20,
                                           sigma=0.3*max(ubound[currentGroup]-lbound[currentGroup])))
      predictedBest <- predicted$par
    }

    currentVector <- cv
    currentVector[currentGroup] <- predictedBest
    newVal <- EvaluateLowerLevel(fun= mainfun,
                                 currentVector = currentVector,
                                 nextGroup =  group[[2]],
                                 furtherGroup = group[c(-1,-2)],
                                 NP_next = NP[2],
                                 NP_further =  NP[c(-1,-2)],
                                 lbound=lbound,ubound=ubound,
                                 nEval,eval_interval,bestval=Inf,infill=infill,
                                 ...)
    print(paste('neval:',nEval,newVal$value))
    if(newVal$value < globalBest){
      cv[group[[2]]] <<- newVal$par
      cv[group[[1]]] <<- predictedBest
      globalBestPop <- cv
      globalBest <- funValue$value
      gb <<- globalBest
    }
    valCollection <- append(valCollection,newVal$value)
    pop <- rbind(pop,predictedBest)
    save(valCollection,pop,file='training.Rdata')
    print(bestval)
    print(t(pop[,1:4]))
    print(predicted$value)
    print(valCollection)

    # rf_model <- randomForest::randomForest(pop,valCollection,ntree=1000)#grf::regression_forest((pop),(valCollection),num.trees=1000)
    rf_model <- SPOT::buildKriging((pop),matrix(valCollection))
  }

  return(list(par=globalBestPop,value=globalBest))
}

bestPredModel.randomForest <- function(x,model){
  if(is.vector(x)){
    x <- matrix(x,nrow=1)
  }
  pred <- predict(model,x)
  return(pred)
}

bestPredModel.kriging <- function(x,model){
  if(is.vector(x)){
    x <- matrix(x,nrow=1)
  }
  pred <- predict(model,x)$y
  return(pred)
}

EIModel.kriging <- function(x,model,minVal){
  if(is.vector(x)){
    x <- matrix(x,nrow=1)
  }
  model$target <- c('y','ei')
  pred <- predict(model,x)$ei
  #EI <- SPOT:::expectedImprovement(pred$predictions,pred$variance,minVal)
  return(pred)
}

setmapModel <- function(x,model){
  pred <- predict(model,x,estimate.variance=T)
  EI <- SPOT:::expectedImprovement(pred$predictions,pred$variance,minVal)
  return(EI)
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
