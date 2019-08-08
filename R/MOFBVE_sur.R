#' Main function for the cooperative coevolution MOFBVE with multilevel framework
#'
#' @title Cooperative coevolution.
#' @param population Initial population
#' @param fun The objective function object to be solved
#' @param phaseSolver Optimizer/solver for phase-optimization step.
#' @param group Vector of list. Each list contains a group of non-separable variables. Available method: \code{differential_grouping}.
#' @param grouping_control Grouping parameters. The members depends on which grouping being used.
#' @param nCycle Number of cycles to be run.
#' @param ... Further arguments passed to \code{fun}.
#' @examples
#' optimum <- rep(13,1000)
#' func <- f1cec
#' ctrl <- list(lbound=rep(-100,1000),ubound=rep(100,1000),delta=rep(20,1000))
#' MOFBVE_surrogate(nVar = 1000,fun=func,phaseSolver = cmaes,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
MOFBVE_surrogate <- function(contextVector=NULL,nVar,group=NULL,fun,phaseSolver=cmaes,budget=1000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),nLevel=4,...){
  doParallel::registerDoParallel()
  print(c('Ncores=',foreach::getDoParWorkers() ))
  nEval <- 0
  if(is.null(contextVector)){
    nEval <- nEval + 10000
    population <- (randtoolbox::sobol(nEval,nVar,scrambling = 3))*(ubound-lbound)+lbound
    # Evaluate the whole population
    print('Evaluating initial population...')
    objectiveValue <- fun(population,...)
    bestPopIndex <- which.min(objectiveValue)
    bestPop <- population[bestPopIndex,]
    bestObj <- min(objectiveValue)

    contextVector <- bestPop
  }

  #if(is.null(group)){
    r<- 2
    a <- sensitivity::morris(model=fun,
                             factor=nVar,
                             r = r,
                             design = list(type='oat',levels=8,grid.jump=1),
                             binf=lbound,
                             bsup=ubound,
                             scale=F,...)

    nEval <- nEval + r*(nVar+1)
    mu.star <- apply(a$ee, 2, function(a) mean(abs(a)))
    sigma <- apply(a$ee, 2, sd)

    print('Using k-means...') # clustering
    assignedGroup <- kmeans(cbind(mu.star,sigma),nLevel)$cluster

    EF <- NULL
    group <- NULL
    sigmaMuAll <- sum(mu.star)
    for(groupIndex in 1:nLevel){
      sigmaMu <- sum(mu.star[which(assignedGroup==groupIndex)])
      EF <- append(EF,sigmaMu/sigmaMuAll)
      group <- append(group,list(which(assignedGroup==groupIndex)))
    }
    clusterOrder <- order(EF,decreasing = T)
    save(group,file='MOFBVE_sur_group.Rdata')
#  }else{

#  }
  print('Groups assigned')
  print(group[clusterOrder])
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval

  while((budget-nEval)>0 ){
    for(groupIndex in clusterOrder[1:length(group)]) {
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]

      pop <- InitializePopulationLHS(groupSize,groupSize*10,minVal=lbound[groupMember],maxVal=ubound[groupMember])
      obj <- subfunction(pop,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...)
      nEval <- nEval + groupSize*10

      # optimize surrogate model
      nCMAESiter <- 10*groupSize
      for(cmaesIter in 1:nCMAESiter){
        randomForestModel <- SPOT::buildRandomForest(pop,obj)
        best<- cmaes::cma_es(par=contextVector[groupMember],fn=modelPrediction,model=randomForestModel,upper= ubound[groupMember],lower= lbound[groupMember],control = list())
        pop <- rbind(pop,best$par)
        obj <- append(obj,subfunction(best$par,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...))
      }
      rm(randomForestModel)
      nEval <- nEval + nCMAESiter

      print('updating context vector for interconnection step...')
      contextVector[groupMember] <- best$par
      obj <- fun(contextVector,...)
      nEval <- nEval + 1
      if(obj < bestObj){
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          bestPop <- contextVector
          bestObj <- obj
        }
      }

      # Interconnection
      best <- cmaes::cma_es(par = contextVector,fn = fun,...,lower = lbound,upper=ubound,control = list(maxit=groupSize*5,mu=groupSize,lambda=groupSize))
      nEval <- nEval + best$counts[1]
      print('Interconnection step finished, updating context vector...')
      contextVector <- best$par
      obj <- fun(contextVector,...)
      nEval <- nEval + 1
      if(obj < bestObj){
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          bestPop <- contextVector
          bestObj <- obj
        }
      }
    }
    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget))
  }
  return(list(x=bestPop,y=bestObj))
}
