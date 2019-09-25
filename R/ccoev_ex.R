#' Main function for the cooperative coevolution procedure.
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
#' cc_ex(nVar = 1000,fun=func,phaseSolver = cmaes,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
cc_ex <- function(contextVector=NULL,nVar,fun,phaseSolver=cmaes,budget=1000000,lbound=rep(-Inf,length(population)),ubound=rep(Inf,length(population)),nLowerLevel=round(0.1*nVar),...){
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

  # print(c(bestPop,bestObj))
  nCycle <- 9
  inactiveVar <- 1:nVar
  while( (length(inactiveVar>1) && (budget-nEval>0)) ){
    print('Throwing LASSO...') # regularization
    predictlasso<-glmnet::glmnet((population), objectiveValue)
    stopIndex <- max(which(predictlasso$df<nLowerLevel))
    stopLambda <- predictlasso$lambda[stopIndex]
    indices <- coef(predictlasso,s=stopLambda)
    indices <- as.matrix(indices)
    activeVariable <- indices[-1]
    activeVariable <- which(activeVariable!=0)

    nActive <- length(activeVariable)
    print(c('nActive=',nActive))
    nInactive <- nVar - nActive
    print(c('nInactive remaining=',nInactive))
    activeCopy <- activeVariable
    activeVariable <- inactiveVar[activeVariable]
    inactiveVar <- inactiveVar[-activeCopy]
    # end regularization


    print(c('Evaluating',activeVariable))
    print('Grouping...')

    grouping_control <- list(ubound= ubound[activeVariable]/2,lbound= lbound[activeVariable],delta=(ubound[activeVariable]-lbound[activeVariable])*0.3)
    grouping_result <- differential_grouping(nActive,fun=subfunction,grouping_control,contextVector=bestPop,groupMember=activeVariable,mainfun=fun,...)
    group <- grouping_result$group
    nEval <- nEval + grouping_result$nEval

    new_group <-NULL
    for(groupIndex in 1:length(group)) new_group <- append(new_group,list( activeVariable[group[[groupIndex]]]))

    group <- new_group
    # error checking on groups
    if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
    if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')

    for(cycleIndex in 1:nCycle){
      print('starting new cycle:')
      print(cycleIndex)
      # Optimize for each group
      groupBest <- NULL
      for(groupIndex in 1:length(group)) {
        nSubEval <- 0
        groupSize <- length(group[[groupIndex]])
        groupMember <- group[[groupIndex]]
        #print(c('Evaluating Group:',groupMember))

        # build surrogate model on regularized data
        # sample points for regularized data
        pop <- InitializePopulationLHS(groupSize,groupSize*1000,minVal=lbound[groupMember],maxVal=ubound[groupMember])
        obj <- subfunction(pop,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...)
        nSubEval <- nSubEval + groupSize*1000

        # optimize surrogate model
        start_time <- Sys.time()
        nCMAESiter <- 1000*groupSize
        for(cmaesIter in 1:nCMAESiter){
          randomForestModel <- SPOT::buildRandomForest(t(pop),matrix(obj,ncol=1))
          best<- cmaes::cma_es(par=contextVector[groupMember],fn=modelPrediction,model=randomForestModel,upper= ubound[groupMember],lower= lbound[groupMember],control = list())
          pop <- cbind(pop,best$par)
          obj <- append(obj,subfunction(best$par,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...))
        }
        rm(randomForestModel)
        timediff <-  Sys.time() - start_time
        nSubEval <- nSubEval + nCMAESiter
        groupBest <- cbind(groupBest,list(best=best,nEval=nSubEval))
      }

      # update the context vector for active variables
      print('updating context vector...')
      for(groupIndex in 1:length(group)){
        groupSize <- length(group[groupIndex])
        groupMember <- group[[groupIndex]]
        contextVector[groupMember] <- as.matrix(groupBest)[,groupIndex]$best$par
        nEval <- nEval + as.matrix(groupBest)[,groupIndex]$nEval
      }
      obj <- EvaluateIndividual(contextVector,fun=fun,...)
      nEval <- nEval + 1
      if(obj < bestObj){
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          bestPop <- contextVector
          bestObj <- obj
        }
      }
    }
    nVar <- nInactive

    print('drawing new population')
    population <- (randtoolbox::sobol(10000,nVar,scrambling = 3))*(ubound[inactiveVar]-lbound[inactiveVar])+lbound[inactiveVar]
    # Evaluate the whole population
    print('Evaluating population of inactive vars, subject to last active vars...')
    objectiveValue <- subfunction(population,contextVector=contextVector,groupMember=inactiveVar,mainfun=fun,...)
    nEval <- nEval + 10000
    bestPopIndex <- which.min(objectiveValue)
    # update the context vector for inactive variables
    contextVector[inactiveVar] <- population[bestPopIndex,]
    filename <- paste('remaining_inactive_is_',nInactive,'.Rdata',sep ='')
  }
  return(list(x=bestPop,y=bestObj))
}

