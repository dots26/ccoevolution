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
#' population <- (randtoolbox::sobol(5,100))
#' optimum <- rep(0,100)
#' func <- f1cec
#' ctrl <- list(lbound=rep(-100,100),ubound=rep(100,100),delta=rep(100,100))
#' group <- differential_grouping(ncol(population),func,ctrl,o=optimum)
#' cc(nVar = 100,fun=func,phaseSolver = cmaes,group=group,nCycle = 2,o=optimum,lbound=rep(-100,100),ubound=rep(100,100))
#' @export
cc <- function(contextVector=NULL,nVar,fun,phaseSolver=cmaes,budget=1000000,group=NULL,grouping_control=list(),nCycle=9,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),...){
  doParallel::registerDoParallel()
  print(c('Ncores=',foreach::getDoParWorkers() ))
  nEval <- 0
  if(is.null(contextVector)){
    nEval <- nEval + 10000      #10000
    population <- (randtoolbox::sobol(nEval,nVar,scrambling = 3))*(ubound-lbound)+lbound
    # Evaluate the whole population
    print('Evaluating initial population...')
    objectiveValue <- fun(population,...)
    bestPopIndex <- which.min(objectiveValue)
    bestPop <- population[bestPopIndex,]
    bestObj <- min(objectiveValue)

    contextVector <- bestPop
  }else{
    bestPop <- contextVector
    bestObj <- fun(bestPop,...)
  }

  # grouping
  if(is.null(group)){ # if no group is supplied, then the group is determined by differential grouping
    print('Grouping...')
    dg <- differential_grouping(ncol(population),fun,grouping_control,...)
    nEval <- nEval + dg$nEval #510501
    group <- dg$group
    print(c('nEval',nEval))
    save(group, file='group.Rdata')
  }
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')

  nCycle <- 9

  newContextVector <- contextVector
  for(cycleIndex in 1:nCycle){
    print(c('Cycle Number:',cycleIndex))
    # Optimize for each group
    for(groupIndex in 1:length(group)){
      groupSize <- length(group[groupIndex])
      groupMember <- group[[groupIndex]]

      # generate population for current phase
      best<-cmaes::cma_es(par = contextVector[groupMember],
                          fn=subfunction,
                          contextVector=contextVector,
                          groupMember = groupMember,mainfun=fun,...,
                          upper= ubound[groupMember],
                          lower= lbound[groupMember],
                          control = list(mu = 50,lambda=50,maxit=1080))

      nEval <- nEval + best$counts[1]
      print(c('best count',best$counts[1]))
      newContextVector[groupMember] <- best$par
    }
    print('Updating context vector')
    contextVector <- newContextVector
    obj <- fun(newContextVector,...)
    nEval <- nEval + 1

    if(obj < bestObj){
      if((budget-nEval)>0){ # only update if it doesnt exceed budget
        bestPop <- newContextVector
        bestObj <- obj
      }
    }

    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget,budget,nEval))
#    save(list=ls(),file=paste('datacc_',seed,'.Rdata',sep=''))
  }
  return(list(x=bestPop,y=bestObj))
}



modelPrediction <- function(newPoint,model){
  predict(model,matrix(newPoint,nrow=1))$y
}
