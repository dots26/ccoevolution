#' Main function for the cooperative coevolution MOFBVE with multilevel framework
#'
#' @title Cooperative coevolution.
#' @param population Initial population
#' @param fun The objective function object to be solved
#' @param group Vector of list. Each list contains a group of non-separable variables. Available method: \code{differential_grouping}.
#' @param grouping_control Grouping parameters. The members depends on which grouping being used.
#' @param nCycle Number of cycles to be run.
#' @param ... Further arguments passed to \code{fun}.
#' @examples
#' optimum <- rep(13,1000)
#' func <- f1cec
#' ctrl <- list(lbound=rep(-100,1000),ubound=rep(100,1000),delta=rep(20,1000))
#' SACC(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
SACC <- function(contextVector=NULL,nVar,fun,group=NULL,budget=1000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),nLevel=4,evalInterval=100000,...){
  doParallel::registerDoParallel()
  #groupSize <- 20
  # print(c('Ncores=',foreach::getDoParWorkers() ))
  nEval <- 0
  convergence_history <- NULL
  # print('sens')
  nEval <- nEval + 10000
  population <- (randtoolbox::sobol(nEval,nVar,scrambling = 3))*(ubound-lbound)+lbound
  # Evaluate the whole population
  # print('Evaluating initial population...')
  objectiveValue <- fun(population,...)
  bestPopIndex <- which.min(objectiveValue)
  bestPop <- population[bestPopIndex,]
  bestObj <- min(objectiveValue)

  contextVector <- bestPop
  # r<- 20
  # a <- sensitivity::morris(model=fun,
  #                          factor=nVar,
  #                          r = r,
  #                          design = list(type='oat',levels=8,grid.jump=1),
  #                          binf=lbound,
  #                          bsup=ubound,
  #                          scale=F,...)
  # bestPopIndex <- which.min(a$y)
  # bestPop <- a$X[bestPopIndex,]
  # bestObj <- min(a$y)
  prevLevel <- NULL
  #group <- NULL
  if(is.null(group)){
    for(groupingIndex in 1:(nLevel-1)){
      nLevelVar <- floor(nVar/nLevel)
      print(paste0('Throwing LASSO...#',groupingIndex)) # regularization
      predictlasso<-glmnet::glmnet((population), objectiveValue)
      stopIndex <- max(which(predictlasso$df<nLevelVar*groupingIndex))
      stopLambda <- predictlasso$lambda[stopIndex]
      indices <- coef(predictlasso,s=stopLambda)
      indices <- as.matrix(indices)
      activeVariable <- indices[-1]
      activeVariable <- which(activeVariable!=0)
      activeVariable <- (setdiff(activeVariable,prevLevel))
      prevLevel <- unique(c(prevLevel,activeVariable))
      group <- append(group,list(activeVariable))
    }
    lastGroup <- setdiff(1:nVar,prevLevel)
    group <- append(group,list(lastGroup))
  }
  print(group)
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval
  convergence_history <- append(convergence_history,bestObj)
  #saved
  while((budget-nEval)>0 ){
    for(groupIndex in 1:length(group)) {
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]

      print(c('optimizing group',groupIndex,'with',groupSize,'members'))
      # group optimization
      best<- cma_es(contextVector[groupMember],
                    fn = subfunctionCMA,
                    contextVector = contextVector,
                    groupMember = groupMember,mainfun=fun,...,
                    lower = lbound[groupMember],
                    upper=ubound[groupMember],
                    control = list(vectorized=T,
                                   mu=groupSize,lambda=groupSize,
                                   maxit=900,
                                   sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                   diag.value=T))
      nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
      if(nlogging_this_layer>0){
        for(i in 1:nlogging_this_layer){
          nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
          nGeneration_to_consider <- floor(nEval_to_logging/groupSize)
          bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
          convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
          # print(convergence_history)
        }
      }
      nEval <- nEval + best$counts[1]
      print('updating context vector for interconnection step...')
      if((budget-nEval)>0){ # only update if it doesnt exceed budget
        if(!is.null(best$par)){
          contextVector[groupMember] <- best$par
          obj <- best$value
          if(obj < bestObj){
            bestPop <- contextVector
            bestObj <- obj
            # print('Update:')
            # print(bestObj)
          }
        }
      }else{
        break
      }
      # Interconnection
      best <- cma_es(contextVector,
                     fn = subfunctionCMA,
                     contextVector=contextVector,
                     groupMember=1:nVar,
                     mainfun=fun,
                     ...,
                     lower = lbound,
                     upper=ubound,
                     control = list(vectorized=T,mu=groupSize,lambda=groupSize,
                                    maxit=90,
                                    sigma=0.3*max(ubound-lbound),
                                    diag.value=T))
      nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
      if(nlogging_this_layer>0){
        for(i in 1:nlogging_this_layer){
          nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
          nGeneration_to_consider <- floor(nEval_to_logging/groupSize)
          bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
          convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
          # print(convergence_history)
        }
      }
      nEval <- nEval + best$counts[1]

      print('Interconnection step finished, updating context vector...')
      if((budget-nEval)>0){ # only update if it doesnt exceed budget
        if(!is.null(best$par)){
          contextVector <- best$par
          obj <- best$value
          if(obj < bestObj){
            bestPop <- contextVector
            bestObj <- obj
            # print('Update:')
            # print(bestObj)
          }
        }
      }else{
        break
      }
      leftBudget <- budget - nEval
      print(c('Comp budget left:',leftBudget,budget,nEval))
    #  # print(paste('dataMOFBVE_',seed,'.Rdata',sep=''))
     # save(list=ls(),file=paste('dataMOFBVE_',seed,'.Rdata',sep=''))
    }
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
