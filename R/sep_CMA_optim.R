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
#' optim_sepCMA(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
optim_sepCMA <- function(contextVector=NULL,nVar,fun,...,
                         group=NULL,
                         budget=3000000,
                         lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),evalInterval=100000){

  contextVector <- runif(nVar)
  convergence_history <- NULL
  nEval <- 0
  groupSize <- nVar
  groupMember <- 1:nVar
  bestObj <- Inf

  while((budget-nEval)>10000){ # restart if at least enough for 10 iterations
    groupSize <- nVar
    groupMember <- 1:nVar
    CMAES_control <- list(vectorized=T,
                          mu=groupSize,lambda=groupSize,
                          maxit=(budget-nEval)/1000-1,
                          sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                          diag.value=T)
    # group optimization
    best<- sep_cma_es(contextVector[groupMember],
                      fn = subfunctionCMA,
                      contextVector = contextVector,
                      groupMember = groupMember,mainfun=fun,...,
                      lower = lbound[groupMember],
                      upper = ubound[groupMember],
                      control = CMAES_control)

    nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
    if(nlogging_this_layer>0){
      for(i in 1:nlogging_this_layer){
        nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
        nGeneration_to_consider <- floor(nEval_to_logging/groupSize)
        bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
        convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
      }
    }
    nEval <- nEval + best$counts[1]
    print('updating context vector...')
    if((budget-nEval)>0){ # only update if it doesnt exceed budget
      if(!is.null(best$par)){
        contextVector[groupMember] <- best$par
        obj <- best$value
        if(obj < bestObj){
          bestPop <- contextVector
          bestObj <- obj
        }
      }
    }else{
      break
    }

    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget,budget,nEval))
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
