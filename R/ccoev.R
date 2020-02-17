#' Main function for the cooperative coevolution procedure.
#'
#' @title Cooperative coevolution.
#' @param population Initial population
#' @param fun The objective function object to be solved
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
#' cc(nVar = 100,fun=func,group=group,nCycle = 2,o=optimum,lbound=rep(-100,100),ubound=rep(100,100))
#' @export
cc <- function(contextVector=NULL,nVar,
               fun,
               budget=1000000,
               group=NULL,grouping_control=list(),
               nCycle=9,
               lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),
               evalInterval=100000,...){
  #doParallel::registerDoParallel()
  convergence_history <- NULL
  nEval <- 0
  #groupSize <- 20
  if(is.null(contextVector)){
    nEval <- nEval + 10000     #10000
    population <- (randtoolbox::sobol(nEval,nVar,scrambling = 3))*(ubound-lbound)+lbound
    # Evaluate the whole population
    ## print('Evaluating initial population...')
    objectiveValue <- fun(population,...)
    bestPopIndex <- which.min(objectiveValue)
    bestPop <- population[bestPopIndex,]
    bestObj <- min(objectiveValue)

    contextVector <- bestPop
  }else{
    bestPop <- contextVector
    bestObj <- fun(bestPop,...)
  }
  convergence_history <- append(convergence_history,bestObj)

  # grouping
  if(is.null(group)){ # if no group is supplied, then the group is determined by differential grouping
    ## print('Grouping...')
    dg <- DG2(length(bestPop),fun,grouping_control,...)
    nEval <- nEval + dg$nEval
    group <- dg$group

  }
  nlogging_this_layer <- floor(nEval/evalInterval)
  if(nlogging_this_layer>0){
    for(i in 1:nlogging_this_layer){
      nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
      convergence_history <- append(convergence_history,bestObj)
      # print(c('conv',convergence_history))
    }
  }
  print(nlogging_this_layer)
  print(nEval)
  print(group)
  print(dg$separable)

  # error checking on groups
  # if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  # if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')

  newContextVector <- contextVector
  leftBudget <- budget - nEval
  while(leftBudget>0){
    leftBudget <- budget - nEval
    ## print(c('Comp budget left:',leftBudget,budget,nEval))
    # Optimize for each group (non-sep)
    for(groupIndex in 1:(length(group))) {
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]
      if(groupSize>0){
        # generate population for current phase
        best<-cma_es(par = contextVector[groupMember],
                     fn=subfunctionCMA,
                     contextVector=contextVector,
                     groupMember = groupMember,mainfun=fun,...,
                     upper= ubound[groupMember],
                     lower= lbound[groupMember],
                     control = list(vectorized=T,
                                    mu=(groupSize),lambda=groupSize,
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
            # print(c('conv',convergence_history))
          }
        }
        nEval <- nEval + best$counts[1]
        ## print(c('best count',best$counts[1]))
        if(!is.null(best$par))
          newContextVector[groupMember] <- best$par

        leftBudget <- budget - nEval
        print(c('Comp budget left:',leftBudget,budget,nEval))
        print(best$value)
        if(best$value < bestObj){
          if((budget-nEval)>0){
            bestPop <- newContextVector
            bestObj <- best$value
            # print('Update:')
            # print(bestObj)
          }else{
            break
          }
        }
        # save(list=ls(),file=paste('datacc_',seed,'.Rdata',sep=''))
      }
    }

    # separable part
    groupSize <- length(dg$separable)
    groupMember <- dg$separable
    if(groupSize>0){
      best<-sep_cma_es(par = contextVector[groupMember],
                       fn=subfunctionCMA,
                       contextVector=contextVector,
                       groupMember = groupMember,mainfun=fun,...,
                       upper= ubound[groupMember],
                       lower= lbound[groupMember],
                       control = list(vectorized=T,
                                      mu=(groupSize),lambda=(groupSize),
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
          # print(c('conv',convergence_history))
        }
      }
      nEval <- nEval + best$counts[1]
      ## print(c('best count',best$counts[1]))
      if(!is.null(best$par)){
        newContextVector[groupMember] <- best$par
      }

      leftBudget <- budget - nEval
      print(c('Comp budget left:',leftBudget,budget,nEval))
      print(c( bestObj,best$value))

      if(best$value < bestObj){
        if((budget-nEval)>0){
          bestPop <- newContextVector
          bestObj <- best$value
          # print('Update:')

        }else{
          break
        }
      }
    }
    contextVector <- newContextVector
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}


