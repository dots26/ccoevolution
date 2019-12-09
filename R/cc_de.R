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
#' cc_de(nVar = 100,fun=func,group=group,nCycle = 2,o=optimum,lbound=rep(-100,100),ubound=rep(100,100))
#' @export
cc_de <- function(contextVector=NULL,nVar,fun,budget=1000000,group=NULL,grouping_control=list(),nCycle=9,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),evalInterval=100000,...){
  doParallel::registerDoParallel()
  print('new code')
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
  print(group)
  print(dg$separable)
  nGroup <- length(group)
  # error checking on groups
  # if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  # if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  DE_control <- vector(mode = "list",length=nGroup+1)

  if(nGroup>0){
    for(groupIndex in 1:nGroup){
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]
      DE_control[[groupIndex]] <- list(ccm=0.5,itermax=50)
    }
  }
  DE_control[[nGroup+1]]  <- list(ccm=0.5,itermax=50)


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
        NP <- 50
        # generate population for current phase
        pop <- t(InitializePopulationLHS(NP,groupSize,lbound[groupMember],ubound[groupMember]))
        best <-  sansde(pop=pop,
                        bestmem=bestPop[groupMember],
                        bestval=bestObj,
                        fname = subfunction,
                        contextVector = contextVector,
                        groupMember = groupMember,mainfun=fun,...,
                        Lbound = lbound[groupMember],
                        Ubound =ubound[groupMember],
                        control = DE_control[[groupIndex]])
        DE_control[[groupIndex]]$ccm <- best$ccm
        nlogging_this_layer <- floor((nEval+best$used_FEs)/evalInterval)-floor(nEval/evalInterval)
        if(nlogging_this_layer>0){
          for(i in 1:nlogging_this_layer){
            nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
            nGeneration_to_consider <- floor(nEval_to_logging/NP)
            bestObj_logging <- min(best$tracerst)
            convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
            # print(c('conv',convergence_history))
          }
        }
        nEval <- nEval + best$used_FEs
        ## print(c('best count',best$counts[1]))



        leftBudget <- budget - nEval
        print(c('Comp budget left:',leftBudget,budget,nEval))

        if(best$bestval < bestObj){
          if((budget-nEval)>0){
            newContextVector[groupMember] <- best$bestmem
            bestPop <- newContextVector
            bestObj <- best$bestval
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
      NP <- 50
      pop <- t(InitializePopulationLHS(NP,groupSize,lbound[groupMember],ubound[groupMember]))
      best<- sansde(pop=pop,
                    bestmem=bestPop[groupMember],
                    bestval=bestObj,
                    fname = subfunction,
                    contextVector = contextVector,
                    groupMember = groupMember,mainfun=fun,...,
                    Lbound = lbound[groupMember],
                    Ubound =ubound[groupMember],
                    control = DE_control[[nGroup+1]])
      nlogging_this_layer <- floor((nEval+best$used_FEs)/evalInterval)-floor(nEval/evalInterval)
      if(nlogging_this_layer>0){
        for(i in 1:nlogging_this_layer){
          nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
          nGeneration_to_consider <- floor(nEval_to_logging/groupSize)
          bestObj_logging <- min(best$tracerst)
          convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
          # print(c('conv',convergence_history))
        }
      }
      nEval <- nEval + best$used_FEs
      ## print(c('best count',best$counts[1]))

      leftBudget <- budget - nEval
      print(c('Comp budget left:',leftBudget,budget,nEval))

      if(best$bestval < bestObj){
        if((budget-nEval)>0){
          newContextVector[groupMember] <- best$bestmem
          bestPop <- newContextVector
          bestObj <- best$bestval
          # print('Update:')
          # print(bestObj)
        }else{
          break
        }
      }
    }
    contextVector <- newContextVector
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}


