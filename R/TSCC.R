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
#' TSCC(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
TSCC <- function(contextVector=NULL,nVar,fun,group=NULL,budget=1000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),nLevel=4,evalInterval=100000,...){
  doParallel::registerDoParallel()
  # print(c('Ncores=',foreach::getDoParWorkers() ))
  #groupSize <- 20
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

  a <- randomForest::randomForest(x=population,y=objectiveValue,importance=T)

  impo_acc <- a$importance[,1]
  impo_ranked <- order(impo_acc,decreasing = T)

  if(is.null(group)){
    interval <- floor(nVar/nLevel)
    for(i in 0:(nLevel-2)){
      start <- i*interval+1
      end <- (i+1)*interval
      group <- append(group,list(impo_ranked[start:end]))
    }
    group <- append(group,list(impo_ranked[((nLevel-1)*interval+1):nVar]))
  }

  # print('Secondary Grouping...')
  new_group <- NULL
  dg <- NULL
  for(i in 1:nLevel){
    gctrl <- list(lbound=lbound[group[[i]]],ubound=ubound[group[[i]]])
    subgroup <- DG2(length(group[[i]]),
                    subfunction,
                    control=gctrl,
                    contextVector=contextVector,
                    groupMember=group[[i]],
                    mainfun=fun,...)
    dg <- cbind(dg,subgroup)
    print(subgroup)
    nEval <- nEval + dg[,i]$nEval
  }

  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval
  convergence_history <- append(convergence_history,bestObj)
  #saved
  while((budget-nEval)>0 ){
    for(clusterIndex in 1:nLevel){
      # print(paste('Optimizing Variable level',clusterIndex))
      this.cluster <- group[[clusterIndex]]
      currentClusterGrouping <- dg[,clusterIndex]$group
      sep <- group[[clusterIndex]][dg[,clusterIndex]$separable]

      # optimize separable
      groupMember <- sep
      groupSize <- length(groupMember)

      if(groupSize>0){
       # # print(c('optimizing separable variables'))
        # group optimization
        best<- sep_cma_es(contextVector[groupMember],
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
        print('Interconnection step...')
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
      }

      # optimize non-separable
      print('Optimize current group non-seps...')
      if(length(currentClusterGrouping)>0){
        for(groupIndex in 1:length(currentClusterGrouping)) {
          groupMember <- this.cluster[currentClusterGrouping[[groupIndex]]]
          groupSize <- length(groupMember)

          best<- cma_es(contextVector[groupMember],
                        fn = subfunctionCMA,
                        contextVector = contextVector,
                        groupMember = groupMember,
                        mainfun=fun,...,
                        lower = lbound[groupMember],
                        upper=ubound[groupMember],
                        # control = list(mu=groupSize,lambda=groupSize,maxit=2000))
                        control = list(mu=groupSize,lambda=groupSize,
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

          print('updating context vector (non-separable) for interconnection step...')
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
            }else{
        #      # print('is null')
            }
          }else{
            break
          }
        }
        # Interconnection
        best <- cma_es(par = contextVector,
                       fn = subfunctionCMA,
                       contextVector = contextVector,
                       groupMember = 1:nVar,
                       mainfun=fun,...,
                       lower = lbound,
                       upper=ubound,
                       # control = list(vectorized=T,maxit=1000,mu=20,lambda=20))
                       control = list(vectorized=T,
                                      maxit=90,
                                      mu=groupSize,lambda=groupSize,
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


        # print('Interconnection step (nonseparable) finished, updating context vector...')
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
          }else{
        #    # print('is null')
          }
        }else{
          break
        }
      }
      leftBudget <- budget - nEval
      print(c('Comp budget left:',leftBudget,budget,nEval))
    #  save(list=ls(),file=paste('dataTSCC_',seed,'.Rdata',sep=''))
    }
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
