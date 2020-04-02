#' Main function for the cooperative coevolution MOFBVE with multilevel framework
#'
#' @title Cooperative coevolution.
#' @param population Initial population
#' @param fun The objective function object to be solved
#' @param phaseSolver Optimizer/solver for phase-optimization step.
#' @param group Vector of list. Each list contains a group of non-separable variables. Available method: \code{DG2}.
#' @param grouping_control Grouping parameters. The members depends on which grouping being used.
#' @param nCycle Number of cycles to be run.
#' @param ... Further arguments passed to \code{fun}.
#' @examples
#' optimum <- rep(13,1000)
#' func <- f1cec
#' ctrl <- list(lbound=rep(-100,1000),ubound=rep(100,1000),delta=rep(20,1000))
#' mofbve_grouped(nVar = 1000,fun=func,phaseSolver = cmaes,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
mofbve_grouped <- function(contextVector=NULL,nVar,group=NULL,fun,phaseSolver=cmaes,budget=1000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),nLevel=4,...){
  doParallel::registerDoParallel()
  print(c('Ncores=',foreach::getDoParWorkers() ))
  nEval <- 0
  print(nVar)
  r<- 20
  a <- sensitivity::morris(model=fun,
                           factor=nVar,
                           r = r,
                           design = list(type='oat',levels=8,grid.jump=1),
                           binf=lbound,
                           bsup=ubound,
                           scale=F,...)


  bestPopIndex <- which.min(a$y)
  bestPop <- a$X[bestPopIndex,]
  bestObj <- min(a$y)

  contextVector <- bestPop

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
  group <- group[clusterOrder]

  print('Grouping...')
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
    nEval <- nEval + dg[,i]$nEval
  }

  print(dg)

  print('Groups assigned')
  print(c('Budget left:',budget-nEval))
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval

  while((budget-nEval)>0 ){
    for(clusterIndex in clusterOrder){
      print(paste('Optimizing cluster no.',clusterIndex))
      this.cluster <- group[[clusterIndex]]
      currentClusterGrouping <- dg[,clusterIndex]$group
      sep <- group[[clusterIndex]][dg[,clusterIndex]$separable]

      # optimize separable
      groupMember <- sep
      groupSize <- length(groupMember)
      if(groupSize>0){
        best<- sep_cma_es(contextVector[groupMember],
                          fn = subfunctionCMA,
                          contextVector = contextVector,
                          groupMember = groupMember,mainfun=fun,...,
                          lower = lbound[groupMember],
                          upper=ubound[groupMember],
                          #control = list(mu=20,lambda=20,maxit=2000))
                          control = list(mu=10,
                                         lambda=10,
                                         maxit=1500,
                                         sigma=0.3*max(ubound[groupMember]-lbound[groupMember])))
        nEval <- nEval + best$counts[1]

        print('updating context vector (separable) for interconnection step...')
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          if(!is.null(best$par)){
            contextVector[groupMember] <- best$par
            obj <- best$value
            if(obj < bestObj){
              bestPop <- contextVector
              bestObj <- obj
              print('Update:')
              print(bestObj)
            }
          }else{
            print('is null')
          }
        }else{
          break
        }


        best <- cma_es(par = contextVector,
                       fn = subfunctionCMA,
                       contextVector = contextVector,
                       groupMember = 1:nVar,
                       mainfun=fun,...,
                       lower = lbound,
                       upper=ubound,
                       #control = list(vectorized=T,maxit=1000,mu=20,lambda=20))
                       control = list(vectorized=T,
                                      maxit=150,
                                      mu=10,lambda=10,
                                      sigma=0.3*max(ubound-lbound)))
        nEval <- nEval + best$counts[1]
        print('Interconnection step (separable) finished, updating context vector...')
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          if(!is.null(best$par)){
            contextVector <- best$par
            obj <- best$value
            if(obj < bestObj){
              bestPop <- contextVector
              bestObj <- obj
              print('Update:')
              print(bestObj)
            }
          }else{
            print('is null')
          }
        }else{
          break
        }
      }
      # optimize non-separable
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
                        # control = list(mu=20,lambda=20,maxit=2000))
                        control = list(mu=10,
                                       lambda=10,
                                       maxit=1500,
                                       sigma=0.3*max(ubound[groupMember]-lbound[groupMember])))
          nEval <- nEval + best$counts[1]

          print('updating context vector (non-separable) for interconnection step...')
          if((budget-nEval)>0){ # only update if it doesnt exceed budget
            if(!is.null(best$par)){
              contextVector[groupMember] <- best$par
              obj <- best$value
              if(obj < bestObj){
                bestPop <- contextVector
                bestObj <- obj
                print('Update:')
                print(bestObj)
              }
            }else{
              print('is null')
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
                                      maxit=150,mu=10,
                                      lambda=10,
                                      sigma=0.3*max(ubound-lbound)))
        nEval <- nEval + best$counts[1]

        print('Interconnection step (nonseparable) finished, updating context vector...')
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          if(!is.null(best$par)){
            contextVector <- best$par
            obj <- best$value
            if(obj < bestObj){
              bestPop <- contextVector
              bestObj <- obj
              print('Update:')
              print(bestObj)
            }
          }else{
            print('is null')
          }
        }else{
          break
        }
        leftBudget <- budget - nEval
        print(c('Comp budget left:',leftBudget,budget,nEval))
        save(list=ls(),file=paste('dataMOFBVE_grouped_',seed,'.Rdata',sep=''))
      }
    }
  }
  return(list(x=bestPop,y=bestObj))
}
