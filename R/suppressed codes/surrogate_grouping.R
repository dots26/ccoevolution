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
#' surrogate_grouping(nVar = 1000,fun=func,phaseSolver = cmaes,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
surrogate_grouping <- function(contextVector=NULL,nVar,group=NULL,fun,phaseSolver=cmaes,budget=1000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),nLevel=4,...){
  doParallel::registerDoParallel()
  print(c('Ncores=',foreach::getDoParWorkers() ))
  nEval <- 0
  # if(is.null(contextVector)){
  #   nEval <- nEval + 10000
  #   population <- (randtoolbox::sobol(nEval,nVar,scrambling = 3))*(ubound-lbound)+lbound
  #   # Evaluate the whole population
  #   print('Evaluating initial population...')
  #   objectiveValue <- fun(population,...)
  #   bestPopIndex <- which.min(objectiveValue)
  #   bestPop <- population[bestPopIndex,]
  #   bestObj <- min(objectiveValue)
  #
  #   contextVector <- bestPop
  # }

  #if(is.null(group)){
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
  print(c('clusters',group))

  print('Grouping...')
  new_group <- NULL
  dg <- NULL
  for(i in 1:nLevel){
    gctrl <- list(lbound=lbound[group[[i]]],ubound=ubound[group[[i]]])
    dg <- cbind(dg,DG2(length(group[[i]]),
                       subfunction,
                       control=gctrl,
                       contextVector=contextVector,
                       groupMember=group[[i]],
                       mainfun=fun,...))
    nEval <- nEval + dg[,i]$nEval
    #newGroupLength <- length(dg$group)
    #for(j in 1:newGroupLength){
    #  new_group <- append(new_group,list(group[[i]][dg$group[[j]]]))
    #}
    #separable <- append(separable,list(group[[i]][dg$separable[[1]]]))
    print(c('nEval',nEval))
  }
  #group <- new_group
  # print(c('subclusters',group))

  print('Groups assigned')
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval

  while((budget-nEval)>0 ){
    for(clusterIndex in clusterOrder){
      print(paste('Optimizing cluster no.',clusterIndex))
      this.cluster <- group[[clusterIndex]]
      currentClusterGrouping <- dg[,clusterIndex]$group
      sep <- group[[clusterIndex]][dg[,clusterIndex]]$separable

      # optimize separable
      groupMember <- sep
      groupSize <- length(groupMember)
      pop <- InitializePopulationLHS(groupSize,200,minVal=lbound[groupMember],maxVal=ubound[groupMember])
      obj <- subfunction(pop,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...)
      nEval <- nEval + 200

      # optimize surrogate model
      nCMAESiter <- 20
      for(cmaesIter in 1:nCMAESiter){
        print(paste('building and optimizing RF model for group:',groupIndex,'iter',cmaesIter,'/',nCMAESiter))
        randomForestModel <- SPOT::buildRandomForest(pop,obj)
        best<- sep_cma_es(par=contextVector[groupMember],
                                 fn=modelPrediction,
                                 model=randomForestModel,
                                 upper= ubound[groupMember],
                                 lower= lbound[groupMember],
                                 control = list(vectorized=T,maxit=100,mu=20,lambda=20))
        if(!is.null(best$par)){
          pop <- rbind(pop,best$par)
          obj <- append(obj,subfunction(best$par,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...))
        }
      }
      rm(randomForestModel)
      nEval <- nEval + nCMAESiter

      print('updating context vector (separable) for interconnection step...')
      contextVector[groupMember] <- best$par
      this.obj <- fun(contextVector,...)
      nEval <- nEval + 1
      if(this.obj < bestObj){
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          bestPop <- contextVector
          bestObj <- this.obj
          print('Update:')
          print(bestPop)
          print(bestObj)
        }else{
          break
        }
      }

      # Interconnection
      best <- cmaes::cma_es(par = contextVector,
                            fn = fun,...,
                            lower = lbound,
                            upper=ubound,
                            control = list(vectorized=T,maxit=1000,mu=20,lambda=20))
      nEval <- nEval + best$counts[1]
      print('Interconnection (separable phase) step finished, updating context vector...')
      if(!is.null(best$par)){
        contextVector <- best$par
        this.obj <- fun(contextVector,...)
        nEval <- nEval + 1
        if(this.obj < bestObj){
          if((budget-nEval)>0){ # only update if it doesnt exceed budget
            bestPop <- contextVector
            bestObj <- this.obj
            print('Update:')
            print(bestPop)
            print(bestObj)
          }else{
            break
          }
        }
      }

      # optimize non-sep
      for(groupIndex in 1:length(currentClusterGrouping)) {
        groupMember <- this.cluster[currentClusterGrouping[[groupIndex]]]
        groupSize <- length(groupMember)

        pop <- InitializePopulationLHS(groupSize,200,minVal=lbound[groupMember],maxVal=ubound[groupMember])
        obj <- subfunction(pop,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...)
        nEval <- nEval + 200

        # optimize surrogate model
        nCMAESiter <- 20
        for(cmaesIter in 1:nCMAESiter){
          print(paste('building and optimizing RF model for group:',groupIndex,'iter',cmaesIter,'/',nCMAESiter))
          randomForestModel <- SPOT::buildRandomForest(pop,obj)
          best<- cmaes::cma_es(par=contextVector[groupMember],
                               fn=modelPrediction,
                               model=randomForestModel,
                               upper= ubound[groupMember],
                               lower= lbound[groupMember],
                               control = list(maxit=100,mu=20,lambda=20))
          if(!is.null(best$par)){
            pop <- rbind(pop,best$par)
            obj <- append(obj,subfunction(best$par,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...))
          }
        }
        rm(randomForestModel)
        nEval <- nEval + nCMAESiter

        print('updating context vector (non-sep) for interconnection step...')
        contextVector[groupMember] <- best$par
        this.obj <- fun(contextVector,...)
        nEval <- nEval + 1
        if(this.obj < bestObj){
          if((budget-nEval)>0){ # only update if it doesnt exceed budget
            bestPop <- contextVector
            bestObj <- this.obj
            print('Update:')
            print(bestPop)
            print(bestObj)
          }else{
            break
          }
        }

        # Interconnection
        best <- cmaes::cma_es(par = contextVector,fn = fun,...,lower = lbound,upper=ubound,control = list(maxit=500,mu=20,lambda=20))
        nEval <- nEval + best$counts[1]
        print('Interconnection (non-sep phase) step finished, updating context vector...')
        if(!is.null(best$par)){
          contextVector <- best$par
          this.obj <- fun(contextVector,...)
          nEval <- nEval + 1
          if(this.obj < bestObj){
            if((budget-nEval)>0){ # only update if it doesnt exceed budget
              bestPop <- contextVector
              bestObj <- this.obj
              print('Update:')
              print(bestPop)
              print(bestObj)
            }else{
              break
            }
          }
        }
        leftBudget <- budget - nEval
        print(c('Comp budget left:',leftBudget,budget,nEval))
        save(list=ls(),file=paste('dataMOFBVE_grouped_sur_',seed,'.Rdata',sep=''))
      }
    }
  }
  return(list(x=bestPop,y=bestObj))
}


