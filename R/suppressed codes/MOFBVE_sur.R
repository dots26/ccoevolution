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
  save(group,file='MOFBVE_sur_group.Rdata')
  #  }else{

  #  }
  print('Groups assigned')
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval
  print(c('start optim, budget:',leftBudget))

  while((budget-nEval)>0 ){
    for(groupIndex in clusterOrder[1:length(group)]) {
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]

      pop <- t(InitializePopulationLHS(200,groupSize,minVal=lbound[groupMember],maxVal=ubound[groupMember]))
      print('evaluating population')
      obj <- subfunction(pop,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...)
      nEval <- nEval + 200

      # optimize surrogate model
      nCMAESiter <- 20
      print(length(groupMember))
      for(cmaesIter in 1:nCMAESiter){
        print(paste('building RF model for group:',groupIndex,'iter',cmaesIter,'/',nCMAESiter))
        start_time <- Sys.time()
        randomForestModel <- SPOT::buildRandomForest(pop,obj)
        print( Sys.time()-start_time)
        print('optimizing...')
        start_time <- Sys.time()
        best<- cmaes::cma_es(par=contextVector[groupMember],
                             fn=modelPrediction,
                             model=randomForestModel,
                             upper= ubound[groupMember],
                             lower= lbound[groupMember],
                             control = list(vectorized=T,maxit=100,mu=20,lambda=20))
        print(Sys.time()-start_time)
        if(!is.null(best$par)){
          pop <- rbind(pop,best$par)
          obj <- append(obj,subfunction(best$par,contextVector=contextVector,groupMember=groupMember,mainfun=fun,...))
        }
      }
      rm(randomForestModel)
      nEval <- nEval + nCMAESiter

      print('updating context vector for interconnection step...')
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
                            control = list(vectorized=T,maxit=10,mu=20,lambda=20))
      nEval <- nEval + best$counts[1]
      print('Interconnection step finished, updating context vector...')
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
      save(list=ls(),file=paste('dataMOFBVE_sur_',seed,'.Rdata',sep=''))
    }
  }
  return(list(x=bestPop,y=bestObj))
}
