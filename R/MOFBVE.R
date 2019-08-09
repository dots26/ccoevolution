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
#' MOFBVE(nVar = 1000,fun=func,phaseSolver = cmaes,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
MOFBVE <- function(contextVector=NULL,nVar,fun,group=NULL,phaseSolver=cmaes,budget=1000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),nLevel=4,...){
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

  #  if(is.null(group)){

  print('sens')
  r<- 2
  a <- sensitivity::morris(model=fun,
                           factor=nVar,
                           r = r,
                           design = list(type='oat',levels=8,grid.jump=1),
                           binf=lbound,
                           bsup=ubound,
                           scale=F,...)
  bestPopIndex <- which.min(a$y)
  bestPop <- a$x[bestPopIndex,]
  bestObj <- min(a$y)
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
  save(group,file='MOFBVE_group.Rdata')
  # }
  print('Groups assigned')
  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')
  leftBudget <- budget - nEval
  #saved
  while((budget-nEval)>0 ){
    for(groupIndex in clusterOrder[1:length(group)]) {
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]

      # group optimization
      best<- cmaes::cma_es(contextVector[groupMember],
                           fn = subfunction,
                           contextVector = contextVector,
                           groupMember = groupMember,mainfun=fun,...,
                           lower = lbound[groupMember],
                           upper=ubound[groupMember],
                           control = list(mu=groupSize,lambda=groupSize,maxit=10*groupSize))
      nEval <- nEval + best$counts[1]

      print('updating context vector for interconnection step...')
      contextVector[groupMember] <- best$par
      obj <- fun(contextVector,...)
      nEval <- nEval + 1
      if(obj < bestObj){
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          bestPop <- contextVector
          bestObj <- obj
        }
      }

      # Interconnection
      best <- cmaes::cma_es(contextVector,fn = fun,...,lower = lbound,upper=ubound,control = list(mu=groupSize,lambda=groupSize,maxit=5*groupSize))
      nEval <- nEval + best$counts[1]
      print('Interconnection step finished, updating context vector...')
      contextVector <- best$par
      obj <- fun(contextVector,...)
      nEval <- nEval + 1
      if(obj < bestObj){
        if((budget-nEval)>0){ # only update if it doesnt exceed budget
          bestPop <- contextVector
          bestObj <- obj
        }
      }
    }
    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget,budget,nEval))
    #save(list=ls(),file=paste('dataMOFBVE_sur_',seed,'.Rdata',sep=''))
  }
  return(list(x=bestPop,y=bestObj))
}
