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
#' SACC_DE(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
SACC_DE <- function(contextVector=NULL,nVar,fun,...,
                 group=NULL,
                 budget=3000000,
                 lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),
                 nLevel=4,evalInterval=100000,
                 SA_method=c('morris_mu','morris_k','rf','sobol')){

  SA_method <- SA_method[1]
  nEval <- 0
  convergence_history <- NULL
  vars <- 1:nVar
  prevLevel <- NULL
  groupWeight <- NULL
  #group <- NULL

  if(is.null(group)){
    if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
      r <- 20
      a <- sensitivity::morris(model=fun,
                               factor=nVar,
                               r = r,
                               design = list(type='oat',levels=8,grid.jump=4),
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

      if(SA_method == 'morris_k'){
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
      }
      if(SA_method == 'morris_mu'){
        ranks <- order(mu.star,decreasing = T)

        reg_coefficient <- mu.star
        totalCoef <- sum(reg_coefficient)

        rankedVars <- vars[ranks]
        nLevelVar <- floor(nVar/nLevel)
        for(groupingIndex in 1:(nLevel-1)){
          start <- (groupingIndex-1)*nLevelVar+1
          end <- groupingIndex*nLevelVar
          unorderedCurrentGroup <- rankedVars[start:end]
          orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
          group <- append(group,list(orderedGroup))
          groupWeight <- append(groupWeight,sum(reg_coefficient[orderedGroup]))
        }
        start <- (nLevel-1)*nLevelVar+1
        end <- nVar
        unorderedCurrentGroup <- rankedVars[start:end]
        orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
        groupWeight <- append(groupWeight,sum(reg_coefficient[orderedGroup]))
        group <- append(group,list(orderedGroup))
        clusterOrder <- 1:nLevel
      }
    }
    if(SA_method == 'rf') {
      # cl <- makeCluster(10, type="SOCK")
      # dopara::registerDoSNOW(cl)
      doParallel::registerDoParallel(3)
      nEval <- nEval + 20000

      population <- lhs::randomLHS(nEval,nVar)*(ubound-lbound)+lbound
      objectiveValue <- fun(population,...)
      bestPopIndex <- which.min(objectiveValue)
      bestPop <- population[bestPopIndex,]
      bestObj <- min(objectiveValue)
      currentPop <- list(objectiveValue=objectiveValue,population=population)
      contextVector <- bestPop
      prevLevel <- NULL

      group <- NULL

      nWorker <- foreach::getDoParWorkers()
      treePerCluster <- ceiling(500/nWorker)
      rf <- foreach::foreach(ntree=rep(treePerCluster, nWorker),
                             .combine=randomForest::combine,
                             .packages=c('randomForest')) %dopar% {
                               randomForest::randomForest(population, objectiveValue, ntree=ntree,mtry=100)
                             }

      impo <- abs(rf$importance)
      totalCoef <- sum(impo)
      ranks <- order(impo,decreasing = T)
      rankedVars <- vars[ranks]
      nLevelVar <- floor(nVar/nLevel)
      for(groupingIndex in 1:(nLevel-1)){
        start <- (groupingIndex-1)*nLevelVar+1
        end <- groupingIndex*nLevelVar
        unorderedCurrentGroup <- rankedVars[start:end]
        orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
        groupWeight <- append(groupWeight,sum(impo[orderedGroup]))
        group <- append(group,list(orderedGroup))
      }
      start <- (nLevel-1)*nLevelVar+1
      end <- nVar
      unorderedCurrentGroup <- rankedVars[start:end]
      orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
      groupWeight <- append(groupWeight,sum(impo[orderedGroup]))
      group <- append(group,list(orderedGroup))
      clusterOrder <- 1:nLevel
    }
    if(SA_method == 'sobol'){
      nEvalGroup <- 10
      population1 <- lhs::randomLHS(nEvalGroup,nVar)*(ubound-lbound)+lbound
      population2 <- lhs::randomLHS(nEvalGroup,nVar)*(ubound-lbound)+lbound
      nEval <- nEval + 2*nEvalGroup*(nVar)

      a <- sensitivity::soboljansen(model=fun,
                                    X1=data.frame(population1),
                                    X2=data.frame(population2),
                                    nboot = 0,...)

      bestPopIndex <- which.min(a$y)
      bestPop <- data.matrix(a$X[bestPopIndex,])
      bestObj <- min(a$y)

      contextVector <- bestPop
      prevLevel <- NULL

      sobolTotalIndex <- abs(data.matrix(a$T))
      ranks <- order(sobolTotalIndex,decreasing = T)

      reg_coefficient <- sobolTotalIndex
      totalCoef <- sum(reg_coefficient)

      rankedVars <- vars[ranks]
      nLevelVar <- floor(nVar/nLevel)
      for(groupingIndex in 1:(nLevel-1)){
        start <- (groupingIndex-1)*nLevelVar+1
        end <- groupingIndex*nLevelVar
        unorderedCurrentGroup <- rankedVars[start:end]
        orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
        groupWeight <- append(groupWeight,sum(reg_coefficient[orderedGroup]))
        group <- append(group,list(orderedGroup))
      }
      start <- (nLevel-1)*nLevelVar+1
      end <- nVar
      unorderedCurrentGroup <- rankedVars[start:end]
      orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
      groupWeight <- append(groupWeight,(sum(reg_coefficient[orderedGroup])))
      group <- append(group,list(orderedGroup))
      clusterOrder <- 1:nLevel
    }

    # checkMachineLimit <- log(groupWeight/minWeight)
    # while(is.infinite(checkMachineLimit)

    minWeight <- min(groupWeight)
    groupPortion <- vector(length=length(groupWeight))
    for(groupIndex in 1:length(groupWeight)){
      if(log(groupWeight[groupIndex])>0){
        groupPortion[groupIndex] <- 1 + log(groupWeight[groupIndex])
      }else{
        groupPortion[groupIndex] <- 1
      }
    }
    totalPortion <- sum(groupPortion)
    print(groupPortion)
  }

  CMAES_control <- list()
  DE_control <- list()
  pop <- list()

  for(groupIndex in 1:nLevel){
    groupSize <- length(group[[groupIndex]])
    groupMember <- group[[groupIndex]]
    currentGroupPortion <- groupPortion[groupIndex]
    CMAES_control[[groupIndex]] <- list(vectorized=T,
                                        mu=groupSize,lambda=groupSize,
                                        maxit=2000*round(currentGroupPortion/totalPortion),
                                        sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                        diag.value=T)
    DE_control[[groupIndex]] <- list(ccm=0.5,itermax=round(currentGroupPortion))
    pop[[groupIndex]] <- t(InitializePopulationLHS(NP,groupSize,lbound[groupMember],ubound[groupMember]))

  }
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')

  leftBudget <- budget - nEval
  nRepeat <- floor(nEval/evalInterval)
  for(i in 1:nRepeat){
    convergence_history <- append(convergence_history,bestObj)
  }
  while((budget-nEval)>0){
    for(groupIndex in 1:length(group)) {
      groupSize <- length(group[[groupIndex]])
      NP <- 50
      groupMember <- group[[groupIndex]]
      currentGroupPortion <- groupPortion[groupIndex]
      #print(paste0('optimizing group ',groupIndex,' with ',groupSize,' members, ',currentGroupPortion/totalPortion*100,'% portion'))
      # group optimization
      best <-  sansde(pop=pop[[groupIndex]],
                      bestmem=bestPop[groupMember],
                      bestval=bestObj,
                      fname = subfunction,
                      contextVector = contextVector,
                      groupMember = groupMember,mainfun=fun,...,
                      Lbound = lbound[groupMember],
                      Ubound =ubound[groupMember],
                      control = DE_control[[groupIndex]])
      pop[[groupIndex]] <- best$pop
      DE_control[[groupIndex]]$ccm <- best$ccm

      nlogging_this_layer <- floor((nEval+best$used_FEs)/evalInterval)-floor(nEval/evalInterval)

      if(nlogging_this_layer>0){
        for(i in 1:nlogging_this_layer){
          nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
          nGeneration_to_consider <- floor(nEval_to_logging/NP)
          bestObj_logging <- min(best$tracerst)
          convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
        }
      }
      nEval <- nEval + best$used_FEs

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
    }
    # Interconnection
    # groupMember <- 1:nVar
    # groupSize <- nVar
    # NP <- 50
    # pop <- t(InitializePopulationLHS(NP,groupSize,lbound[groupMember],ubound[groupMember]))
    # best <-  sansde(pop=pop,
    #                 bestmem=contextVector,
    #                 bestval=bestObj,
    #                 fname = subfunction,
    #                 contextVector = contextVector,
    #                 groupMember = groupMember,mainfun=fun,...,
    #                 Lbound = lbound[groupMember],
    #                 Ubound =ubound[groupMember],
    #                 control = list(itermax=1,ccm=0.5))
    # nlogging_this_layer <- floor((nEval+best$used_FEs)/evalInterval)-floor(nEval/evalInterval)
    #
    # if(nlogging_this_layer>0){
    #   for(i in 1:nlogging_this_layer){
    #     nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
    #     nGeneration_to_consider <- floor(nEval_to_logging/NP)
    #     bestObj_logging <- min(best$tracerst)
    #     convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
    #   }
    # }
    # nEval <- nEval + best$used_FEs
    #
    # if((budget-nEval)>0){ # only update if it doesnt exceed budget
    #   if(!is.null(best$par)){
    #     contextVector <- best$par
    #     obj <- best$value
    #     if(obj < bestObj){
    #       bestPop <- contextVector
    #       bestObj <- obj
    #     }
    #   }
    # }else{
    #   break
    # }
    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget,budget,nEval))
    print(bestObj)
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
