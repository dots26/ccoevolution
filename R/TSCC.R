#' Two stage cooperative coevolution with CMA-ES
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
TSCC <- function(contextVector=NULL,nVar,
                 fun,...,
                 group=NULL,
                 budget=3000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),
                 nLevel=4,evalInterval=100000,
                 SA_method=c('morris_mu','morris_k','rf','sobol'),
                 keepCovariance=F){
  SA_method <- SA_method[1]
  minGroupSize <- 50
  groupWeight <- NULL
  nEval <- 0
  convergence_history <- NULL
  group <- NULL
  vars <- 1:nVar
  #group <- NULL
  if(is.null(group)){
    if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
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
    minWeight <- min(groupWeight)
    groupPortion <- 1 + log(groupWeight)-log(minWeight)
    totalPortion <- sum(groupPortion)
  }

  print('ordered group')
  print(group)

  print('Secondary Grouping...')
  new_group <- NULL
  dg <- NULL
  for(i in 1:nLevel){
    print(paste0('Level ',i))
    gctrl <- list(lbound=lbound[group[[i]]],ubound=ubound[group[[i]]])
    subgroup <- DG2(length(group[[i]]),
                    subfunction,
                    control=gctrl,
                    contextVector=contextVector,
                    groupMember=group[[i]],
                    mainfun=fun,...)

    dg <- cbind(dg,subgroup)
    nEval <- nEval + subgroup$nEval
  }
  dg<- dg[,clusterOrder]
  nlogging_this_layer <- floor(nEval/evalInterval)
  if(nlogging_this_layer>0){
    for(i in 1:nlogging_this_layer){
      nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
      convergence_history <- append(convergence_history,bestObj)
      # print(c('conv',convergence_history))
    }
  }
  leftBudget <- budget - nEval

  CMAES_control <- vector(mode = "list",length=nLevel)
  for(clusterIndex in 1:nLevel){
    cluster_member <- group[[clusterIndex]]
    currentClusterGrouping <- dg[,clusterIndex]$group
    sep <- group[[clusterIndex]][dg[,clusterIndex]$separable]
    groupSize <- length(sep)
    groupMember <- sep
    currentGroupPortion <- groupPortion[[clusterIndex]]
    CMAES_control[[clusterIndex]]$sep <- list(vectorized=T,
                                            mu=groupSize,lambda=groupSize,
                                            maxit=round(3600*currentGroupPortion/totalPortion),
                                            sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                            diag.value=T)
    CMAES_control[[clusterIndex]]$nonsep <- list()

    if(length(currentClusterGrouping)>0){
      for(groupIndex in 1:length(currentClusterGrouping)) {
        groupMember <- cluster_member[currentClusterGrouping[[groupIndex]]]
        groupSize <- length(groupMember)
        print('g')
        CMAES_control[[clusterIndex]]$nonsep[[groupIndex]] <- list(mu=groupSize,lambda=groupSize,
                                                                 maxit=round(3600*currentGroupPortion/totalPortion),
                                                                 sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                                                 diag.value=T)
      }
    }
  }

  nRepeat <- floor(nEval/evalInterval)
  for(i in 1:nRepeat){
    convergence_history <- append(convergence_history,bestObj)
  }
  #saved
  print(nEval)
  while((budget-nEval)>0 ){
    for(clusterIndex in 1:nLevel){
      print(paste('Optimizing Variable level',clusterIndex))
      cluster_member <- group[[clusterIndex]]
      currentClusterGrouping <- dg[,clusterIndex]$group
      sep <- group[[clusterIndex]][dg[,clusterIndex]$separable]

      currentGroupPortion <- groupPortion[[clusterIndex]]

      # optimize separable
      groupMember <- sep
      groupSize <- length(groupMember)

      if(groupSize>0){
        print(c('optimizing separable variables'))
        # group optimization
        best<- sep_cma_es(contextVector[groupMember],
                          fn = subfunctionCMA,
                          contextVector = contextVector,
                          groupMember = groupMember,mainfun=fun,...,
                          lower = lbound[groupMember],
                          upper=ubound[groupMember],
                          control = CMAES_control[[clusterIndex]]$sep)
        if(keepCovariance){
          CMAES_control[[clusterIndex]]$sep$cov <- best$cov
          CMAES_control[[clusterIndex]]$sep$sigma <- best$sigma * 2
        }
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

      # optimize non-separable
      print('Optimize current group non-seps...')
      if(length(currentClusterGrouping)>0){
        for(groupIndex in 1:length(currentClusterGrouping)) {
          print(paste0('subgroup ',groupIndex))
          groupMember <- cluster_member[currentClusterGrouping[[groupIndex]]]
          groupSize <- length(groupMember)

          best<- cma_es(contextVector[groupMember],
                        fn = subfunctionCMA,
                        contextVector = contextVector,
                        groupMember = groupMember,
                        mainfun=fun,...,
                        lower = lbound[groupMember],
                        upper=ubound[groupMember],
                        # control = list(mu=groupSize,lambda=groupSize,maxit=2000))
                        control = CMAES_control[[clusterIndex]]$nonsep[[groupIndex]])
          if(keepCovariance){
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$cov <- best$cov
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma <- best$sigma * 2
          }
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

          print('updating context vector...')
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
      }
      leftBudget <- budget - nEval
      print(c('Comp budget left:',leftBudget,budget,nEval))
    }
    # Interconnection
    print('Interconnection')
    groupMember <- 1:nVar
    best <- cma_es(par = contextVector[groupMember],
                   fn = subfunctionCMA,
                   contextVector = contextVector,
                   groupMember = groupMember,
                   mainfun=fun,...,
                   lower = lbound,
                   upper=ubound,
                   # control = list(vectorized=T,maxit=1000,mu=20,lambda=20))
                   control = list(vectorized=T,
                                  maxit=90,
                                  mu=100,lambda=100,
                                  sigma=0.3*max(ubound-lbound),
                                  diag.value=T))
    nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
    if(nlogging_this_layer>0){
      for(i in 1:nlogging_this_layer){
        nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
        nGeneration_to_consider <- floor(nEval_to_logging/100)
        print(c(nEval_to_logging,groupSize,nGeneration_to_consider,length(best$diagnostic$value)))
        bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
        convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
        # print(convergence_history)
      }
    }
    nEval <- nEval + best$counts[1]

    if((budget-nEval)>0){ # only update if it doesnt exceed budget
      if(!is.null(best$par)){
        contextVector <- best$par
        obj <- best$value
        if(obj < bestObj){
          bestPop <- contextVector
          bestObj <- obj
        }
      }else{

      }
    }else{
      break
    }
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
