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
                 group=NULL, ## ignored
                 budget=3000000,lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),
                 nLevel=4,evalInterval=100000,
                 SA_method=c('morris_mu','morris_k','rf','sobol'),
                 keepCovariance=F,limit_sigma=F){

  ########### initialization ###############
  SA_method <- SA_method[1]
  groupWeight <- NULL
  nEval <- 0
  convergence_history <- NULL
  group <- NULL
  vars <- 1:nVar

  ########### grouping ###############
  group <- NULL
  if(is.null(group)){ ## TODO: what should be done if group is supplied
    if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
      r<- 20
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
  }
  print('ordered group')
  print(group)
  print(c(groupPortion,totalPortion))
  print('Secondary Grouping...')
  new_group <- NULL
  dg <- NULL


  for(i in 1:nLevel){
    print(paste0('Level ',i))
    gctrl <- list(lbound=lbound[group[[i]]],ubound=ubound[group[[i]]])
    subgroup <- DG2(length(group[[i]]),
                    subfunction,
                    control=gctrl,
                    contextVector=lbound,
                    groupMember=group[[i]],
                    mainfun=fun,...)
    # join small groups here
    nSubgroup <- length(subgroup)
    subgroupSize <- NULL
    subgroupOrderedBySize <- subgroup
    for(subgroupIndex in 1:nSubgroup){
      subgroupSize <- append(subgroupSize,length(subgroup[[subgroupIndex]]))
    }
    groupOrder <- order(subgroupSize,decreasing = F)
    for(subgroupIndex in 1:nSubgroup){
      current <- groupOrder[subgroupIndex]
      subgroupOrderedBySize[[subgroupIndex]] <- subgroup[[current]]
    }
    dg <- cbind(dg,subgroup)
    print('groups')
    print(subgroup$group)
    print('separable')
    print(subgroup$separable)
    nEval <- nEval + subgroup$nEval
  }
  dg<- dg[,clusterOrder,drop=F]
  nlogging_this_layer <- floor(nEval/evalInterval)
  if(nlogging_this_layer>0){
    for(i in 1:nlogging_this_layer){
      nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
      convergence_history <- append(convergence_history,bestObj)
    }
  }
  leftBudget <- budget - nEval

  ########### initialization of CMAES control param ###############
  CMAES_control <- vector(mode = "list",length=nLevel)
  for(clusterIndex in 1:nLevel){
    cluster_member <- group[[clusterIndex]]
    currentClusterGrouping <- dg[,clusterIndex]$group
    sep <- group[[clusterIndex]][dg[,clusterIndex]$separable]
    groupSize <- length(sep)
    groupMember <- sep
    currentGroupPortion <- groupPortion[clusterIndex]

    if(groupSize>0){
      CMAES_control[[clusterIndex]]$sep <- list(vectorized=T,
                                                mu=groupSize,lambda=groupSize,
                                                maxit=round(2000*(currentGroupPortion/totalPortion)),
                                                sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                                diag.value=T)
    }
    CMAES_control[[clusterIndex]]$nonsep <- list()
    if(length(currentClusterGrouping)>0){
      for(groupIndex in 1:length(currentClusterGrouping)) {
        groupMember <- cluster_member[currentClusterGrouping[[groupIndex]]]
        groupSize <- length(groupMember)
        print(paste('level',clusterIndex,'group',groupIndex,round(500*(currentGroupPortion/totalPortion))))
        CMAES_control[[clusterIndex]]$nonsep[[groupIndex]] <- list(vectorized=T,mu=groupSize,lambda=groupSize,
                                                                   maxit=round(2000*(currentGroupPortion/totalPortion)),
                                                                   sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                                                   diag.value=T)
      }
    }
  }

  nRepeat <- floor(nEval/evalInterval)
  for(i in 1:nRepeat){
    convergence_history <- append(convergence_history,bestObj)
  }
  ########### end initialization& grouping  ###############

  ########### begin cooperative coevolution ###############
  while((budget-nEval)>0 ){
    for(clusterIndex in 1:nLevel){
      print(paste('Optimizing Variable level',clusterIndex))
      cluster_member <- group[[clusterIndex]]
      currentClusterGrouping <- dg[,clusterIndex]$group
      sep <- group[[clusterIndex]][dg[,clusterIndex]$separable]

      currentGroupPortion <- groupPortion[clusterIndex]

      # optimize separable
      groupMember <- sep
      groupSize <- length(groupMember)

      if(groupSize>0){
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
          CMAES_control[[clusterIndex]]$sep$sigma <- best$sigma
          CMAES_control[[clusterIndex]]$sep$ps <- best$ps
          CMAES_control[[clusterIndex]]$sep$pc <- best$pc
        }
        nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
        if(nlogging_this_layer>0){
          for(i in 1:nlogging_this_layer){
            nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
            nGeneration_to_consider <- floor(nEval_to_logging/CMAES_control[[clusterIndex]]$sep$mu)
            #print(nGeneration_to_consider,nlogging_this_layer)
            #print(best$diagnostic)
            if(!is.matrix(best$diagnostic$value)){
              best$diagnostic$value <- matrix(best$diagnostic$value)
            }
            bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
            convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
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

      if(length(currentClusterGrouping)>0){
        for(groupIndex in 1:length(currentClusterGrouping)) {
          groupMember <- cluster_member[currentClusterGrouping[[groupIndex]]]
          groupSize <- length(groupMember)
          # print(contextVector[groupMember])
          # print(groupMember)
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
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma <- best$sigma
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$ps <- best$ps
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$pc <- best$pc
          }

          nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
          if(nlogging_this_layer>0){
            for(i in 1:nlogging_this_layer){
              nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
              nGeneration_to_consider <- floor(nEval_to_logging/CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$mu)
              #print(nGeneration_to_consider,nlogging_this_layer)
              #print(best$diagnostic)
              if(!is.matrix(best$diagnostic$value)){
                best$diagnostic$value <- matrix(best$diagnostic$value)
              }
              bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
              convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))

            }
          }
          nEval <- nEval + best$counts[1]

          # print('curbest')
          # print(bestPop)
          # print(bestObj)
          # print('new')
          # print(best$par)
          # print(best$value)
          if((budget-nEval)>0){ # only update if it doesnt exceed budget
            if(!is.null(best$par)){
              contextVector[groupMember] <- best$par
              obj <- best$value
              if(obj < bestObj){
                bestPop <- contextVector
                bestObj <- obj
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
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
