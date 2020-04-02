#' Main function for the cooperative coevolution with multilevel framework.
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
#' ML_LSGO(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#'
ML_LSGO <- function(nVar,fun,...,
                    group=NULL,variable_importance=NULL,
                    budget=3000000,
                    lbound=NULL,ubound=NULL,
                    nLevel=4,eval_interval=100000,
                    SA_method=c('morris_mu','morris_k','rf','sobol'),
                    maxSeqIter=0,
                    fname_prefix='datalogf1_',infill=c("best","ei"),mode='optim',repair=F){
  infill <- infill[1]
  bestPop <- NULL

  if(is.null(lbound))
    lbound <- rep(-Inf,nVar)
  if(is.null(ubound))
    ubound <- rep(Inf,nVar)
  SA_method <- SA_method[1]
  convergence_history <- NULL
  vars <- 1:nVar
  prevLevel <- NULL
  groupWeight <- NULL

  print(group)
  print('grouping...')
  if(is.null(group)){
    if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
      r <- 20
      a <- sensitivity::morris(model=fun,
                               factor=nVar,
                               r = r,
                               design = list(type='oat',levels=10,grid.jump=4),
                               binf=lbound,
                               bsup=ubound,
                               scale=T,...)

      bestPopIndex <- which.min(a$y)
      bestPop <- a$X[bestPopIndex,]
      bestObj <- min(a$y)

      cv <<- bestPop
      nEval <<- nEval + r*(nVar+1)
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
        clusterOrder <- nLevel:1
        group <- group[clusterOrder]
      }
    }
    if(SA_method == 'rf') {
      # cl <- makeCluster(10, type="SOCK")
      # dopara::registerDoSNOW(cl)
      doParallel::registerDoParallel(3)
      nEval <<- nEval + 20000

      population <- lhs::randomLHS(nEval,nVar)*(ubound-lbound)+lbound
      objectiveValue <- fun(population,...)
      bestPopIndex <- which.min(objectiveValue)
      bestPop <- population[bestPopIndex,]
      bestObj <- min(objectiveValue)
      currentPop <- list(objectiveValue=objectiveValue,population=population)
      #cv <<- bestPop
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
      nEval <<- nEval + 2*nEvalGroup*(nVar)

      a <- sensitivity::soboljansen(model=fun,
                                    X1=data.frame(population1),
                                    X2=data.frame(population2),
                                    nboot = 0,...)

      bestPopIndex <- which.min(a$y)
      bestPop <- data.matrix(a$X[bestPopIndex,])
      bestObj <- min(a$y)

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
    groupPortion <- 1 + floor((log(groupWeight)-log(minWeight)))
    totalPortion <- sum(groupPortion)
  }else{
    if(is.null(variable_importance)){
      reg_coefficient <- rep(1,nVar) # set same coefficient for all
    }else{
      reg_coefficient <- variable_importance # copy from argument
    }
    ## determine group weight from variable importance
    for(i in 1:nLevel){
      groupMember <- group[[i]]
      groupWeight <- append(groupWeight,(sum(reg_coefficient[groupMember])))
    }
    ## determine computational resource portion from group weight
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

  if(is.null(bestPop)){
    bestPop <- runif(nVar)*(-lbound+ubound)+lbound
    bestObj <- fun(bestPop,...)
    nEval <- nEval + 1
  }
  cv <<- bestPop

  print('start optim...')
  print(group)
  budget_left <- budget-nEval
  while(budget_left > 0){
    nGroup <- length(group)
    infill <- infill[1]

    nVar <- length(cv)
    groupLength <- pracma::zeros(1,nGroup)
    for(i in 1:nGroup){
      groupLength[i] <- length(group[[i]])
    }

    # build and collect surrogates
    # by optimizing lower levels
    budget_left <- budget-nEval
    maxit_lower <- 250
    maxit_upper <- 10

    cycle <- 0
    for(groupIndex in 1:nLevel){
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]
      currentGroupPortion <- groupPortion[groupIndex]

      CMAES_control[[groupIndex]] <<- list(vectorized=F,
                                           mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                                           sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                           diag.value=T,diag.pop=T,noImprove=25)
    }

    while(budget_left > 0){
      for(levelIndex in nGroup:2){
        print(paste('level',levelIndex))

        valCollection <- NULL
        upperLevel <- levelIndex - 1
        upperGroup <- group[[upperLevel]]
        currentGroup <- group[[levelIndex]]
        currentVector <- cv

        # optimize upper level
        # inside, optimize the lower level for each suggested upper level vector
        lambda_upper <- CMAES_control[[levelIndex]]$lambda
        CMAES_control[[upperLevel]]$maxit <<- maxit_upper
        CMAES_control[[levelIndex]]$maxit <<- maxit_lower
        CMAES_control[[upperLevel]]$lambda <<- CMAES_control[[upperLevel]]$mu

        bestval_lower <<- bestval
        funValue <- cma_es(par=currentVector[upperGroup],
                           fn = EvaluateUpperLevel_bilevel,
                           upperGroup=upperGroup,
                           lower = lbound[upperGroup],
                           upper = ubound[upperGroup],
                           fun=fun,
                           levelIndex = upperLevel,
                           lowerGroup = group[[levelIndex]],
                           lbound=lbound,
                           ubound=ubound,
                           control=CMAES_control[[upperLevel]],
                           ...)


        CMAES_control[[upperLevel]]$cov <<- funValue$cov
        CMAES_control[[upperLevel]]$sigma <<- funValue$sigma
        CMAES_control[[upperLevel]]$pc <<- funValue$pc
        CMAES_control[[upperLevel]]$ps <<- funValue$ps

        CMAES_control[[upperLevel]]$lambda <<- lambda_upper

        print(c(funValue$value,bestObj))
        if(funValue$value < bestObj){
          print('update best obj')
          cv[upperGroup] <<- funValue$par
          bestObj <- funValue$value
          bestval <<- bestObj
        }
        print(nEval)
        print(bestObj)
      }
      # print('level 1')
      # currentGroup <- group[[1]]
      # currentVector <- cv
      # CMAES_control[[1]]$maxit <<- maxit_lower
      # CMAES_control[[1]]$vectorized <<- T
      # funValue <- cma_es(par=currentVector[currentGroup],
      #                    fn=subfunctionCMA,
      #                    contextVector=cv,
      #                    groupMember=currentGroup,
      #                    mainfun=fun,...,
      #                    lower = lbound[currentGroup],
      #                    upper = ubound[currentGroup],
      #                    control=CMAES_control[[1]],repair=repair)
      # CMAES_control[[1]]$vectorized <<- F
      # nEval <<- nEval + funValue$counts[1]
      # CMAES_control[[1]]$cov <<- funValue$cov
      # CMAES_control[[1]]$sigma <<- funValue$sigma
      # CMAES_control[[1]]$pc <<- funValue$pc
      # CMAES_control[[1]]$ps <<- funValue$ps


      print(c(funValue$value,bestObj))
      if(funValue$value < bestObj){
        cv[currentGroup] <<- funValue$par
        bestObj <- funValue$value
      }

      budget_left <- budget-nEval
      cycle <- cycle + 1
    }
  }

  return(list(x=cv,y=bestObj))
}
