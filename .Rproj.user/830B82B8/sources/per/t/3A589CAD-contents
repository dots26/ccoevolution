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
ML_LSGO2 <- function(nVar,fun,...,
                     group=NULL,variable_importance=NULL,
                     budget=3000000,
                     lbound=NULL,ubound=NULL,
                     nLevel=4,eval_interval=100000,
                     grouping_method=c('SA','DG2','RDG3'),
                     SA_method=c('morris_mu','morris_k','rf','sobol'),
                     maxSeqIter=0,disableIPOP=F,maxit_upper=0,
                     fname_prefix='datalogf1_',infill=c("best","ei"),mode='optim',repair=F){
  infill <- infill[1]
  nEval <- 0
  bestPop <- NULL
  learning_period <- 3000
  nEvalHist <- NULL
  bestObjHist <- NULL
  if(is.null(lbound))
    lbound <- rep(-Inf,nVar)
  if(is.null(ubound))
    ubound <- rep(Inf,nVar)
  SA_method <- SA_method[1]
  convergence_history <- NULL

  vars <- 1:nVar
  prevLevel <- NULL
  groupWeight <- NULL
  allcov <- diag(nVar)
  pc <- numeric(nVar)
  ps <- numeric(nVar)
  print(group)
  print('grouping...')
  if(is.null(group)){
    if(grouping_method=='SA'){
      if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
        r <- 20
        a <- sensitivity::morris(model=fun,
                                 factors=nVar,
                                 r = r,
                                 design = list(type='oat',levels=8,grid.jump=4),
                                 binf=lbound,
                                 bsup=ubound,
                                 scale=T,...)

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
          clusterOrder <- nLevel:1
          group <- group[clusterOrder]
          groupWeight <- groupWeight[clusterOrder]
          save(reg_coefficient,a,ranks,unorderedCurrentGroup,orderedGroup,groupWeight,file='rc2.Rdata')
        }
      }
      if(SA_method == 'rf') {
        doParallel::registerDoParallel(3)
        nEval <- nEval + 20000

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
        clusterOrder <- nLevel:1
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
      nGroup <- nLevel
    }else if(grouping_method=="RDG3"){
      reg_coefficient <- rep(1,nVar) # set same coefficient for all
      RDG3_result <- RDG3(fun = fun,...,
                          lbound =lbound,ubound=ubound)
      group <- RDG3_result$group
      group <- append(group, list(RDG3_result$separable))
      nEval <- nEval + RDG3_result$nEval
      bestPopIndex <- which.min(RDG3_result$y)
      bestPop <- RDG3_result$x[bestPopIndex,]
      bestObj <- min(RDG3_result$y)
      contextVector <- bestPop
      nGroup <- length(group)

      for(i in 1:nGroup){
        groupMember <- group[[i]]
        groupWeight <- append(groupWeight,(sum(reg_coefficient[groupMember])))
      }
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
    for(i in 1:nGroup){
      groupMember <- group[[i]]
      groupWeight <- append(groupWeight,(sum(reg_coefficient[groupMember])))
    }
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


  if(is.null(bestPop)){
    bestPop <- runif(nVar)*(-lbound+ubound)+lbound
    bestObj <- fun(bestPop,...)
    nEval <- nEval + 1
  }
  # cv <<- bestPop

  print('start optim...')
  print(group)
  budget_left <- budget-nEval
  while(budget_left > 0){
    nGroup <- length(group)
    infill <- infill[1]

    nVar <- length(contextVector)
    groupLength <- pracma::zeros(1,nGroup)
    for(i in 1:nGroup){
      groupLength[i] <- length(group[[i]])
    }

    budget_left <- budget-nEval


    cycle <- 0
    CMAES_control <- list()
    for(groupIndex in 1:nGroup){
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]
      currentGroupPortion <- groupPortion[groupIndex]
      maxit <- round(learning_period*(currentGroupPortion/totalPortion))*nGroup
      CMAES_control[[groupIndex]] <- list(vectorized=T,
                                          mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                                          sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                          diag.value=F,diag.pop=F,noImprove=25,maxit=maxit,disableIPOP=disableIPOP)
    }

    while(budget_left > 0){
      for(levelIndex in nGroup:1){
        budget_left <- budget-nEval
        if(budget_left<=0)
          break
        print(paste('level',levelIndex))
        currentGroupPortion <-  groupPortion[levelIndex]
        currentGroup <- group[[levelIndex]]
        contextVector <- bestPop
        message(paste0('optimizing group ',levelIndex,' ',currentGroupPortion/totalPortion*100,'% portion'))

        if(maxit_upper>0){
          if(levelIndex!=1)
          {
            upperLevel <- levelIndex - 1
            upperLength1 <- length(group[[upperLevel]])
            upperLength2 <-  length(group[[levelIndex]])
            upperGroup <- append(group[[upperLevel]],group[[levelIndex]])
            upper_s <- CMAES_control[[upperLevel]]$sigma
            lower_s <- CMAES_control[[levelIndex]]$sigma
            pc <- append(CMAES_control[[upperLevel]]$pc,CMAES_control[[levelIndex]]$pc)
            ss <- min(c(upper_s,lower_s))
            groupSize <- length(upperGroup)

            # upper level CMAES
            cov <- diag(upperLength1+upperLength2)
            if(!is.null( CMAES_control[[upperLevel]]$cov))
              cov[1:upperLength1,1:upperLength1] <- CMAES_control[[upperLevel]]$cov
            if(!is.null(CMAES_control[[levelIndex]]$cov))
              cov[upperLength1+(1:upperLength2),upperLength1+(1:upperLength2)] <- CMAES_control[[levelIndex]]$cov

            # make sure budget is not exceeded
            budget_left <- budget-nEval
            if(maxit_upper*CMAES_control[[upperLevel]]$lambda>budget_left){
              maxit_upper <- floor(budget_left/CMAES_control[[upperLevel]]$lambda)
              if(maxit_upper<=0){
                nEval <- budget
                break
              }
            }

            CMAES_c <- list(vectorized=T, cov=cov, pc=pc,
                            mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                            sigma=ss,
                            diag.value=F,diag.pop=F,noImprove=25,maxit=maxit_upper)

            funValue <- cma_es(par=contextVector[upperGroup],
                               fn=subfunctionCMA,
                               contextVector=contextVector,
                               groupMember=upperGroup,
                               mainfun=fun,...,
                               lower = lbound[upperGroup],
                               upper = ubound[upperGroup],
                               disableIPOP=disableIPOP,
                               control=CMAES_c,
                               repair=repair)

            nEval <- nEval + funValue$counts[1]
            print(nEval)
            CMAES_control[[upperLevel]]$cov <- funValue$cov[1:upperLength1,1:upperLength1]
            CMAES_control[[upperLevel]]$pc <- funValue$pc[1:upperLength1]
            CMAES_control[[upperLevel]]$sigma <- funValue$sigma
            # CMAES_control[[upperLevel]]$ps <- funValue$ps[1:upperLength1]

            CMAES_control[[levelIndex]]$cov <- funValue$cov[upperLength1+(1:upperLength2),upperLength1+(1:upperLength2)]
            CMAES_control[[levelIndex]]$pc <- funValue$pc[upperLength1+(1:upperLength2)]
            CMAES_control[[levelIndex]]$sigma <- funValue$sigma
            # CMAES_control[[levelIndex]]$ps <- funValue$ps[upperLength1+(1:upperLength2)]

            if(funValue$value < bestObj){
              print('update best obj')
              print(c(funValue$value,bestObj))
              contextVector[upperGroup] <- funValue$par
              bestObj <- funValue$value
              nEvalHist <- append(nEvalHist,nEval)
              bestObjHist <- append(bestObjHist,bestObj)
              save(nEvalHist,bestObjHist,contextVector,file='datanoupper.Rdata')
            }
          }
        }

        # lower level optimization
        budget_left <- budget-nEval
        if(CMAES_control[[levelIndex]]$maxit*CMAES_control[[levelIndex]]$lambda>budget_left){
          CMAES_control[[levelIndex]]$maxit <- floor(budget_left/CMAES_control[[levelIndex]]$lambda)
          if(CMAES_control[[levelIndex]]$maxit<=0){
            nEval <- budget
            break
          }
        }


        CMAES_c <- CMAES_control[[levelIndex]]

        funValue <- cma_es(par=contextVector[currentGroup],
                           fn=subfunctionCMA,
                           contextVector=contextVector,
                           groupMember=currentGroup,
                           mainfun=fun,...,
                           lower = lbound[currentGroup],
                           upper = ubound[currentGroup],
                           disableIPOP=disableIPOP,
                           control=CMAES_c,
                           repair=repair,
                           inputScaleFactor=1,
                           inputScaleShift=0)
        # allcov[currentGroup,currentGroup] <- funValue$cov
        CMAES_control[[levelIndex]]$cov <- funValue$cov
        CMAES_control[[levelIndex]]$sigma <- funValue$sigma
        CMAES_control[[levelIndex]]$ps <- funValue$ps
        CMAES_control[[levelIndex]]$pc <- funValue$pc


        nEval <- nEval + funValue$counts[1]
        print(nEval)

        budget_left <- budget-nEval
        if(budget_left<=10)
          break
        if(!is.null(funValue$par)){
          contextVector[currentGroup] <- funValue$par

          if(funValue$value < bestObj){
            print('update best obj')
            print(c(funValue$value,bestObj))
            bestObj <- funValue$value
            bestPop <- contextVector
            nEvalHist <- append(nEvalHist,nEval)
            bestObjHist <- append(bestObjHist,bestObj)
            save(nEvalHist,bestObjHist,contextVector,file='datanoupper.Rdata')
          }
        }
      }

      budget_left <- budget-nEval
      if(budget_left<=10)
        break
      cycle <- cycle + 1
    }
  }

  return(list(x=bestPop,y=bestObj))
}
