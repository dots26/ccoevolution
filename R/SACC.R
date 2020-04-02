#' Main function for the cooperative coevolution with sensitivity-analysis-based grouping
#' The CMA-ES is used as solver and its
#' step size is set based on  \code{lbound} and  \code{ubound}.
#' This means \code{lbound} and  \code{ubound}
#' must be finite.
#'
#' @title Cooperative coevolution with sensitivity-analysis-based grouping.
#' @param population Initial population. Row major order.
#' @param fun The objective function object to be solved
#' @param lbound Lower bound of the decision vector. A vector with finite values.
#' @param ubound upper bound of the decision vector. A vector with finite values.
#' @param group Vector of list. Each list contains members of the group. If not supplied, the grouping will use the chosen SA_method.
#' @param variable_importance A vector. Used if group is supplied. If NULL, all variables will have same portion.
#' @param nLevel Number of level for Morris method. Used if group is not supplied
#' @param SA_method sensitivity analysis method. It is recommended to use the morris method: \code{'morris_mu'}.
#' Although other SA_method should not cause an error, the package is mostly tested with \code{'morris_mu'}.
#' @param keepCovariance Set whether the covariance matrix is persistent througout the run.
#' If false, the Cov matrix will reset to identity matrix in each cycle.
#' @param evalInterval Interval for data logging.
#' @param scale If \code{TRUE}, the input vector will be scaled to [0,1] for processing,
#' otherwise the it will be in the range \code{[lbound,ubound]}
#' @param ... Further arguments passed to \code{fun}.
#' @examples
#' optimum <- rep(13,1000)
#' func <- f1cec
#' ctrl <- list(lbound=rep(-100,1000),ubound=rep(100,1000),delta=rep(20,1000))
#' SACC(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
SACC <- function(contextVector=NULL,nVar,fun,...,
                 group=NULL,variable_importance=NULL,
                 budget=3000000,
                 lbound=rep(0,nVar),ubound=rep(1,nVar),
                 nLevel=4,evalInterval=100000,
                 SA_method=c('morris_mu','morris_k','rf','sobol'),
                 keepCovariance = T,scale=T){
  ############ initialization ###############

  if(!is.null(contextVector))
    nVar <- length(contextVector)

  SA_method <- SA_method[1]
  nEval <- 0
  convergence_history <- NULL
  vars <- 1:nVar
  prevLevel <- NULL
  groupWeight <- NULL
  learning_period <- 500
  #group <- NULL

  ########### grouping ################
  if(is.null(group)){
    if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
      r <- 20
      print(fun)
      print(nVar)
      a <- sensitivity::morris(model=fun,
                               factors=nVar,
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
      save(a,file='ee.Rdata')
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

  if(is.null(contextVector)){
    contextVector <- matrix(runif(nVar)*(ubound-lbound)+lbound,nrow=1)
    bestObj <- fun(contextVector,...)
    nEval <- nEval + 1
  }

  ########### initialize CMAES control parameter ############
  CMAES_control <- list()
  for(groupIndex in 1:nLevel){
    groupSize <- length(group[[groupIndex]])
    groupMember <- group[[groupIndex]]
    currentGroupPortion <- groupPortion[groupIndex]
    CMAES_control[[groupIndex]] <- list(vectorized=T,
                                        # mu=groupSize,lambda=groupSize,
                                        mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                                        maxit=round(learning_period*(currentGroupPortion/totalPortion)),
                                        sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                        diag.value=T)
  }
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(!all(unlist(lapply(group,is.vector)))) stop('Sublist of group is of wrong mode, all of them should also be a vector')

  leftBudget <- budget - nEval
  nRepeat <- floor(nEval/evalInterval)
  if(nRepeat>0){
    for(i in 1:nRepeat){
      convergence_history <- append(convergence_history,bestObj)
    }
  }

  ######### end initialization and grouping ##################
  if(scale){
    scaling_shift <- lbound
    scaling_factor <- 1/(ubound-lbound)
    ubound <- rep(1,nVar)
    lbound <- rep(0,nVar)
  }else{
    scaling_shift <- rep(0,nVar)
    scaling_factor <- rep(1,nVar)
  }
  ######### begin cooperative coevoloution ####################
  while((budget-nEval)>0){
    term_code_total <- 0
    term_code_zero <- 0
    groupIndex <- 1
    while(groupIndex <= length(group)) {
      groupSize <- length(group[[groupIndex]])
      groupMember <- group[[groupIndex]]
      currentGroupPortion <- groupPortion[groupIndex]
      message(paste0('optimizing group ',groupIndex,' with ',groupSize,' members, ',currentGroupPortion/totalPortion*100,'% portion'))

      # group optimization
      # set.seed(100)
      # print(contextVector[groupMember])
      # print(groupMember)
      best<- cma_es(contextVector[groupMember],
                    fn = subfunctionCMA,
                    contextVector = contextVector,
                    groupMember = groupMember,
                    mainfun=fun,...,
                    lower = lbound[groupMember],
                    upper = ubound[groupMember],
                    control = CMAES_control[[groupIndex]],
                    inputScaleFactor=scaling_factor[groupMember],
                    inputScaleShift=scaling_shift[groupMember])

      termination_code <- best$termination_code
      term_code_total <- term_code_total + 1

      # termination codes:
      # 0: no error -> do nothing/update covariance matrix
      # 1: sigma divergence -> create smaller groups
      # 2: no improvement for some generations -> increase pop size
      # 3: sigma is too small to have effect on covariance matrix -> increase pop size
      # 4: sigma is too small to have effect on mean x -> increase pop size
      # 5: condition number too large -> increase pop size

      if(termination_code==0){
        term_code_zero <- term_code_zero + 1
        if(keepCovariance){
          CMAES_control[[groupIndex]]$cov <- best$cov
          CMAES_control[[groupIndex]]$sigma <- best$sigma
          CMAES_control[[groupIndex]]$ps <- best$ps
          CMAES_control[[groupIndex]]$pc <- best$pc
        }
      }
      lambda <- CMAES_control[[groupIndex]]$lambda
      if(F){
      if(termination_code==1){
        print('separating group')
        # print(paste0('index: ', groupIndex))
        if(groupSize>3){
          # create smaller groups!
          # modify group
          original_group_index <- groupIndex
          original_group <- group[[original_group_index]]
          original_group_size <- length(original_group) # number of variable in the group
          original_group_count <- length(group) # number of group in the cluster

          modified_group_set <- group # make a copy to do a shift for the groups
          coeff <- reg_coefficient[groupMember]
          coeff_order <- order(coeff,decreasing = T)
          first_group_indices <- coeff_order[1:round(original_group_size/2)]
          second_group_indices <- coeff_order[(round(original_group_size/2)+1):original_group_size]

          first_group <- groupMember[first_group_indices]
          second_group <- groupMember[second_group_indices]

          modified_group_set[[original_group_index]] <- first_group # top half of the group
          modified_group_set[[original_group_index+1]] <- second_group# bottom half of the group

          # create CMAES control for the new groups
          # copy the control for current cluster
          CMAES_control_tmp <- CMAES_control
          CMAES_control_tmp[[original_group_index]] <- CMAES_control[[original_group_index]]
          CMAES_control_tmp[[original_group_index+1]] <- CMAES_control[[original_group_index]]

          CMAES_control_tmp[[original_group_index]]$mu <- max(c(4,round(CMAES_control_tmp[[original_group_index]]$mu/2)))
          CMAES_control_tmp[[original_group_index+1]]$mu <- max(c(4,round(CMAES_control_tmp[[original_group_index+1]]$mu/2)))
          CMAES_control_tmp[[original_group_index]]$lambda <- max(c(4,round(CMAES_control_tmp[[original_group_index]]$lambda/2)))
          CMAES_control_tmp[[original_group_index+1]]$lambda <- max(c(4,round(CMAES_control_tmp[[original_group_index+1]]$lambda/2)))

          # expand group portion
          groupPortion_tmp <- groupPortion
          groupPortion_tmp[original_group_index] <- groupPortion[original_group_index]
          groupPortion_tmp[original_group_index+1] <- groupPortion[original_group_index]

          # separate the covariance matrix, ps, and pc if exists
          if(!is.null(CMAES_control_tmp[[original_group_index]]$cov)){
            first_cov_size <- length(first_group_indices)

            tmp_cov_matrix <- matrix(,nrow=round(original_group_size/2),ncol=round(original_group_size/2))
            for(ii in 1:first_cov_size){
              for(jj in 1:first_cov_size){
                tmp_cov_matrix[ii,jj] <- CMAES_control[[original_group_index]]$cov[first_group_indices[ii],first_group_indices[jj]]
              }
            }
            CMAES_control_tmp[[original_group_index]]$cov <- tmp_cov_matrix

            # start_second_cov <- round(original_group_size/2)+1
            second_cov_size <- length(second_group_indices)
            tmp_cov_matrix <- matrix(,nrow=second_cov_size,ncol=second_cov_size)
            for(ii in 1:second_cov_size){
              for(jj in 1:second_cov_size){
                tmp_cov_matrix[ii,jj] <- CMAES_control[[original_group_index]]$cov[second_group_indices[ii],second_group_indices[jj]]
              }
            }
            CMAES_control_tmp[[original_group_index+1]]$cov <- tmp_cov_matrix

            # separate pc
            CMAES_control_tmp[[original_group_index]]$pc <- CMAES_control[[groupIndex]]$pc[first_group_indices]
            CMAES_control_tmp[[original_group_index+1]]$pc <- CMAES_control[[groupIndex]]$pc[second_group_indices]

            # separate ps
            CMAES_control_tmp[[original_group_index]]$ps <- CMAES_control[[groupIndex]]$ps[first_group_indices]
            CMAES_control_tmp[[original_group_index+1]]$ps <- CMAES_control[[groupIndex]]$ps[second_group_indices]
          }

          # shift the rest
          if((original_group_index+1)<=original_group_count)
            for(mod_index in (original_group_index+1):original_group_count){
              modified_group_set[[mod_index+1]] <- group[[mod_index]]
              CMAES_control_tmp[[mod_index+1]] <- CMAES_control[[mod_index]]
              groupPortion_tmp[mod_index+1] <- groupPortion[mod_index]
            }


          # replace the grouping and control with the added group
          group <- modified_group_set
          CMAES_control <- CMAES_control_tmp
          groupPortion <- groupPortion_tmp
          # print(groupPortion)
          groupIndex <- groupIndex+1 # a jump to skip the newly created group
          # print(paste0('index: ', groupIndex))
        }
      }else if(any(termination_code==c(2,3,4,5))){
        print('expanding population')
        # IPOP-CMAES -> double pop size
        CMAES_control[[groupIndex]]$mu <- 2*CMAES_control[[groupIndex]]$mu
        CMAES_control[[groupIndex]]$lambda <- 2*CMAES_control[[groupIndex]]$lambda

        #CMAES_control[[groupIndex]]$maxit <- round(CMAES_control[[groupIndex]]$maxit/2)
        # if(any(termination_code==c(3,4))){
        #   CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma * 10
        # }
        # groupIndex <- groupIndex-1
      }
        }


      nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)

      if(nlogging_this_layer>0){
        for(i in 1:nlogging_this_layer){
          nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
          nGeneration_to_consider <- floor(nEval_to_logging/lambda)

          if(!is.matrix(best$diagnostic$value)){
            best$diagnostic$value <- matrix(best$diagnostic$value)
          }
          bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
          convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
        }
        save(convergence_history,file='conv.Rdata')
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
      groupIndex <- groupIndex+1
      print(bestObj)
    }

    successful_termination <- term_code_zero/term_code_total
    maxit_multiplier <- (successful_termination-0.5)/2+1
    message('======Adapting Learning Period========')
    for(groupIndex in 1:nLevel){
      CMAES_control[[groupIndex]]$maxit <- ceiling(CMAES_control[[groupIndex]]$maxit*maxit_multiplier)
    }

    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget,budget,nEval))
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
