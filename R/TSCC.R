#' Two stage cooperative coevolution with CMA-ES
#' The CMA-ES is used as solver and its
#' step size is set based on  \code{lbound} and  \code{ubound}.
#' This means \code{lbound} and  \code{ubound}
#' must be finite.
#'
#' @title Cooperative coevolution with two stage grouping.
#' @param population Initial population
#' @param fun The objective function object to be solved
#' @param lbound Lower bound of the decision vector. A vector with finite values.
#' @param ubound upper bound of the decision vector. A vector with finite values.
#' @param first_grouping Vector of list. Each list should contain the group members.
#' If not supplied, the grouping will use the chosen SA_method.
#' This group would be subject to the method: \code{DG2}.
#' @param variable_importance Used if first_grouping is supplied. If NULL, all variables will have same importance.
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
#' TSCC(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
TSCC <- function(contextVector=NULL,nVar,
                 fun,...,
                 first_grouping=NULL, variable_importance=NULL,
                 budget=3000000,lbound=rep(0,nVar),ubound=rep(1,nVar),
                 nLevel=4,evalInterval=100000,
                 SA_method=c('morris_mu','morris_k','rf','sobol'),
                 keepCovariance=F,scale=T,disableIPOP=F){

  ########### initialization ###############
  recordConv <- NULL
  recordNEval <- NULL
  if(!is.null(contextVector))
    nVar <- length(contextVector)
  SA_method <- SA_method[1]
  groupWeight <- NULL
  nEval <- 0
  convergence_history <- NULL
  group <- NULL
  vars <- 1:nVar
  learning_period <- 500
  clusterOrder <- 1:nLevel

  ########### grouping ###############
  group <- first_grouping
  if(is.null(group)){ ## TODO: what should be done if group is supplied
    if(SA_method == 'morris_mu' || SA_method == 'morris_k') {
      r<- 20
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
          #orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
          group <- append(group,list(unorderedCurrentGroup))
          groupWeight <- append(groupWeight,sum(reg_coefficient[unorderedCurrentGroup]))
        }
        start <- (nLevel-1)*nLevelVar+1
        end <- nVar
        unorderedCurrentGroup <- rankedVars[start:end]
        #orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
        groupWeight <- append(groupWeight,sum(reg_coefficient[unorderedCurrentGroup]))
        group <- append(group,list(unorderedCurrentGroup))
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
  print('ordered group')
  print(group)
  print(c(groupPortion,totalPortion))
  print('Secondary Grouping...')
  new_group <- NULL
  dg <- NULL

  if(is.null(contextVector)){
    contextVector <- matrix(runif(nVar)*(ubound-lbound)+lbound,nrow=1)
    bestObj <- fun(contextVector,...)
    nEval <- nEval + 1
  }

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
    copy_sep <- subgroup$separable
    subgroup$separable <- list()
    subgroup$separable[[1]] <- copy_sep
    dg <- cbind(dg,subgroup)
    print('groups')
    print(subgroup$group)
    print('separable')
    print(subgroup$separable)
    nEval <- nEval + subgroup$nEval
  }
  dg<- dg[,clusterOrder,drop=F]
  leftBudget <- budget - nEval

  ########### initialization of CMAES control param ###############
  CMAES_control <- vector(mode = "list",length=nLevel)
  for(clusterIndex in 1:nLevel){
    cluster_member <- group[[clusterIndex]]
    currentClusterGrouping <- dg[,clusterIndex]$group
    currentClusterSep <- dg[,clusterIndex]$separable
    currentGroupPortion <- groupPortion[clusterIndex]
    CMAES_control[[clusterIndex]]$sep <- list()
    CMAES_control[[clusterIndex]]$nonsep <- list()

    if(length(currentClusterSep)>0){
      groupSize <- length(currentClusterSep[[1]])
      groupMember <- currentClusterSep[[1]]
      CMAES_control[[clusterIndex]]$sep[[1]] <- list(vectorized=T,
                                                     # mu=groupSize,lambda=groupSize,
                                                     mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                                                     maxit=round(learning_period*(currentGroupPortion/totalPortion)),
                                                     sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                                     diag.value=T)

    }
    if(length(currentClusterGrouping)>0){
      for(groupIndex in 1:length(currentClusterGrouping)) {
        groupMember <- cluster_member[currentClusterGrouping[[groupIndex]]]
        groupSize <- length(groupMember)
        CMAES_control[[clusterIndex]]$nonsep[[groupIndex]] <- list(vectorized=T,
                                                                   # mu=groupSize,lambda=groupSize,
                                                                   mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                                                                   maxit=round(learning_period*(currentGroupPortion/totalPortion)),
                                                                   sigma=0.3*max(ubound[groupMember]-lbound[groupMember]),
                                                                   diag.value=T)
      }
    }
  }

  nRepeat <- floor(nEval/evalInterval)
  if(nRepeat>0){
    for(i in 1:nRepeat){
      convergence_history <- append(convergence_history,bestObj)
    }
  }
  print(nEval)
  print(bestObj)
  ########### end initialization& grouping  ###############

  if(scale){
    scaling_shift <- lbound
    scaling_factor <- 1/(ubound-lbound)
    ubound <- rep(1,nVar)
    lbound <- rep(0,nVar)
  }else{
    scaling_shift <- rep(0,nVar)
    scaling_factor <- rep(1,nVar)
  }

  ########### begin cooperative coevolution ###############
  while((budget-nEval)>0 ){
    term_code_total <- 0
    term_code_zero <- 0
    for(clusterIndex in 1:nLevel){
      print(paste('Optimizing Variable level',clusterIndex))
      cluster_member <- group[[clusterIndex]]
      currentClusterGrouping <- dg[,clusterIndex]$group
      currentClusterSep <- dg[,clusterIndex]$separable

      currentGroupPortion <- groupPortion[clusterIndex]

      # optimize separable
      if(length(currentClusterSep)>0){
        print('OPT SEP')
        groupIndex <- 1
        while(groupIndex <= length(currentClusterSep)) {
          groupMember <- cluster_member[currentClusterSep[[groupIndex]]]
          groupSize <- length(groupMember)

          best<- sep_cma_es(contextVector[groupMember],
                            fn = subfunctionCMA,
                            contextVector = contextVector,
                            groupMember = groupMember,
                            mainfun=fun,...,
                            lower = lbound[groupMember],
                            upper=ubound[groupMember],
                            # control = list(mu=groupSize,lambda=groupSize,maxit=2000))
                            control = CMAES_control[[clusterIndex]]$sep[[groupIndex]],
                            inputScaleFactor=scaling_factor[groupMember],
                            inputScaleShift=scaling_shift[groupMember])
          termination_code <- best$termination_code
          term_code_total <- term_code_total + 1
          # termination codes:
          # 0: no error -> do nothing
          # 1: sigma divergence -> create smaller groups
          # 2: no improvement for 20 generations -> increase pop size
          # 3: sigma is too small to have effect on covariance matrix -> increase pop size
          # 4: sigma is too small to have effect on mean x -> increase pop size
          # 5: condition number too large -> increase pop size

          if(termination_code==0){
            term_code_zero <- term_code_zero + 1
            if(keepCovariance){
              CMAES_control[[clusterIndex]]$sep[[groupIndex]]$cov <- best$cov
              CMAES_control[[clusterIndex]]$sep[[groupIndex]]$sigma <- best$sigma
              CMAES_control[[clusterIndex]]$sep[[groupIndex]]$ps <- best$ps
              CMAES_control[[clusterIndex]]$sep[[groupIndex]]$pc <- best$pc
            }
          }
          lambda <- CMAES_control[[clusterIndex]]$sep[[groupIndex]]$lambda
          if(termination_code==1){
            message("separating group")
            # create smaller groups!
            # modify dg[,clusterIndex]$group
            original_group_index <- groupIndex
            original_group <- dg[,clusterIndex]$separable[[original_group_index]]
            original_group_size <- length(original_group) # number of variable in the group
            original_group_count <- length(dg[,clusterIndex]$separable) # number of group in the cluster

            modified_group_set <- dg[,clusterIndex]$separable # make a copy to do a shift for the groups
            group_coef <- reg_coefficient[groupMember]
            group_coef_order <- order(group_coef,decreasing = T)
            first_group_indices <- group_coef_order[1:round(original_group_size/2)]
            second_group_indices <- group_coef_order[(round(original_group_size/2)+1):original_group_size]

            modified_group_set[[original_group_index]] <- dg[,clusterIndex]$separable[[original_group_index]][1:round(original_group_size/2)]
            modified_group_set[[original_group_index+1]] <- dg[,clusterIndex]$separable[[original_group_index]][(round(original_group_size/2)+1):original_group_size]

            # create CMAES control for the new groups
            # copy the control for current cluster
            CMAES_control_tmp <- CMAES_control[[clusterIndex]]$sep
            CMAES_control_tmp[[original_group_index]] <- CMAES_control[[clusterIndex]]$sep[[original_group_index]]
            CMAES_control_tmp[[original_group_index+1]] <- CMAES_control[[clusterIndex]]$sep[[original_group_index]]

            CMAES_control_tmp[[original_group_index]]$mu <- max(c(4,round(CMAES_control[[clusterIndex]]$sep[[original_group_index]]$mu/2)))
            CMAES_control_tmp[[original_group_index+1]]$mu <- max(c(4,round( CMAES_control[[clusterIndex]]$sep[[original_group_index]]$mu/2 )))
            CMAES_control_tmp[[original_group_index]]$lambda <- max(c(4,round(CMAES_control[[clusterIndex]]$sep[[original_group_index]]$lambda/2)))
            CMAES_control_tmp[[original_group_index+1]]$lambda <- max(c(4,round(CMAES_control[[clusterIndex]]$sep[[original_group_index]]$lambda/2)))

            # separate the covariance matrix, ps, and pc if exists
            if(!is.null(CMAES_control_tmp[[original_group_index]]$cov)){
              first_cov_size <- length(first_group_indices)

              tmp_cov_matrix <- matrix(,nrow=round(original_group_size/2),ncol=round(original_group_size/2))
              for(ii in 1:first_cov_size){
                for(jj in 1:first_cov_size){
                  tmp_cov_matrix[ii,jj] <- CMAES_control[[clusterIndex]]$sep[[original_group_index]]$cov[first_group_indices[ii],first_group_indices[jj]]
                }
              }
              CMAES_control_tmp[[original_group_index]]$cov <- tmp_cov_matrix

              # start_second_cov <- round(original_group_size/2)+1
              second_cov_size <- length(second_group_indices)
              tmp_cov_matrix <- matrix(,nrow=second_cov_size,ncol=second_cov_size)
              for(ii in 1:second_cov_size){
                for(jj in 1:second_cov_size){
                  tmp_cov_matrix[ii,jj] <- CMAES_control[[clusterIndex]]$sep[[original_group_index]]$cov[second_group_indices[ii],second_group_indices[jj]]
                }
              }
              CMAES_control_tmp[[original_group_index+1]]$cov <- tmp_cov_matrix

              # separate pc
              CMAES_control_tmp[[original_group_index]]$pc <- CMAES_control[[clusterIndex]]$sep[[groupIndex]]$pc[first_group_indices]
              CMAES_control_tmp[[original_group_index+1]]$pc <- CMAES_control[[clusterIndex]]$sep[[groupIndex]]$pc[second_group_indices]

              # separate ps
              CMAES_control_tmp[[original_group_index]]$ps <- CMAES_control[[clusterIndex]]$sep[[groupIndex]]$ps[first_group_indices]
              CMAES_control_tmp[[original_group_index+1]]$ps <- CMAES_control[[clusterIndex]]$sep[[groupIndex]]$ps[second_group_indices]
            }

            # shift the rest
            if((original_group_index+1)<=original_group_count)
              for(mod_index in (original_group_index+1):original_group_count){
                modified_group_set[[mod_index+1]] <- dg[,clusterIndex]$separable[[mod_index]]
                CMAES_control_tmp[[mod_index+1]] <- CMAES_control[[clusterIndex]]$sep[[mod_index]]
              }

            dg[,clusterIndex]$separable <- modified_group_set
            CMAES_control[[clusterIndex]]$sep <- CMAES_control_tmp
            groupIndex <- groupIndex+1 # a jump to skip the newly created group
            currentClusterSep <- dg[,clusterIndex]$separable # update the grouping for local function
          }else if(any(termination_code==c(2,3,4,5))){
            message('expanding population')
            # IPOP-CMAES -> double pop size
            CMAES_control[[clusterIndex]]$sep[[groupIndex]]$mu <- 2*CMAES_control[[clusterIndex]]$sep[[groupIndex]]$mu
            CMAES_control[[clusterIndex]]$sep[[groupIndex]]$lambda <- 2*CMAES_control[[clusterIndex]]$sep[[groupIndex]]$lambda
            if(any(termination_code==c(3,4))){
              #     CMAES_control[[clusterIndex]]$sep[[groupIndex]]$sigma <- CMAES_control[[clusterIndex]]$sep[[groupIndex]]$sigma *10
            }
          }

          # nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
          # if(nlogging_this_layer>0){
          #   for(i in 1:nlogging_this_layer){
          #     nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
          #     nGeneration_to_consider <- floor(nEval_to_logging/lambda)
          #     # print(c('record sep',nGeneration_to_consider,nlogging_this_layer,mu,dim(best$diagnostic$value)))
          #     #print(best$diagnostic)
          #     if(!is.matrix(best$diagnostic$value)){
          #       best$diagnostic$value <- matrix(best$diagnostic$value)
          #     }
          #
          #     bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
          #     convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
          #
          #   }
          # }
          nEval <- nEval + best$counts[1]


          if((budget-nEval)>0){ # only update if it doesnt exceed budget
            if(!is.null(best$par)){
              contextVector[groupMember] <- best$par
              obj <- best$value
              if(obj < bestObj){
                bestPop <- contextVector
                bestObj <- obj

                recordConv <- append(recordConv,bestObj)
                recordNEval <- append(recordNEval,nEval)
              }
            }else{
              #      # print('is null')
            }
          }else{
            break
          }
          groupIndex <- groupIndex+1
          print(bestObj)
        }
      }

      # optimize non-separable

      if(length(currentClusterGrouping)>0){
        print('OPT NONSEP')
        groupIndex <- 1
        while(groupIndex <= length(currentClusterGrouping)) {
          groupMember <- cluster_member[currentClusterGrouping[[groupIndex]]]
          groupSize <- length(groupMember)

          # message(paste0('optimizing group ',groupIndex,' with ',groupSize,' members, ',currentGroupPortion/totalPortion*100,'% portion'))

          best<- cma_es(contextVector[groupMember],
                        fn = subfunctionCMA,
                        contextVector = contextVector,
                        groupMember = groupMember,
                        mainfun=fun,...,
                        lower = lbound[groupMember],
                        upper=ubound[groupMember],
                        control = CMAES_control[[clusterIndex]]$nonsep[[groupIndex]],
                        inputScaleFactor=scaling_factor[groupMember],
                        inputScaleShift=scaling_shift[groupMember],
                        disableIPOP=disableIPOP)
          termination_code <- best$termination_code
          term_code_total <- term_code_total + 1

          # termination codes:
          # 0: no error -> do nothing
          # 1: sigma divergence -> create smaller groups
          # 2: no improvement for 20 generations -> increase pop size
          # 3: sigma is too small to have effect on covariance matrix -> increase pop size
          # 4: sigma is too small to have effect on mean x -> increase pop size
          # 5: condition number too large -> increase pop size

          if(termination_code==0){
            term_code_zero <- term_code_zero + 1
            if(keepCovariance){
              CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$cov <- best$cov
              CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma <- best$sigma
              CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$ps <- best$ps
              CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$pc <- best$pc
            }
          }
          lambda <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$lambda
          if(termination_code==1){
            message("separating group")
            # if(groupSize>1){
            # create smaller groups!
            # modify dg[,clusterIndex]$group
            original_group_index <- groupIndex
            original_group <- dg[,clusterIndex]$group[[original_group_index]]
            original_group_size <- length(original_group) # number of variable in the group
            original_group_count <- length(dg[,clusterIndex]$group) # number of group in the cluster

            modified_group_set <- dg[,clusterIndex]$group # make a copy to do a shift for the groups

            group_coef <- reg_coefficient[groupMember]
            group_coef_order <- order(group_coef,decreasing = T)
            first_group_indices <- group_coef_order[1:round(original_group_size/2)]
            second_group_indices <- group_coef_order[(round(original_group_size/2)+1):original_group_size]

            modified_group_set[[original_group_index]] <- dg[,clusterIndex]$group[[original_group_index]][1:round(original_group_size/2)]
            modified_group_set[[original_group_index+1]] <- dg[,clusterIndex]$group[[original_group_index]][(round(original_group_size/2)+1):original_group_size]

            # create CMAES control for the new groups
            # copy the control for current cluster
            CMAES_control_tmp <- CMAES_control[[clusterIndex]]$nonsep
            CMAES_control_tmp[[original_group_index]] <- CMAES_control[[clusterIndex]]$nonsep[[original_group_index]]
            CMAES_control_tmp[[original_group_index+1]] <- CMAES_control[[clusterIndex]]$nonsep[[original_group_index]]

            CMAES_control_tmp[[original_group_index]]$mu <- max(c(4,round( CMAES_control_tmp[[original_group_index]]$mu /2 )))
            CMAES_control_tmp[[original_group_index+1]]$mu <- max(c(4,round( CMAES_control_tmp[[original_group_index+1]]$mu/2 )))
            CMAES_control_tmp[[original_group_index]]$lambda <- max(c(4,round(CMAES_control_tmp[[original_group_index]]$lambda/2)))
            CMAES_control_tmp[[original_group_index+1]]$lambda <- max(c(4,round(CMAES_control_tmp[[original_group_index+1]]$lambda/2)))

            # separate the covariance matrix, ps, and pc if exists
            if(!is.null(CMAES_control_tmp[[original_group_index]]$cov)){
              first_cov_size <- length(first_group_indices)

              tmp_cov_matrix <- matrix(,nrow=round(original_group_size/2),ncol=round(original_group_size/2))
              for(ii in 1:first_cov_size){
                for(jj in 1:first_cov_size){
                  tmp_cov_matrix[ii,jj] <- CMAES_control[[clusterIndex]]$nonsep[[original_group_index]]$cov[first_group_indices[ii],first_group_indices[jj]]
                }
              }
              CMAES_control_tmp[[original_group_index]]$cov <- tmp_cov_matrix

              # start_second_cov <- round(original_group_size/2)+1
              second_cov_size <- length(second_group_indices)
              tmp_cov_matrix <- matrix(,nrow=second_cov_size,ncol=second_cov_size)
              for(ii in 1:second_cov_size){
                for(jj in 1:second_cov_size){
                  tmp_cov_matrix[ii,jj] <- CMAES_control[[clusterIndex]]$nonsep[[original_group_index]]$cov[second_group_indices[ii],second_group_indices[jj]]
                }
              }
              CMAES_control_tmp[[original_group_index+1]]$cov <- tmp_cov_matrix

              # separate pc
              CMAES_control_tmp[[original_group_index]]$pc <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$pc[first_group_indices]
              CMAES_control_tmp[[original_group_index+1]]$pc <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$pc[second_group_indices]

              # separate ps
              CMAES_control_tmp[[original_group_index]]$ps <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$ps[first_group_indices]
              CMAES_control_tmp[[original_group_index+1]]$ps <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$ps[second_group_indices]
            }

            # shift the rest
            if((original_group_index+1)<=original_group_count){
              for(mod_index in (original_group_index+1):original_group_count){
                modified_group_set[[mod_index+1]] <- dg[,clusterIndex]$group[[mod_index]]
                CMAES_control_tmp[[mod_index+1]] <- CMAES_control[[clusterIndex]]$nonsep[[mod_index]]
              }
            }

            # replace the grouping and control with the added group
            dg[,clusterIndex]$group <- modified_group_set
            CMAES_control[[clusterIndex]]$nonsep <- CMAES_control_tmp
            groupIndex <- groupIndex+1 # a jump to skip the newly created group
            currentClusterGrouping <- dg[,clusterIndex]$group # update the grouping for local function
            # }
          }else if(any(termination_code==c(2,3,4,5))){
            message('expanding population')

            # IPOP-CMAES -> double pop size
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$mu <- 2*CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$mu
            CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$lambda <- 2*CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$lambda
            if(any(termination_code==c(3,4))){
              #CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma <- CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$sigma *10
            }
          }

          nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
          if(nlogging_this_layer>0){
            for(i in 1:nlogging_this_layer){
              # print(c('pre-rec',nEval,best$counts[1]))
              nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
              nGeneration_to_consider <- floor(nEval_to_logging/lambda)
              # print(c('record ns',nGeneration_to_consider,nlogging_this_layer,lambda,dim(best$diagnostic$value)))
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
                print('UPDATE!')
                bestPop <- contextVector
                bestObj <- obj
                recordConv <- append(recordConv,bestObj)
                recordNEval <- append(recordNEval,nEval)
              }
            }else{

            }
          }else{
            break
          }
          groupIndex <- groupIndex+1
          print(bestObj)
        }
      }
      leftBudget <- budget - nEval
      print(c('Comp budget left:',leftBudget,budget,nEval))
    }
    successful_termination <- term_code_zero/term_code_total
    maxit_multiplier <- (successful_termination-0.5)/2+1

    print('======Adapting Learning Period========')
    print(maxit_multiplier)
    for(clusterIndex in 1:nLevel){
      currentClusterGrouping <- dg[,clusterIndex]$group
      currentClusterSep <- dg[,clusterIndex]$separable
      if(length(currentClusterSep)>0)
        for(groupIndex in 1:length(currentClusterSep)){
          CMAES_control[[clusterIndex]]$sep[[groupIndex]]$maxit <- ceiling(CMAES_control[[clusterIndex]]$sep[[groupIndex]]$maxit*maxit_multiplier)
        }
      if(length(currentClusterGrouping)>0)
        for(groupIndex in 1:length(currentClusterGrouping)){
          CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$maxit <- ceiling(CMAES_control[[clusterIndex]]$nonsep[[groupIndex]]$maxit*maxit_multiplier)
        }
    }
  }
  return(list(x=bestPop,y=bestObj,record=list(nEval=recordNEval,conv=recordConv)))
}
