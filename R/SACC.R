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
#' SACC(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)
#' @export
SACC <- function(contextVector=NULL,nVar,fun,...,
                 group=NULL,
                 budget=3000000,
                 lbound=rep(-Inf,nVar),ubound=rep(Inf,nVar),
                 nLevel=4,evalInterval=100000,
                 SA_method=c('morris_mu','morris_k','rf','sobol'),
                 keepCovariance = F,limit_sigma=F){

  ############ initialization ###############
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
                    control = CMAES_control[[groupIndex]])

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
      if(termination_code==1){
        print('separating group')
        # print(paste0('index: ', groupIndex))
        if(groupSize>1){
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
    # if(successful_termination>0.66){ # most subproblems are successful, increase learning period
      #learning_period <- round(learning_period*1.5)
      message('======Extending Learning Period========')
      for(groupIndex in 1:nLevel){
        CMAES_control[[groupIndex]]$maxit <- ceiling(CMAES_control[[groupIndex]]$maxit*maxit_multiplier)
      }
    # }
    # if(successful_termination<0.33){ # most subproblems are failing, reduce learning period
    #   # learning_period <- round(learning_period/1.5)
    #   message('======Reducing Learning Period========')
    #   for(groupIndex in 1:nLevel){
    #     CMAES_control[[groupIndex]]$maxit <- ceiling(CMAES_control[[groupIndex]]$maxit/1.5)
    #   }
    # }
    leftBudget <- budget - nEval
    print(c('Comp budget left:',leftBudget,budget,nEval))
  }
  return(list(x=bestPop,y=bestObj,conv=convergence_history))
}
