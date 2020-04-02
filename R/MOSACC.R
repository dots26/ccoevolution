#' Main function for the cooperative coevolution with sensitivity-analysis-based grouping
#' for multiobjective optimization. The CMA-ES is used as solver and its
#' step size is set based on  \code{lbound} and  \code{ubound}. This means \code{lbound} and  \code{ubound}
#' must be finite.
#'
#' @title Cooperative coevolution multiobj.
#' @param contextVectors Initial top level population
#' @param objs Initial top level objective values
#' @param fun The objective function object to be solved
#' @param group Vector of list. Each list contains members of the group. If not supplied, the grouping will use the chosen SA_method.
#' @param variable_importance A matrix, each row containing variable importance w.r.t. an objective,
#' i.e. row 1 is variable importance for objective 1.
#' Used if group is supplied. If NULL, all variables will have same portion.
#' @param nLevel Number of level for Morris method. Used if group is not supplied
#' @param lbound Lower bound of the decision vector. A vector with finite values.
#' @param ubound upper bound of the decision vector. A vector with finite values.
#' @param SA_method sensitivity analysis method. It is recommended to use the morris method: \code{'morris_mu'}.
#' Although other SA_method should not cause an error, the package is mostly tested with \code{'morris_mu'}.
#' @param keepCovariance Set whether the covariance matrix is persistent througout the run.
#' If false, the Cov matrix will reset to identity matrix in each cycle.
#' @param evalInterval Interval for data logging.
#' @param inputScale If \code{TRUE}, the input vector will be scaled to [0,1] for processing,
#' otherwise the it will be in the range \code{[lbound,ubound]}
#' @param ... Further arguments passed to \code{fun}.
#' @examples
#' optimum <- rep(13,1000)
#' func <- f1cec
#' ctrl <- list(lbound=rep(-100,1000),ubound=rep(100,1000),delta=rep(20,1000))
#' MOSACC(nVar = 1000,fun=func,lbound=rep(-100,1000),ubound=rep(100,1000),o=optimum)

MOSACC_ori <- function(contextVectors=NULL,objs=NULL,popSize=50,fun,...,
                   group=NULL,variable_importance=NULL,
                   budget=3000000,
                   lbound=rep(0,nVar),ubound=rep(1,nVar),
                   nLevel=4,evalInterval=100000,
                   SA_method=c('morris_mu','morris_k','rf','sobol'),
                   keepCovariance = T,inputScale=T,control=NULL,vectorized=F, repair=F,
                   scaleObjective=T,hist_save=F,hist_filename="hist.Rdata"){

  ############ initialization ###############

  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v))
      return (default)
    else
      return (v)
  }
  if(!is.null(contextVectors)){
    nVar<- nrow(contextVectors)
    nObj <- nrow(objs)
    contextVectorCount <- ncol(contextVectors)
  }else{
    nVar <- controlParam("nVar", 1000)
    nObj <- controlParam("nObj", 2)
    contextVectorCount <- 0
  }
  learning_period <- controlParam("learning_period",500*m)

  SA_method <- SA_method[1]
  nEval <- 0
  convergence_history <- NULL
  nEvalHistory <- NULL
  vars <- 1:nVar
  groupWeight <- NULL
  m <- nObj
  bestHV <- 0
  # contextVectors <- NULL
  #group <- NULL


  ########### grouping ################
  if(is.null(group)){
    if(SA_method == 'morris_mu') {
      ## TODO: fill this from MOCC_main
    }
  }else{
    if(is.null(variable_importance)){
      reg_coefficient <- matrix(rep(1,nVar*m),nrow=m) # set same coefficient for all
    }else{
      reg_coefficient <- variable_importance # copy from argument
    }
    ## determine group weight from variable importance
    groupIndex <- 1
    for(i in 1:nLevel){
      groupMember <- group[[groupIndex]]
      groupWeight <- append(groupWeight,(sum(reg_coefficient[groupMember] )))
      groupIndex <- groupIndex + 1
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

  budget_left <- budget-nEval

  print('init pop')
  while(contextVectorCount < popSize){
    # print(contextVectorCount)
    x <- MaOEA::InitializePopulationLHS(numberOfIndividuals = popSize-contextVectorCount,
                                        chromosomeLength = nVar,
                                        minVal = lbound,maxVal = ubound)
    save(x,file='initx.Rdata')
    objs <- cbind(objs,fun(x,...))
    contextVectors <- cbind(contextVectors,x)
    contextVectorCount <- contextVectorCount + ncol(x)
    nEval <- nEval+popSize-contextVectorCount
  }

  print(group)
  if(inputScale){
    scaling_shift <- lbound
    scaling_factor <- 1/(ubound-lbound)
    ubound <- rep(1,nVar)
    lbound <- rep(0,nVar)
  }else{
    scaling_shift <- rep(0,nVar)
    scaling_factor <- rep(1,nVar)
  }

  CMAES_control <- list()
  nGroup <- length(group)
  for(groupIndex in 1:nGroup){
    groupSize <- length(group[[groupIndex]])
    groupMember <- group[[groupIndex]]
    currentGroupPortion <- groupPortion[groupIndex]
    CMAES_control[[groupIndex]] <- list(vectorized=vectorized,
                                        # mu=groupSize,lambda=groupSize,
                                        mu=floor(4+floor(3*log(groupSize))/2),lambda=4+floor(3*log(groupSize)),
                                        maxit=round(learning_period*(currentGroupPortion/totalPortion)),
                                        sigma=0.1*min(ubound[groupMember]-lbound[groupMember]),
                                        diag.value=T)
  }

  rnk <- nsga2R::fastNonDominatedSorting(t(objs))
  objs <- objs[,rnk[[1]]]
  contextVectors <- contextVectors[,rnk[[1]]]
  contextVectorCount <- length(rnk[[1]])

  print('start optim')
  # updatedFront <- F
  while(budget_left>0){
    shuffle <- fgpt::fyshuffle(1:contextVectorCount)
    print(shuffle)
    for (CVIndex in shuffle) {
      groupIndex <- 1
      while(groupIndex<=length(group)){
        print(paste('group',groupIndex))
        term_code_total <- 0
        term_code_zero <- 0
        updatedFront <- F
        add_point <- F
        if(contextVectorCount < popSize){
          add_point <- T
          CMAES_control[[groupIndex]]$maxit <- round(0.2*learning_period*(currentGroupPortion/totalPortion))
        }else{
          CMAES_control[[groupIndex]]$maxit <- round(learning_period*(currentGroupPortion/totalPortion))
        }

        currentNadir <- vector()
        for(objIndex in 1:nObj)
          currentNadir[objIndex] <- max(objs[objIndex,])

        reference <- currentNadir*1.1
        scale <- 1
        if(scaleObjective){
          scale <- 1/currentNadir
        }
        bestHV <- MaOEA::GetHypervolume(objs,reference)
        contextVector <- contextVectors[,CVIndex]
        groupMember <- group[[groupIndex]]

        best <- cma_es(par=contextVector[groupMember],contextVector = contextVector,groupMember=groupMember,
                       fn = optimizeHV, nondom_front=objs, reference=reference,
                       removed_point=CVIndex,objectiveScale=scale,
                       mainfun=fun,...,add_point=add_point,
                       lower =  lbound[groupMember], upper=ubound[groupMember],
                       control = CMAES_control[[groupIndex]],repair = repair,
                       inputScaleFactor=scaling_factor[groupMember], inputScaleShift=scaling_shift[groupMember])

        termination_code <- best$termination_code
        print(paste('term',termination_code))

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

        if(any(termination_code==c(2,3,4,5))){
          print('expanding population')
          # IPOP-CMAES -> double pop size
          CMAES_control[[groupIndex]]$mu <- 2*CMAES_control[[groupIndex]]$mu
          CMAES_control[[groupIndex]]$lambda <- 2*CMAES_control[[groupIndex]]$lambda
        }
        nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)

        nEval <- nEval + best$counts[1]
        budget_left <- budget-nEval
        print(budget_left)
        if(budget_left>0){
          if(!is.null(best$value)){
            newPar <- contextVectors[,CVIndex]
            newPar[groupMember] <- best$par

            newCV <- contextVectors
            newObj <- fun(newPar,...)
            newObj_set <- objs

            if(budget_left<=0){ break}

            if(!add_point){
              newObj_set[,CVIndex] <- newObj
              newCV[,CVIndex] <- newPar
              rnk <- nsga2R::fastNonDominatedSorting(t(newObj_set))

              newHV <- MaOEA::GetHypervolume(newObj_set[,rnk[[1]]])
            }else{
              print('adding')
              newObj_set <- cbind(newObj_set,newObj)
              newCV <- cbind(newCV,newPar)
              rnk <- nsga2R::fastNonDominatedSorting(t(newObj_set))

              # graphics.off()
              # x <- (0:1000)/1000
              # y <- (1-x^2)^0.5
              # plot(x,y,xlim=c(0,2),ylim=c(0,2))
              # points(objs[1,]/2,objs[2,]/4,type='p')
              # points(newObj/c(2,4),col='green')

              newHV <- MaOEA::GetHypervolume(newObj_set[,rnk[[1]]])
            }

            nEval <- nEval + 1

            message(paste('new HV:',newHV,bestHV,add_point))
            if(bestHV< newHV){
              # update nondominated front
              if(contextVectorCount==length(rnk[[1]])){
                if(add_point){
                  print('new point dominates')
                  removedPoint <- rnk[[2]][1]
                  print(paste('removed',removedPoint,'while improving',CVIndex,'to add to',contextVectorCount))
                  if(removedPoint<=contextVectorCount){
                    contextVectors[,removedPoint] <- newPar
                    objs[,removedPoint] <- newObj
                  }
                }else{
                  contextVectors <- newCV
                  objs <- newObj_set
                }
              }
              if(contextVectorCount!=length(rnk[[1]])){
                print('front changed!')
                contextVectors <- newCV[,rnk[[1]]]
                objs <- newObj_set[,rnk[[1]]]
                contextVectorCount <- length(rnk[[1]])
                updatedFront<- T
              }

              print(paste('Updating context vector number:',CVIndex))
              print(apply(X=objs/c(2,4),FUN = norm,type='2',MARGIN = 2))
              print(paste('Relative to target:',MaOEA::GetHypervolume(objs,reference=1.1*c(2,4))))
              convergence_history <- append(convergence_history,MaOEA::GetHypervolume(objs,reference=1.1*c(2,4)))
              nEvalHistory <- append(nEvalHistory,nEval)
              hist <- list(x=contextVectors,y=objs,convergence_history=convergence_history,nEvalHistory=nEvalHistory)
              if(hist_save)
                save(hist,file=hist_filename)
            }
          }
        }
        groupIndex <- groupIndex+1  ## move this up if groups are nested

      }
      if(budget_left<=0){ break}
      print(paste('next group:',groupIndex))
    }
  }

  return(list(x=contextVectors,y=objs,convergence_history=convergence_history,nEvalHistory=nEvalHistory))
}

# TODO:go through only non-dom. use SMS-EMOA rule to remove points

optimizeHV <- function(x, contextVector,groupMember,mainfun,...,
                       nondom_front,
                       removed_point,
                       add_point=F,
                       reference=NULL,
                       objectiveScale=rep(1,nrow(nondom_front))){
  complete_x <- matrix(contextVector)
  if(!is.null(reference))
    reference <- reference*objectiveScale
  if(is.vector(x)) {
    popSize <- 1
    x <- matrix(x)
  }else{
    popSize <- ncol(x)
  }
  nondom_front <- nondom_front*objectiveScale
  # print(nondom_front)
  par <- x
  complete_x <- pracma::repmat(complete_x,n=1,m=popSize)
  complete_x[groupMember,] <- x
  newObj <- mainfun(complete_x,...)*objectiveScale

  # original HV when all objs are used
  if(!add_point)
    hv_pre <- MaOEA::GetHypervolume(nondom_front[,-removed_point],reference)
  if(add_point)
    hv_pre <- MaOEA::GetHypervolume(nondom_front,reference)

  graphics.off()
  x <- (0:1000)/1000
  y <- (1-x^2)^0.5
  plot(x,y,xlim=c(0,2),ylim=c(0,2))
  points(nondom_front[1,]/2/objectiveScale[1],nondom_front[2,]/4/objectiveScale[2],col='red')
  points(newObj[1,]/2/objectiveScale[1],newObj[2,]/4/objectiveScale[2],col='blue')

  objectiveValue <- NULL

  for(individualIndex in 1:popSize){
    nondom_front_tmp <- nondom_front
    is.dominated_by_ref <- all(newObj[,individualIndex]>reference)
    if(is.dominated_by_ref){
      new_ref <- vector()
      for(objIndex in 1:nrow(nondom_front)){
        new_ref[objIndex] <- max(newObj[objIndex,individualIndex],reference[objIndex])
      }

      hv_new <- MaOEA::GetHypervolume(nondom_front_tmp,new_ref)
      objectiveValue <- append(objectiveValue,-hv_new)
    }else{
      nondom_front_tmp <- cbind(nondom_front,newObj[,individualIndex])
      if(!add_point)
        nondom_front_tmp <- nondom_front_tmp[,-removed_point,drop=F]

      nPoint <- ncol(nondom_front_tmp)
      rnk <- nsga2R::fastNonDominatedSorting(t(nondom_front_tmp))

      rnkIndex <- integer(ncol(nondom_front_tmp))
      sortingIndex <- 1;
      while (sortingIndex <= length(rnk)) {
        rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
        sortingIndex <- sortingIndex + 1;
      }
      newPoint_rnk <- rnkIndex[nPoint]

      if(newPoint_rnk==1){
        partly_exceed_ref <- any(newObj[,individualIndex]>reference)
        if(partly_exceed_ref){
          new_ref <- vector()
          for(objIndex in 1:nrow(nondom_front_tmp)){
            new_ref[objIndex] <- max(newObj[objIndex,individualIndex]*1.1,reference[objIndex])
          }
          hv_new <- MaOEA::GetHypervolume(nondom_front_tmp[,rnk[[1]]],new_ref)# diversity
        }else{
          hv_new <- MaOEA::GetHypervolume(nondom_front_tmp[,rnk[[1]]],reference) # convergence
        }
      }else{
        hv_new <- MaOEA::GetHypervolume(newObj[,individualIndex],reference)
      }
      objectiveValue <- append(objectiveValue,hv_new-hv_pre)
    }
  }

  # points(newObj[1,which.min(-objectiveValue)]/2/objectiveScale[1],newObj[2,which.min(-objectiveValue)]/4/objectiveScale[2],col='blue')
  # print(min(-objectiveValue))

  return(-objectiveValue)
}

