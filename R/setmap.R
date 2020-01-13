multilevel_setmap <- function(group,contextVector,mainfun,NP=NULL,
                              ubound,lbound,
                              budget=3000000,eval_interval=100000,
                              fname_prefix="datalog_",control=NULL,
                              surrogate_type=2,bestval=Inf,...){
  nGroup <- length(group)
  globalBest <- Inf
  if(is.null(NP)){
    NP <- rep(20,nGroup)
  }
  nVar <- length(cv)
  groupLength <- pracma::zeros(1,nGroup)
  for(i in 1:nGroup){
    groupLength[i] <- length(group[[i]])
  }

  globalBestPop <- contextVector

  # build and collect surrogates
  # by optimizing lower levels
  # and use the built setmap to optimize higher levels
  setmap <- vector(mode="list",nGroup-1)
  rf_surrogate_all_level <- vector(mode='list',nGroup-1)
  for(levelIndex in nGroup:2){
    print(paste('level',levelIndex))
    parCollection <- NULL
    upperLevel <- levelIndex - 1
    upperGroup <- group[[upperLevel]]
    currentGroup <- group[[levelIndex]]
    currentVector <- contextVector
    pop <- InitializePopulationSobol(groupLength[upperLevel],NP[upperLevel],
                                     minVal = lbound[upperGroup],
                                     maxVal = ubound[upperGroup])
    for(i in 1:NP[upperLevel]){
      currentVector[upperGroup] <- pop[i,]
      if(levelIndex < nGroup){
        setmap <- setmap[(levelIndex):(nGroup-1)]
        furtherGroup <- (levelIndex+1):nGroup
      }else{
        setmap <- NULL
        furtherGroup <- NULL
      }

      funValue <- EvaluateLowerLevel_setmap(mainfun,
                                            currentVector,
                                            group[[levelIndex]],group[furtherGroup],
                                            NP[levelIndex],
                                            lbound=lbound,ubound=ubound,
                                            prevnEval = nEval,
                                            eval_interval,bestval=Inf,
                                            setmap=setmap,...)
      parCollection <- rbind(parCollection,funValue$par)
      if(funValue$value < gb){
        #bestval <- funValue$value
        cv[currentGroup] <<- funValue$par
        gb <<- funValue$value
      }
    }
    # build surrogate/mapping for each lower level variable (parCollection)
    rf_surrogate_this_level <- vector(mode='list',length(currentGroup))
    for(varIndex in 1:length(currentGroup)){
      rf_surrogate_this_level[varIndex] <- list(model=randomForest::randomForest(pop,parCollection[,varIndex],ntree=1000))
    }
    save(pop,parCollection,file = paste0('rfbase',levelIndex,'.Rdata'))
    # compound the surrogates into a list, setmap
    rf_surrogate_all_level[levelIndex-1] <- list(models=rf_surrogate_this_level)
  }
  save(rf_surrogate_all_level,file='setmapTop.Rdata')
  # use the built surrogates in the main optimization loop
  furtherGroup <- group[-1]
  topGroup <- group[[1]]

  budget_left <- budget-nEval
  print('start main optim loop')
  while (budget_left > 0 ) {
    pop <- InitializePopulationSobol(groupLength[1],NP[1],
                                     minVal = lbound[1],
                                     maxVal = ubound[1])
    funVal <- sansde(fname=lowerLevelOptim,
                     pop=pop,
                     furtherGroup=furtherGroup,
                     groupMember=topGroup,
                     bestmem=currentVector[topGroup],
                     bestval=bestval,
                     Lbound=lbound[topGroup],
                     Ubound=ubound[topGroup],
                     contextVector=currentVector,
                     mainfun=fun,setmap=rf_surrogate_all_level,
                     ...)

    if(funval$value < gb){
      contextVector[topGroup] <- funVal$par
      gb <<- funVal$value
      print(paste("bestval update:",gb))
    }

    budget_left <- budget-nEval
  }
  # TODO: verify the models
  # TODO: change number of iterations

  return(list(par=globalBestPop,value=globalBest))
}
