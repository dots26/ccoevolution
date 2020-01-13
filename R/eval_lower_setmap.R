EvaluateLowerLevel_setmap <- function(fun,currentVector,
                                      nextGroup,furtherGroup=NULL,
                                      NP_next,
                                      controlNext,controlFurther,
                                      lbound,ubound,
                                      prevnEval=0,eval_interval=100000,
                                      maxIter=5,fname_prefix="datalog_",bestval=Inf,
                                      setmap=NULL,lowersetmap=NULL,...){
  optimPar <- currentVector[nextGroup]
  optimVal <- bestval

  if(is.null(setmap)){ # if no mapping available, i.e. lowest level
    # generate lower level population
    pop <- InitializePopulationLHS(ncol=length(nextGroup),
                                   nrow=NP_next,
                                   minVal = lbound[nextGroup],
                                   maxVal = ubound[nextGroup])
    # evaluate and optimize the real function
    funVal <- cma_es(fn = subfunctionCMA,
                     par=optimPar,
                     groupMember = nextGroup,
                     lower=lbound[nextGroup],
                     upper=ubound[nextGroup],
                     contextVector=currentVector,
                     mainfun=fun,control=list(maxit=2000,
                                              mu=2,lambda=10,
                                              sigma=0.3*max(ubound[nextGroup]-lbound[nextGroup])),
                     ...)
    # funVal <- sansde(fname=subfunction,
    #                  pop=pop,
    #                  groupMember = nextGroup,
    #                  bestmem=currentVector[nextGroup],
    #                  bestval=bestval,
    #                  Lbound=lbound[nextGroup],
    #                  Ubound=ubound[nextGroup],
    #                  contextVector=currentVector,
    #                  mainfun=fun,control=list(itermax=200),
    #                  ...)
    prevnEval <- nEval
    nEval <<- nEval + funVal$counts[1]

    if(funVal$value<bestval){
      optimVal <- funVal$value
      optimPar <- funVal$par
      bestval <- optimVal
    }

    previousLogged <- prevnEval%/%eval_interval
    print(paste('neval',nEval,previousLogged))
    logData <- ((nEval)%/%eval_interval)-previousLogged
    if(logData){
      argList <- list(...)
      save(optimVal,optimPar,argList,file=paste0(fname_prefix,previousLogged+1,'.Rdata'))
    }
  }else{ # when mapping available, i.e. not lowest
    # optimize, but
    # use setmap successively to obtain the best pop for the lower levels
    # evaluate and optimize the real function
    # generate lower level population
    pop <- InitializePopulationLHS(ncol=length(nextGroup),
                                   nrow=NP_next,
                                   minVal = lbound[nextGroup],
                                   maxVal = ubound[nextGroup])

    funVal <- cma_es(fn = lowerLevelOptim,
                     par=currentVector[nextGroup],
                     groupMember = nextGroup,
                     lower=lbound[nextGroup],
                     upper=ubound[nextGroup],
                     contextVector=currentVector,
                     mainfun=fun,control=list(maxit=100,
                                              mu=NP_next,lambda=NP_next,
                                              sigma=0.3*max(ubound[nextGroup]-lbound[nextGroup])),
                     CMA=T,setmap=setmap,...)


    # funVal <- sansde(fname=lowerLevelOptim,
    #                  pop=pop,
    #                  furtherGroup=furtherGroup,
    #                  groupMember=nextGroup,
    #                  bestmem=currentVector[nextGroup],
    #                  bestval=bestval,
    #                  Lbound=lbound[nextGroup],
    #                  Ubound=ubound[nextGroup],
    #                  contextVector=currentVector,
    #                  mainfun=fun,setmap=setmap,control=list(itermax=200),
    #                  ...)
    if(funVal$value<bestval){
      optimVal <- funVal$value
      optimPar <- funVal$par
      bestval <- optimVal
    }
    prevnEval <- nEval
  }

  return(list(par=optimPar,value=optimVal))
}

lowerLevelOptim <- function(population,contextVector,groupMember,furtherGroup=NULL,mainfun,setmap,CMA=F,...){
  if(is.matrix(population)){
    if(CMA)
      population <- t(population)
    popSize <- nrow(population)
  }
  if(is.vector(population))
    popSize <- 1

  nVar <- length(contextVector)
  contextVector <- rep(contextVector,popSize)

  contextVector <- t(matrix(contextVector,nrow=nVar))
  contextVector[,groupMember] <- population
  # setmap length = nlowerlevel groups
  # same as further group
  # each setmap map an upperlevel to furthergroup[]
  upperPop <- population
  if(length(furtherGroup)>0){
    for(lowerIndex in 1:length(furtherGroup)){
      nextGroup <- furtherGroup[[lowerIndex]]
      lowerBestPop <-  NULL
      for(topIndex in 1:nrow(population)){
        save(setmap,file='setmapLower.Rdata')
        lowerBestPop <- rbind(lowerBestPop,predictPop(setmap[lowerIndex],upperPop[topIndex,]))
        contextVector[topIndex,nextGroup] <- lowerBestPop[topIndex,]
      }
      upperPop <- lowerBestPop
    }
  }

  objectiveValue <- mainfun(contextVector,...)
}

predictPop <- function(setmap,individual){
  setmap_length <- length(setmap$models) # number of variables to be predicted
  map_result <- vector(mode ="numeric",setmap_length)
  for(i in 1:setmap_length){
    map_result[i] <-  predict(setmap$models[i]$model,individual)$predictions
  }
  return(map_result)
}
