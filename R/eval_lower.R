EvaluateLowerLevel <- function(fun,currentVector,
                               nextGroup,furtherGroup,
                               NP_next,NP_further,
                               controlNext,controlFurther,
                               lbound,ubound,
                               prevnEval=0,eval_interval=100000,
                               maxIter=5,fname_prefix="datalog_",bestval=Inf,infill=c("ei","best"),...){
  infill <- infill[1]
  optimPar <- currentVector[nextGroup]
  optimVal <- bestval
  if(length(furtherGroup)==0){
    # pop <- InitializePopulationLHS(ncol=length(nextGroup),
    #                                 nrow=NP_next,
    #                                 minVal = lbound[nextGroup],
    #                                 maxVal = ubound[nextGroup])
    #
    # funVal <- sansde(fname=subfunction,
    #                  pop=pop,
    #                  groupMember=nextGroup,
    #                  bestmem=currentVector[nextGroup],
    #                  bestval=bestval,
    #                  Lbound=lbound[nextGroup],
    #                  Ubound=ubound[nextGroup],
    #                  contextVector=currentVector,
    #                  mainfun=fun,
    #                  ...)
    funVal <- cma_es(fn=subfunctionCMA,
                     contextVector=currentVector,
                     par=currentVector[nextGroup],
                     lower=lbound[nextGroup],
                     upper=ubound[nextGroup],
                     mainfun=fun,
                     groupMember=nextGroup,
                     control=list(maxit=1000,vectorized=T,mu=10,lambda=20,
                                  sigma=0.3*max(ubound[nextGroup]-lbound[nextGroup])))
    print('lower level par')
    print(funVal$par)
    prevnEval <- nEval
    nEval <<- nEval + funVal$counts[1]


    if(funVal$value<optimVal){
      currentVector[nextGroup] <- funVal$par
      optimVal <- funVal$value
      optimPar <- funVal$par
    }
    if(funVal$value<gb){
      gb <<- funVal$value
      cv[nextGroup] <<- optimPar
    }

    previousLogged <- prevnEval%/%eval_interval
    print(paste('neval',nEval,previousLogged))
    logData <- ((nEval)%/%eval_interval)-previousLogged
    if(logData){
      argList <- list(...)
      save(optimVal,optimPar,cv,gb,argList,file=paste0(fname_prefix,previousLogged+1,'.Rdata'))
    }
  }else{
    #generate lower level population
    groupLength <- length(nextGroup)
    pop <- InitializePopulationLHS(groupLength,
                                   NP_next,
                                   minVal = lbound[nextGroup],
                                   maxVal = ubound[nextGroup])

    valCollection <- NULL
    for(i in 1:NP_next){
      print(paste('level',length(furtherGroup),'individual number',i,'of',NP_next))
      currentVector[nextGroup] <- pop[i,]
      funVal <- EvaluateLowerLevel(fun,
                                   currentVector,
                                   furtherGroup[[1]],furtherGroup[-1],
                                   NP_further[1],NP_further[-1],
                                   lbound=lbound,ubound=ubound,
                                   prevnEval = nEval,eval_interval = eval_interval,
                                   fname_prefix=fname_prefix,bestval=bestval,infill=infill,
                                   ...)
      valCollection <- append(valCollection,funVal$value)
    }

    # build surrogate here
    print(paste("surrogate level:",length(furtherGroup)))
    # rf_model <- grf::regression_forest((pop),(valCollection),num.trees=1000)
    # rf_model <- randomForest::randomForest((pop),(valCollection),ntree=1000)
    rf_model <- SPOT::buildKriging((pop),matrix(valCollection))

    bestIndex <- which.min(valCollection)
    bestPar <- pop[bestIndex,]

    if(min(valCollection)<bestval){
      optimVal <- min(valCollection)
      optimPar <- bestPar
      bestval <- optimVal
    }
    if(min(valCollection)<gb){
      gb <<- min(valCollection)
      cv[nextGroup] <<- bestPar
    }
    # sequential optimization
    for(iterOptim in 1:maxIter){
      # predictedBest <- ccoevolution::cma_es(par=bestPar,
      #                                       fn=bestPredModel,
      #                                       model=rf_model,
      #                                       lower = lbound[nextGroup],upper=ubound[nextGroup],
      #                                       control=list(mu=10,lambda=10,maxit=30))$par
      print(paste('level',length(furtherGroup),' surrogate, iter number',iterOptim,'of',maxIter))

      # surrogate_pop <- InitializePopulationLHS(ncol=length(nextGroup),
      #                                          nrow=10,
      #                                          minVal = lbound[nextGroup],
      #                                          maxVal = ubound[nextGroup])
      if(infill=="best"){
        fname=bestPredModel
        # predictedBest <- sansde(fname=fname,
        #                         pop=surrogate_pop,
        #                         model=rf_model,
        #                         bestmem = cv[nextGroup],
        #                         bestval = bestval,
        #                         Lbound=lbound[nextGroup],
        #                         Ubound=ubound[nextGroup],
        #                         control=list(itermax=100))$par
        predictedBest <- cma_es(fn=fname,
                                par=currentVector[nextGroup],
                                model=rf_model,
                                lower=lbound[nextGroup],
                                upper=ubound[nextGroup],
                                control=list(maxit=100,mu=2,lambda=10,
                                             sigma=0.3*max(ubound[nextGroup]-lbound[nextGroup])))$par
      }else if(infill=="ei"){
        fname=EIModel
        # predictedBest <- sansde(fname=fname,
        #                         pop=surrogate_pop,
        #                         model=rf_model,
        #                         bestmem = cv[nextGroup],
        #                         bestval = bestval,
        #                         Lbound=lbound[nextGroup],
        #                         Ubound=ubound[nextGroup],
        #                         control=list(itermax=100),minVal=bestval)$par
        predictedBest <- cma_es(fn=fname,
                                par=currentVector[nextGroup],
                                model=rf_model,
                                lower=lbound[nextGroup],
                                upper=ubound[nextGroup],
                                minVal=bestval,
                                control=list(maxit=100,mu=2,lambda=10,
                                             sigma=0.3*max(ubound[nextGroup]-lbound[nextGroup])))$par
      }
      currentVector[nextGroup] <-  predictedBest
      newVal <- EvaluateLowerLevel(fun,
                                   currentVector,
                                   furtherGroup[[1]],furtherGroup[-1],
                                   NP_further[1],NP_further[-1],
                                   lbound=lbound,ubound=ubound,
                                   prevnEval = nEval,bestval=optimVal,
                                   fname_prefix=fname_prefix,infill=infill,
                                   ...)
      valCollection <- append(valCollection,newVal$value)
      pop <- rbind(pop,predictedBest)
      # rf_model <- randomForest::randomForest((pop),(valCollection),ntree=1000)
      rf_model <- SPOT::buildKriging((pop),matrix(valCollection))
    }
    bestIndex <- which.min(valCollection)
    bestPar <- pop[bestIndex,]

    if(min(valCollection)<bestval){
      #currentVector[nextGroup] <-  bestPar
      optimVal <- min(valCollection)
      optimPar <- bestPar
      bestval <- optimVal
    }
    print(paste0('current optim val: ',optimVal))
  }
  return(list(par=optimPar,value=optimVal,nEval=nEval))
}
