EvaluateLowerLevel_bilevel <- function(contextVector,fun,
                                       group,
                                       lbound,ubound,
                                       eval_interval = 100000,
                                       levelIndex,
                                       ...){
  control <- CMAES_control[[levelIndex]]
  prevEval <- nEval

  optimPar <- contextVector[group]
  optimVal <- fun(contextVector,...)

  funVal <- cma_es(fn=subfunctionCMA,
                   contextVector=contextVector,
                   ...,
                   par=optimPar,
                   lower=lbound,
                   upper=ubound,
                   mainfun=fun,
                   groupMember=group,
                   control=control,disableIPOP=T)
  termination_code <- funVal$termination_code

  control$vectorized <- F
  control$sigma <- funVal$sigma
  control$cov <- funVal$cov
  control$pc <- funVal$pc
  control$ps <- funVal$ps

  # CMAES_control[[levelIndex]] <<- control

  groupIndex <- levelIndex

  nEval <<- nEval + funVal$counts[1]

  if(funVal$value<optimVal){
    optimVal <- funVal$value
    optimPar <- funVal$par
  }

  previousLogged <- prevEval%/%eval_interval

  logData <- ((nEval)%/%eval_interval)-previousLogged
  if(logData){
    argList <- list(...)
    print(paste('nlog',previousLogged+1))
    save(optimVal,optimPar,cv,argList,file=paste0('datalog',previousLogged+1,'.Rdata'))
  }
  return(list(lowerPar=optimPar,value=optimVal,newControl=control))
}


EvaluateUpperLevel_bilevel <- function(x,fun,
                                       upperGroup,
                                       lowerGroup,
                                       levelIndex,
                                       lbound,ubound,
                                       prevnEval=0,eval_interval=100000,
                                       fname_prefix="datalog_",
                                       ...){
  if((any(x>ubound))){
    message('constraint violated - ubound')
  }else if((any(x<lbound))){
    message('constraint violated - lbound')
  }

  cv_x <- cv # copy the context vector
  cv_x[upperGroup] <- x # replace the upper group with the current x values.

  # get the objective value by optimizing the lower level
  res <- EvaluateLowerLevel_bilevel(contextVector = cv_x,
                                    fun=fun,
                                    group=lowerGroup,
                                    lbound=lbound[lowerGroup],
                                    ubound=ubound[lowerGroup],
                                    levelIndex = levelIndex+1,
                                    ...)


  if(res$value < bestval_lower){
    print(paste('lower level optimized:',res$value))
    print(nEval)
    record$nEval <<- append(record$nEval ,nEval)
    record$conv <<- append(record$conv ,res$value)

    save(history,res,cv_x,file='history_MLLSGO_4group_BLCMA.Rdata')

    bestval_lower <<- res$value
    cv[lowerGroup] <<- res$lowerPar
    CMAES_control[[levelIndex+1]] <<- res$newControl
  }
  # res$value is the optimum value for the current level, subject to the optimum lower level
  return(res$value)
}
