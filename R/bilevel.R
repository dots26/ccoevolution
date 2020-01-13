EvaluateLowerLevel_bilevel <- function(fun,currentVector,
                                       nextGroup,
                                       controlNext,
                                       lbound,ubound,
                                       prevnEval=0,eval_interval=100000,
                                       maxIter=300,fname_prefix="datalog_",
                                       bestval=Inf,
                                       ...){
  optimPar <- currentVector[nextGroup]
  optimVal <- bestval

  funVal <- cma_es(fn=subfunctionCMA,
                   contextVector=currentVector,
                   ...,
                   par=currentVector[nextGroup],
                   lower=lbound[nextGroup],
                   upper=ubound[nextGroup],
                   mainfun=fun,
                   groupMember=nextGroup,
                   logFeasible=T,
                   control=list(maxit=maxIter,vectorized=T,mu=2,lambda=5,
                                sigma=0.3*max(ubound[nextGroup]-lbound[nextGroup])))
  prevnEval <- nEval
  nEval <<- nEval + funVal$counts[1]

  if(funVal$value<optimVal){
    optimVal <- funVal$value
    optimPar <- funVal$par
  }
  if(funVal$value<gb){
    gb <<- funVal$value
    cv <<- currentVector
    cv[nextGroup] <<- optimPar
    print(paste('result',gb))
  }

  previousLogged <- prevnEval%/%eval_interval
  #
  logData <- ((nEval)%/%eval_interval)-previousLogged
  if(logData){
    argList <- list(...)
    print(paste('nlog',previousLogged+1))
    save(optimVal,optimPar,cv,gb,argList,file=paste0(fname_prefix,previousLogged+1,'.Rdata'))
  }

  return(list(par=optimPar,value=optimVal))
}


EvaluateUpperLevel_bilevel <- function(upperVector,upperGroup,fun,#currentVector,
                                       nextGroup,
                                       controlNext,
                                       lbound,ubound,
                                       prevnEval=0,eval_interval=100000,
                                       maxIter=300,fname_prefix="datalog_",
                                       bestval=Inf,
                                       ...){
  currentVector <- cv
  currentVector[upperGroup] <- upperVector
  funValue <- EvaluateLowerLevel_bilevel(fun=fun,
                                         currentVector=currentVector,
                                         nextGroup=nextGroup,
                                         lbound=lbound,
                                         ubound=ubound,
                                         prevnEval = prevnEval,
                                         eval_interval=eval_interval,bestval=bestval,
                                         maxIter=maxIter,
                                         fname_prefix=fname_prefix,
                                         ...)


  return(funValue$value)
}
