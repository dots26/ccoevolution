# EvaluateLowerLevel_bilevel <- function(contextVector,fun,
#                                        group,
#                                        lbound,ubound,
#                                        eval_interval = 100000,
#                                        levelIndex,
#                                        ...){
#   control <- CMAES_control[[levelIndex]]
#   prevEval <- nEval
#
#   optimPar <- contextVector[group]
#   optimVal <- fun(contextVector,...)
#
#   funVal <- cma_es(fn=subfunctionCMA,
#                    contextVector=contextVector,
#                    ...,
#                    par=optimPar,
#                    lower=lbound,
#                    upper=ubound,
#                    mainfun=fun,
#                    groupMember=group,
#                    control=control,disableIPOP=T)
#   termination_code <- funVal$termination_code
#
#   control$vectorized <- F
#   control$sigma <- funVal$sigma
#   control$cov <- funVal$cov
#   control$pc <- funVal$pc
#   control$ps <- funVal$ps
#
#   # CMAES_control[[levelIndex]] <<- control
#
#   groupIndex <- levelIndex
#
#   nEval <<- nEval + funVal$counts[1]
#
#   if(funVal$value<optimVal){
#     optimVal <- funVal$value
#     optimPar <- funVal$par
#   }
#
#   previousLogged <- prevEval%/%eval_interval
#
#   logData <- ((nEval)%/%eval_interval)-previousLogged
#   if(logData){
#     argList <- list(...)
#     print(paste('nlog',previousLogged+1))
#     save(optimVal,optimPar,cv,argList,file=paste0('datalog',previousLogged+1,'.Rdata'))
#   }
#   return(list(lowerPar=optimPar,value=optimVal,newControl=control))
# }
#
#
# EvaluateUpperLevel_bilevel <- function(x,fun,
#                                        upperGroup,
#                                        prevnEval=0,eval_interval=100000,
#                                        fname_prefix="datalog_",
#                                        ...){
#   if((any(x>ubound))){
#     message('constraint violated - ubound')
#   }else if((any(x<lbound))){
#     message('constraint violated - lbound')
#   }
#
#   cv_x <- cv # copy the context vector
#   cv_x[upperGroup] <- x # replace the upper group with the current x values.
#
#   res <- fun(cv_x,...)
#   # res$value is the optimum value for the current level, subject to the optimum lower level
#   return(res)
# }
