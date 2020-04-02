# Interconnection
# groupMember <- 1:nVar
# mu<-100
# best <- cma_es(contextVector[groupMember],
#                fn = subfunctionCMA,
#                contextVector=contextVector,
#                groupMember=groupMember,
#                mainfun=fun,
#                ...,
#                lower = lbound,
#                upper=ubound,
#                control = list(vectorized=T,
#                               mu=mu,lambda=mu,
#                               maxit=90,
#                               sigma=0.3*max(ubound-lbound),
#                               diag.value=T))
# nlogging_this_layer <- floor((nEval+best$counts[1])/evalInterval)-floor(nEval/evalInterval)
#
# if(nlogging_this_layer>0){
#   for(i in 1:nlogging_this_layer){
#     nEval_to_logging <- (evalInterval*i) - nEval%%evalInterval
#     nGeneration_to_consider <- floor(nEval_to_logging/mu)
#     print(nGeneration_to_consider,nlogging_this_layer)
#     print(best$diagnostic)
#     if(!is.matrix(best$diagnostic$value)){
#       best$diagnostic$value <- matrix(best$diagnostic$value)
#     }
#     bestObj_logging <- min(best$diagnostic$value[1:nGeneration_to_consider,])
#     convergence_history <- append(convergence_history,min(bestObj_logging,convergence_history[length(convergence_history)],bestObj))
#   }
# }
# nEval <- nEval + best$counts[1]
#
# if((budget-nEval)>0){ # only update if it doesnt exceed budget
#   if(!is.null(best$par)){
#     contextVector <- best$par
#     obj <- best$value
#     if(obj < bestObj){
#       bestPop <- contextVector
#       bestObj <- obj
#     }
#   }
# }else{
#   break
# }
