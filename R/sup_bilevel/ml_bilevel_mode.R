multilevel_bilevel <- function(group,contextVector,mainfun,NP=NULL,
                               ubound,lbound,
                               budget=3000000,eval_interval=100000,
                               fname_prefix="datalog_",control=NULL,
                               bestval=Inf,maxSeqIter=10,infill=c('best','ei'),
                               mode='optim',...){
  nGroup <- length(group)
  infill <- infill[1]
  if(is.null(NP)){
    NP <- rep(20,nGroup)
  }
  nVar <- length(contextVector)
  groupLength <- pracma::zeros(1,nGroup)
  for(i in 1:nGroup){
    groupLength[i] <- length(group[[i]])
  }
  optimVal <- bestval
  optimPar <- contextVector

  # build and collect surrogates
  # by optimizing lower levels
  budget_left <- budget-nEval
  while(budget_left > 0){
    for(levelIndex in nGroup:2){
      print(paste('level',levelIndex))
      maxIter <- 100

      valCollection <- NULL
      upperLevel <- levelIndex - 1
      upperGroup <- group[[upperLevel]]
      currentGroup <- group[[levelIndex]]
      currentVector <- contextVector

      if(mode=='optim'){
        pop <- NULL
        maxit_upper <- 10
        mu_upper <- 5
        lambda_upper <- 5
        # optimize upper level
        # inside, optimize the lower level for each suggested upper level vector
        funValue <- cma_es(par=currentVector[upperGroup],
                           fn = EvaluateUpperLevel_bilevel,
                           upperGroup=upperGroup,
                           lower = lbound[upperGroup],
                           upper = ubound[upperGroup],
                           fun=mainfun,
                           nextGroup = group[[levelIndex]],
                           lbound=lbound,
                           ubound=ubound,
                           prevnEval = nEval,
                           eval_interval=eval_interval,
                           maxIter=maxIter,
                           bestval=Inf,
                           control=list(maxit=maxit_upper,
                                        sigma=0.3*max((lbound-ubound)[upperGroup]), #local search
                                        mu=mu_upper,lambda=lambda_upper,
                                        diag.pop=T,diag.value=T),
                           ...)

        print(funValue$value)
        save(funValue,cv,file = 'TT.Rdata')

        # collect all individual and value in the upper level optimization
        for(iter_upper in 1:maxit_upper){
          for(mu_upper_index in 1:mu_upper){
            valCollection <- append(valCollection,funValue$diagnostic$value[iter_upper,mu_upper_index])
            pop <- rbind(pop,funValue$diagnostic$pop[,mu_upper_index,iter_upper])
          }
        }

        if(funValue$value < optimVal){
          optimPar[upperGroup] <- funValue$par
          optimVal <- funValue$value
          currentVector <- cv
        }
      }else{
        range <- ubound[upperLevel]-lbound[upperLevel]
        dist_to_ubound <- ubound[upperLevel]-currentVector[upperGroup]
        dist_to_lbound <- currentVector[upperGroup]-lbound[upperLevel]

        pop <- InitializePopulationSobol(groupLength[upperLevel],NP[upperLevel],
                                         minVal = max(currentVector[upperGroup]-range*0.3,lbound[upperLevel]),
                                         maxVal = min(currentVector[upperGroup]+range*0.3,lbound[upperLevel]))
        valCollection <- NULL
        for(i in 1:NP[upperLevel]){
          print(paste('level',levelIndex-1,'individual',i))
          currentVector[upperGroup] <- pop[i,]
          # optimize the lower level for each upper level
          funValue <- EvaluateLowerLevel_bilevel(mainfun,
                                                 currentVector,
                                                 group[[levelIndex]],
                                                 lbound=lbound,ubound=ubound,
                                                 prevnEval = nEval,
                                                 eval_interval,bestval=Inf,maxIter=maxIter,
                                                 ...)
          valCollection <- append(valCollection,funValue$value)
          if(funValue$value < optimVal){
            optimPar[upperGroup] <- pop[i,]
            optimPar[currentGroup] <- funValue$par
            optimVal <- funValue$value
            currentVector <- cv
          }
        }
      }
      # build surrogate/mapping based on valCollection and pop
      if(maxSeqIter>0){
        print('seq optim')
        for(sequentialIter in 1:maxSeqIter){
          # print(paste('Iter',sequentialIter,nEval))
          rf_model <- SPOT::buildKriging(pop,matrix(valCollection))
          # rf_model <- randomForest::randomForest(pop,valCollection,ntree=1000)
          if(infill=="best"){
            funname <- bestPredModel.randomForest
          }else if(infill=='ei'){
            funname <- EIModel.kriging
          }
          # predict the best design for upper level
          predictedBest <- cma_es(fn=funname,
                                  par=currentVector[upperGroup],
                                  model=rf_model,
                                  lower=lbound[upperGroup],
                                  upper=ubound[upperGroup],
                                  control=list(maxit=200,mu=2,lambda=10,
                                               sigma=0.3*max(ubound[upperGroup]-lbound[upperGroup])))$par

          currentVector[upperGroup] <-  predictedBest
          # verify
          newVal <- EvaluateLowerLevel_bilevel(mainfun,
                                               currentVector,
                                               currentGroup,
                                               lbound=lbound,ubound=ubound,
                                               prevnEval = nEval,bestval=optimVal,
                                               fname_prefix=fname_prefix,maxIter=maxIter,
                                               ...)
          if(newVal$value<optimVal){
            optimPar[upperGroup] <- predictedBest
            optimPar[currentGroup] <- newVal$par
            optimVal <- newVal$value
          }
          pop <- rbind(pop,predictedBest)
          valCollection <- append(valCollection,newVal$value)
        }
      }
    }
    budget_left <- budget-nEval
  }

  return(list(par=optimPar,value=optimVal))
}
