rm(list=ls())

library('ccoevolution2')
nVar <- 1000
optimum <- rep(13,nVar)
func_set <- c(f4cec,f5cec,f6cec,f7cec)
nfunc <- length(func_set)
ctrl <- list(lbound=rep(-100,nVar),ubound=rep(100,nVar),delta=rep(20,nVar))
budget <- 3000000
lbound <- rep(-100,nVar)
ubound <- rep(100,nVar)
delta <- rep(20,nVar)

group_size <- c(50,25,25,100,50,25,25,700)
for(func_index in 1){
  func <- func_set[[func_index]]
  res <- NULL
  resVal <- NULL
  for(qqIndex in 1:length(qq_mat)){
    qq <- 1#qq_mat[qqIndex]
    methods <- 1
    seed <- ((qq-1)%/%4)*1000
    groupStart <- 1
    group <- NULL
    groupEnd <- 0
    permutation <- NULL
    set.seed(seed)

    rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
    # weight <- ccoevolution:::cec_weight(7)
    weight <- c(10,1000,1,1,50,10,2)
    weight <- append(weight,0.000001)


    for(i in 1:8){
      group <- append(group,list(groupEnd+(1:group_size[i])))
      permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
      groupEnd <- groupEnd + group_size[i]
    }
    for(i in 1:7){
      rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
    }
    rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])


    if(methods==1){
      print('SACC')
      a <- SACC(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)
    }else if(methods==2){
      print('CCDG')
      a <- cc(nVar = nVar,fun=func,budget=budget,grouping_control=ctrl,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)
    }else if(methods==3){
      print('MOFBVE standard')
      a <- MOFBVE(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)
    }else if(methods==4){
      print('TSCC')
      a <- TSCC(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)
    }
    res <- cbind(res,a)
    resVal <- cbind(resVal,a$y)
    save(list=ls(),file=paste('algo',qq_mat[1],'_fun',func_index+3,'.Rdata',sep=''))
  }
}



========
  if(SA_method == 'pls') {
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

    pls <- plsr(objectiveValue~population,,data=currentPop)
    ranks <- order(coef(pls),decreasing = T)
    rankedVars <- vars[ranks]
    nLevelVar <- floor(nVar/nLevel)
    for(groupingIndex in 1:(nLevel-1)){
      start <- (groupingIndex-1)*nLevelVar+1
      end <- groupingIndex*nLevelVar
      unorderedCurrentGroup <- rankedVars[start:end]
      orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
      group <- append(group,list(orderedGroup))
    }
    start <- (nLevel-1)*nLevelVar+1
    end <- nVar
    unorderedCurrentGroup <- rankedVars[start:end]
    orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
    group <- append(group,list(orderedGroup))
  }

if(SA_method == 'lasso' || SA_method == 'ridge'){
  if(SA_method == 'lasso' ) alpha=1
  if(SA_method == 'ridge') alpha=0
  nEval <- nEval + 20000
  #population <- (randtoolbox::sobol(nEval,nVar,scrambling = 3))*(ubound-lbound)+lbound
  population <- lhs::randomLHS(nEval,nVar)
  #print(dim(population))
  objectiveValue <- fun(population,...)
  bestPopIndex <- which.min(objectiveValue)
  bestPop <- population[bestPopIndex,]
  bestObj <- min(objectiveValue)

  contextVector <- bestPop
  prevLevel <- NULL

  group <- NULL
  # prediction <- glmnet::cv.glmnet(population, objectiveValue,alpha=alpha,intercept=T)

  # if(SA_method=="lasso"){
  #   D = diag(1,nVar)
  #   prediction <- glmnet::cv.glmnet(population, objectiveValue,alpha=alpha,intercept=T)
  #   for(groupingIndex in 1:(nLevel-1)){
  #     nLevelVar <- floor(nVar/nLevel)
  #     print(paste0(SA_method,'...#',groupingIndex)) # regularization
  #
  #     stopIndex <- max(which(prediction$df<nLevelVar*groupingIndex))
  #     stopLambda <- prediction$lambda[stopIndex]
  #     indices <- coef(prediction,s=stopLambda)
  #     indices <- as.matrix(indices)
  #     activeVariable <- indices[-1]
  #     activeVariable <- which(activeVariable!=0)
  #     activeVariable <- (setdiff(activeVariable,prevLevel))
  #     prevLevel <- unique(c(prevLevel,activeVariable))
  #     group <- append(group,list(activeVariable))
  #   }
  #   lastGroup <- setdiff(1:nVar,prevLevel)
  #   group <- append(group,list(lastGroup))
  #   clusterOrder <- 1:nLevel
  # }
  if(SA_method=="ridge"){
    data <- list(population=population,objectiveValue=objectiveValue)
    prediction <- lmridge::lmridge("objectiveValue ~ population^2",data)
    print(paste0(SA_method,'...')) # regularization
    reg_coefficient <- coef(prediction,s='lambda.min')[-1]
    ranks <- order(reg_coefficient,decreasing = T)
    totalCoef <- sum(reg_coefficient)

    rankedVars <- vars[ranks]
    nLevelVar <- floor(nVar/nLevel)
    for(groupingIndex in 1:(nLevel-1)){
      start <- (groupingIndex-1)*nLevelVar+1
      end <- groupingIndex*nLevelVar
      unorderedCurrentGroup <- rankedVars[start:end]
      orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
      group <- append(group,list(orderedGroup))
    }
    start <- (nLevel-1)*nLevelVar+1
    end <- nVar
    unorderedCurrentGroup <- rankedVars[start:end]
    orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
    group <- append(group,list(orderedGroup))
  }
}
