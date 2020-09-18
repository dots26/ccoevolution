#' @export
WOF_SMSEMOA <- function(population=NULL,
                        nVar=ncol(population),populationSize,fun,growPop=F,
                        maxPopulationSize=10*populationSize,
                        control=list(),nIterInner =10,
                        group_selection=c("roundrobin"),
                        ...){
  last100group <- rep(0,100)
  last100success <- rep(1,100)
  if((is.null(population) && is.null(nVar)))
    stop("Population and nVar cannot be both NULL.")
  
  if(is.vector(population))
    population <- matrix(population,nrow=1)
  
  initialPopulationSize <- populationSize
  hv_log <- NULL
  nEval <- 0
  controlParam <- function(name, default) {
    v <- control[[name]]
    
    if (is.null(v))
      return (default)
    else
      return (v)
  }
  
  ## Parameters:
  vectorized <- controlParam("vectorized",F)
  scaleinput <- controlParam("scaleinput",T)
  group <-  controlParam("group", NULL)
  importance <- controlParam('importance',NULL)
  groupmethod  <- controlParam("groupmethod", 'morris')
  pm <- controlParam("mutationrate", min(c(length(group)/nVar,1)))
  pc <- controlParam("crossoverrate", 1)
  hypervolumeMethod <- controlParam('hvmethod','exact')
  hypervolumeMethodParam <- controlParam('hvmethodparam',list())
  refPoint <- controlParam('refpoint',NULL) # if not NULL and loghv is TRUE, HV history w.r.t to this point will be recorded
  
  loghv <- controlParam('loghv',F) # will log HV if refpoint is supplied
  if(loghv)
    hv_log <- NULL
  # if dynamic, prescribed ref will be used for HVhistory, but the SMSEMOA will adapt the refpoint location
  refpointmethod <- controlParam("refpointmethod",c("dynamic","static"))
  computationRefPoint <- refPoint
  if(refpointmethod[1]=="dynamic")
    computationRefPoint <- NULL
  lbound <- controlParam('lower',rep(0,nVar))
  ubound <- controlParam('upper',rep(1,nVar))
  md <- controlParam( 'mutationdistribution',20)
  cd <- controlParam( 'crossoverdistribution',30)
  alpha <- controlParam( 'alpha',1.5)
  initialpop <-  controlParam( 'initialpop',c('lhs','group'))
  gctrl <- controlParam('groupingcontrol',list(r=20))
  budget <- controlParam('budget',1e5)
  
  gctrl$lbound <- lbound
  gctrl$ubound <- ubound
  
  scale_multip <- 1
  scale_shift <- 0
  
  if(is.null(population)){
    if(initialpop[1]=='lhs'){
      # scaled
      population <- t(InitializePopulationLHS(ncol = populationSize,
                                              nrow = nVar,
                                              minVal = lbound,
                                              maxVal = ubound))
      rescaled_pop <-t (t(population)*scale_multip + scale_shift)
      populationObjective <- fun(rescaled_pop,...)
    }
  }else{
    currentpopulationSize <- nrow(population) # already scaled
    nVar <- ncol(population) # row major!
    
    generateMore <- NULL
    if(populationSize-currentpopulationSize>0)
      generateMore <- t(InitializePopulationLHS(ncol = populationSize-currentpopulationSize,
                                                nrow = nVar,
                                                minVal = lbound,
                                                maxVal = ubound))
    population <- rbind(population,generateMore)
    rescaled_pop <- t (t(population)*scale_multip + scale_shift)
    
    populationObjective <- fun(rescaled_pop,...)
  }
  
  scaledRefPoint <- (refPoint- scale_shift) / scale_multip
  vars <- 1:nVar
  
  # grouping
  # the grouping does not use the "scaleinput" parameter
  if(is.null(group)){
    if(groupmethod=='DG2'){
      a <- DG2(length(group[[i]]),
               fun,
               control=gctrl,...)
      group <- append(a$group,list(a$separable))
      nEval <- nEval + a$nEval
    }
    if(groupmethod=='morris'){
      a <- sensitivity::morris(model=fun,
                               factors=nVar,
                               r = gctrl$r,
                               design = list(type='oat',levels=8,grid.jump=4),
                               binf=lbound,
                               bsup=ubound,
                               scale=F,...)
      
      nEval <- nEval + r*(nVar+1)
      mu.star <- apply(a$ee, 2, function(a) mean(abs(a)))
      sigma <- apply(a$ee, 2, sd)
      
      ranks <- order(mu.star,decreasing = T)
      
      reg_coefficient <- mu.star
      totalCoef <- sum(reg_coefficient)
      
      rankedVars <- vars[ranks]
      nLevelVar <- floor(nVar/nLevel)
      for(groupingIndex in 1:(nLevel-1)){
        start <- (groupingIndex-1)*nLevelVar+1
        end <- groupingIndex*nLevelVar
        unorderedCurrentGroup <- rankedVars[start:end]
        group <- append(group,list(unorderedCurrentGroup))
        groupWeight <- append(groupWeight,sum(reg_coefficient[unorderedCurrentGroup]))
      }
      start <- (nLevel-1)*nLevelVar+1
      end <- nVar
      unorderedCurrentGroup <- rankedVars[start:end]
      groupWeight <- append(groupWeight,sum(reg_coefficient[unorderedCurrentGroup]))
      group <- append(group,list(unorderedCurrentGroup))
      clusterOrder <- 1:nLevel
      importance <- reg_coefficient
    }
  }
  
  if(is.null(importance)){
    importance <- rep(1,nVar)
  }
  
  if(initialpop[1]=='group'){
    nIndivGroup <- nrow(a$x)
    takenSample <- sample(nIndivGroup,populationSize)
    population <- a$x[takenSample,]
    populationObjective <- a$y[takenSample,]
  }
  
  nObjective <- ncol(populationObjective)
  nGroup <- length(group)
  
  budget_left <- budget - nEval
  success_off <- NULL
  total_off <- 0
  total_success <- 0
  nTimesGrow <- 0
  growingPeriod <- budget/round((maxPopulationSize/initialPopulationSize)-1)
  
  groupSuccess <- rep(1,nGroup)
  denom <- rep(1,nGroup)
  groupFitness <- groupSuccess/denom
  groupIndex <- 0
  while(budget_left>0){
    ## main SMS-EMOA -- not WOF
    groupIndex <- groupIndex%%length(group)+1
    
    print("SMS")
    for (iter in 1:nIterInner)
    { 
      # print(iter)
      parentIndex <- sample(1:populationSize,2,replace = FALSE)
      
      # call SMS-EMOA here
      
      optimRes <- MaOEA::SMSEMOA(population = t(population), # transposed because it should be col major
                                 fun = fun,
                                 nObjective = nObjective,
                                 control = control,
                                 ...)
      population <- t(optimRes$population) # here, row major
      popObj <- t(optimRes$objective)
      
      plot(0,0,xlim=c(0,4),ylim=c(0,8))
      
      points(popObj)
      rankss <- nsga2R::fastNonDominatedSorting(popObj)
      points(popObj[rankss[[1]],,drop=F],col="purple")
      # points(popObj[removedPoint,,drop=F],col="blue")
      
      # newPointSurvive <- optimRes$successfulOffspring
     
      if((loghv && (!is.null(refPoint)))){
        hv_log <- append(hv_log,GetHypervolume(t(popObj),
                                               reference = refPoint,
                                               method = hypervolumeMethod))
      }
      
      nEval <- nEval + 1
      budget_left <- budget - nEval
    }
    
    
    rnk <- nsga2R::fastNonDominatedSorting((popObj))
    firstFront <- t(popObj[rnk[[1]],])
    refInFirstFront <- which.max(MaOEA::GetHVContribution(firstFront))
    chosenReference <- rnk[[1]][refInFirstFront]

    # WOF
    # check distance to boundaries
    direction_vector <- (population[chosenReference,]-lbound)
    unit_dir_vector <- direction_vector/norm(direction_vector,"2")
    
    dir_to_ubound <- (ubound-lbound)
    
    w_lbound <- rep(-1,nGroup)
    scale_to_ubound <- NULL
    for(groupIndex in 1:nGroup){
      groupMember <- group[[groupIndex]]
      thisGroupUpper <- min(dir_to_ubound[groupMember]/direction_vector[groupMember])-1
      scale_to_ubound <- append(scale_to_ubound,thisGroupUpper)
    }

    w_ubound <- scale_to_ubound
    
    w_control <- control 
    w_control$lower <- w_lbound
    w_control$upper <- w_ubound
    
    w_population <- t(InitializePopulationLHS(ncol = 10,
                                              nrow = nGroup,
                                              minVal = w_lbound,
                                              maxVal = w_ubound))
    print("WOF")
    
    for (iter in 1:nIterInner)
    {
      # print(iter)
      # rnk <- nsga2R::fastNonDominatedSorting((popObj))
      # objRange <- apply(popObj[,(varNo+1):(varNo+objDim)], 2, max) -
      #   apply(popObj[,(varNo+1):(varNo+objDim)], 2, min)
      
      # chosenReference <- which.max(nsga2R::crowdingDist4frnt(pop=population,
      #                                              rnk = rnk,
      #                                              rng = objRange))
 
      # call SMSEMOA for weight here
      optimRes <- MaOEA::SMSEMOA(population = t(w_population),
                                 fun = weight_transformed_function, 
                                 mainfun=fun,
                                 nObjective = nObjective,
                                 control = w_control,
                                 original_ubound = ubound, original_lbound=lbound,
                                 x=population[chosenReference,],
                                 group=group,
                                 ...)
      w_population <- t(optimRes$population)

      # print(MaOEA::GetHypervolume(optimRes$objective,reference = c(2.2,4.4)))
      nEval <- nEval + 1
      budget_left <- budget - nEval
    }
    print("ending cycle")
    newIndividual <- population[chosenReference,]

    largest_contrib <- which.max(MaOEA::GetHVContribution(optimRes$objective))
    bestObj <- optimRes$objective[,largest_contrib]
    best_w <- w_population[largest_contrib,]
    for(groupIndex in 1:nGroup){
      groupMember <- group[[groupIndex]]
      newIndividual[groupMember] <- (best_w[groupIndex])*
        direction_vector[groupMember]+
        population[chosenReference,groupMember]
    }
    
    combinedPopulation <- rbind(population,newIndividual) # row major
    combinedObj <- rbind((popObj),bestObj) # row major
 
    rnk <- nsga2R::fastNonDominatedSorting((combinedObj))
    worstFrontIndex <- rnk[[length(rnk)]]
    
    if(length(worstFrontIndex)==1){
    }else{
      leastContrib <- GetLeastContributor(t(combinedObj[worstFrontIndex,]))
      removedPoint <- worstFrontIndex[leastContrib]
    }
    population <- combinedPopulation[-removedPoint,]
    populationObjective <- combinedObj[-removedPoint,]
    print(paste("this cycle best HV:",MaOEA::GetHypervolume(t(combinedObj[-removedPoint,]),reference =  c(2.2,4.4))))
    plot(0,0,xlim=c(0,4),ylim=c(0,8))
    
    points(combinedObj)
    rankss <- nsga2R::fastNonDominatedSorting(combinedObj)
    points(combinedObj[rankss[[1]],],col="red")
    points(combinedObj[removedPoint,,drop=F],col="blue")
  }
  rescaled_pop <-(population)*scale_multip + scale_shift ## check scaling
  
  
  return(list(x=t(rescaled_pop),y=populationObjective,hv_log=hv_log))
}

# optimize w, with reference vector x, on function mainfun
# w should be column major, x is a vector
# row_major indicate if the mainfun accept row major input. default to F (col major)
weight_transformed_function <- function(w,x,
                                        original_ubound,
                                        original_lbound,
                                        group,mainfun,row_major=F,...){
  nGroup <- length(group)
  x <- matrix(x,nrow=length(x))
  direction_vector <- matrix(x-original_lbound)
  new_x <- x
  for(groupIndex in 1:nGroup){
    groupMember <- group[[groupIndex]]

    new_x[groupMember,] <- x[groupMember,,drop=F]+
      direction_vector[groupMember,,drop=F]%*%(w[groupIndex,])
  }
  if(row_major)
    new_x <- t(new_x) # new_x is transposed before being used as input to mainfun
  objective <- mainfun(new_x,...) 
} 


## will be used in CC-SMS-EMOA based algorithms
# subfunction <- function(x,mainfun,contextVector,groupMember...){
#   
# }
