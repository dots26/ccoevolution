#' @export
CC_SMSEMOA <- function(population=NULL,
                       nVar=ncol(population),populationSize,fun,growPop=F,
                       maxPopulationSize=10*populationSize,
                       control=list(),nIterInner=1,
                       group_selection=c("adaptive","random","roundrobin","adaptive100"),
                       ...){
  cdlog <- NULL
  mdlog <- NULL
  
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
  scaleinput <- controlParam("scaleinput",F)
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
  lbound <- controlParam('lower',rep(0,nVar))
  ubound <- controlParam('upper',rep(1,nVar))
  md <- controlParam( 'mutationdistribution',20)
  cd <- controlParam( 'crossoverdistribution',30)
  step_adapt <- controlParam('stepadapt',F)
  alpha <- controlParam( 'alpha',1.5)
  initialpop <-  controlParam( 'initialpop',c('lhs','group'))
  gctrl <- controlParam('groupingcontrol',list(r=20))
  budget <- controlParam('budget',1e5)
  
  gctrl$lbound <- lbound
  gctrl$ubound <- ubound
  
  scale_multip <- 1
  scale_shift <- 0
  if(scaleinput){
    scale_shift <- -lbound
    scale_multip <- (ubound-lbound)
    
    ubound <- rep(1,nVar)
    lbound <- rep(0,nVar)
    if(!is.null(population))
      population <- t((t(population) - scale_shift) / scale_multip) # scale available population
  }
  
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
      generateMore <- t(InitializePopulationLHS(numberOfIndividuals = populationSize-currentpopulationSize,
                                                chromosomeLength = nVar,
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
  nIter <- rep(nIterInner,nGroup) # number of SMSEMOA iteration on each group
  
  budget_left <- budget - nEval
  success_off <- NULL
  total_off <- 0
  total_success <- 0
  nTimesGrow <- 0
  growingPeriod <- budget/round((maxPopulationSize/initialPopulationSize)-1)
  # 17.06 introduce group fitness
  groupSuccess <- rep(1,nGroup)
  denom <- rep(1,nGroup)
  groupFitness <- groupSuccess/denom
  groupIndex <- 0
  while(budget_left>0){
    # growing population
    while(nEval%/%growingPeriod>nTimesGrow){
      parentIndex <- (1:populationSize)
      offspring <- nsga2R::boundedSBXover(parent_chromosome = population[parentIndex,],
                                          lowerBounds = lbound,
                                          upperBounds =  ubound,
                                          cprob = pc,mu = cd)
      offspring <- nsga2R::boundedPolyMutation(parent_chromosome = offspring,
                                               lowerBounds = lbound,
                                               upperBounds =  ubound,
                                               mprob=pm,mum = md)
      newPopSize <- populationSize+initialPopulationSize
      offspring <- offspring[sample(1:populationSize,
                                    initialPopulationSize,replace = F),]
      rescaled_offspring <- (offspring*scale_multip + scale_shift)
      offspringObjectiveValue <- fun(rescaled_offspring,...)
      
      population <- rbind(population,offspring)
      populationObjective <- rbind(populationObjective,offspringObjectiveValue)
      populationSize <- populationSize + initialPopulationSize
      nTimesGrow <- nTimesGrow +1
      print(paste("popsize:",populationSize))
    }
    if(group_selection == "adaptive")
      groupIndex <- sample(length(group),1,prob = groupFitness)
    if(group_selection == "adaptive100")
      groupIndex <- sample(length(group),1,prob = groupFitness)
    if(group_selection == "random")
      groupIndex <- sample(length(group),1)
    if(group_selection == "roundrobin")
      groupIndex <- groupIndex%%length(group)+1
    
    last100group <- last100group[-1]
    last100group <- append(last100group,groupIndex)
    for(iter in 1:nIter[groupIndex])
    {  
      xx <- (0:1000)/1000*2
      yy <- 2*(4-xx*xx)^0.5
      
      #TODO determine CV
      groupMember <- group[[groupIndex]]
      groupSize <- length(groupMember)
      
      offspring<-matrix(,ncol = nVar,nrow = 0)
      offspringObjectiveValue<-matrix(,ncol = nObjective,nrow=0)
      
      randomParent <- (runif(1)<1)
      parentIndex <- vector(length = 2)
      if(randomParent){
        # print("random parent")
        parentIndex <- sample(1:populationSize,2,replace = FALSE)
      }else{
        # print("exploring parent")
        bestInDim <- sample(1:nObjective,1,replace=F)
        parentIndex[1] <- which.min(populationObjective[,bestInDim[1]])
        parentIndex[2] <- sample(c(1:populationSize)[-parentIndex[1]],1)
      }
      #Crossover
      offspring <- nsga2R::boundedSBXover(parent_chromosome = population[parentIndex,groupMember],
                                          lowerBounds = lbound[groupMember],
                                          upperBounds =  ubound[groupMember],
                                          cprob = pc,mu = cd)
      offspring <- matrix(offspring[sample(1:2,1),],nrow=1, ncol=groupSize)
      #Mutation
      offspring <- nsga2R::boundedPolyMutation(parent_chromosome = offspring,
                                               lowerBounds = lbound[groupMember],
                                               upperBounds =  ubound[groupMember],
                                               mprob=pm,mum = md)
      # evaluate objective
      this_offspring <- offspring[1,,drop=FALSE]
      
      #one parent
      closestCV <- population[parentIndex[1],]
      farthestCV <- population[parentIndex[2],]
      #rescaled_offspring <- this_offspring*scale_multip + scale_shift
      closestCV[groupMember] <- this_offspring
      farthestCV[groupMember] <- this_offspring
      if(runif(1)<1){
        rescaled_offspring <- t(matrix(closestCV*scale_multip + scale_shift))
      }else{
        rescaled_offspring <- t(matrix(farthestCV*scale_multip + scale_shift))
      }
      offspringObjectiveValue <- fun(rescaled_offspring,...)
      
      combinedPopulation <- rbind(population,closestCV)
      combinedObjectiveValue <- rbind(populationObjective,offspringObjectiveValue)
      
      # nondominated sorting
      rnk <- nsga2R::fastNonDominatedSorting(combinedObjectiveValue)
      
      rnkIndex <- integer(ncol(combinedPopulation))
      sortingIndex <- 1;
      while (sortingIndex <= length(rnk)) {
        rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
        sortingIndex <- sortingIndex + 1;
      }
      
      # fill new population
      newPopulation <- combinedPopulation
      newPopulationObjective <- combinedObjectiveValue
      
      worstFrontIndex <- rnk[[length(rnk)]]
      
      # individualIndexTobeChecked <- NULL
      newPointSurvives <- 1
      
      if(length(rnk[[length(rnk)]])>1){
        # normalize objectives to 0-1, minimization problem, ideal point is at (0,0)
        normalizationList <- MaOEA::Normalize(t(combinedObjectiveValue),refPoint)
        normalizedObjective <- normalizationList$normalizedObjective
        individualIndexTobeChecked <- worstFrontIndex
        
        if(refpointmethod[1]=="static") {
          smallestContributor <- GetLeastContributor(normalizedObjective[,worstFrontIndex],
                                                     reference=normalizationList$transformedReference,
                                                     method = hypervolumeMethod,
                                                     hypervolumeMethodParam = hypervolumeMethodParam)
        }else{
          # print(normalizedObjective[,worstFrontIndex])
          smallestContributor <- GetLeastContributor(normalizedObjective[,worstFrontIndex],
                                                     method = hypervolumeMethod,
                                                     hypervolumeMethodParam = hypervolumeMethodParam)
        }
        
        removedPoint <- individualIndexTobeChecked[smallestContributor]
        
        newPopulation <- newPopulation[-removedPoint,]
        newPopulationObjective <- newPopulationObjective[-removedPoint,]
        
        if(removedPoint == nrow(combinedPopulation))
          newPointSurvives <- 0
      }else{
        removedPoint <- rnk[[length(rnk)]][1]
        if(removedPoint == nrow(combinedPopulation))
          newPointSurvives <- 0
        newPopulation <- newPopulation[-(rnk[[length(rnk)]][1]),]
        newPopulationObjective <- newPopulationObjective[-(rnk[[length(rnk)]][1]),]
      }
      success_off <- append(success_off,newPointSurvives)
      
      if(length(success_off)>100){
        success_off <- success_off[-1]
      }
      total_off <- total_off+1
      denom[groupIndex] <- denom[groupIndex] + 1
      thisGroupInLast100 <- which(last100group==groupIndex)
      denom100 <- length(thisGroupInLast100)
      
      if(step_adapt){
        if(newPointSurvives){
          total_success <- total_success+1
          cd <- max(c((cd +1)/alpha-1,0))
          md <- max(c((md +1)/alpha-1,15))
          last100success <- last100success[-1]
          last100success <- append(last100success,1)
          
          groupSuccess[groupIndex] <- groupSuccess[groupIndex]+1
        }else{
          cd <- min(c((cd +1)*alpha-1,50))
          md <- min(c((md +1)*alpha-1,30))
          last100success <- last100success[-1]
          last100success <- append(last100success,0)
        }
      }
      # print(paste("cd",cd))
      cdlog <- append(cdlog,cd)
      mdlog <- append(mdlog,md)
      checklast100success <- sum(last100success[thisGroupInLast100])
      if(group_selection == "adaptive")
        groupFitness[groupIndex] <- groupSuccess[groupIndex]/denom[groupIndex]
      if(group_selection == "adaptive100")
        groupFitness[groupIndex] <- (checklast100success+1)/denom100
      # print(paste(' Success Rate:',mean(success_off)))
      populationObjective <- newPopulationObjective
      
      if((loghv && (!is.null(refPoint)))){
        hv_log <- append(hv_log,GetHypervolume(t(populationObjective),
                                               reference = refPoint,
                                               method = hypervolumeMethod))
      }
      population <- newPopulation
      nEval <- nEval + 1
      budget_left <- budget - nEval
    }
    # plot(xx,yy,xlim=c(0,6),ylim=c(0,10),type="l")
    
    # points(t(t(populationObjective)))
    # points(combinedObjectiveValue[nrow(combinedObjectiveValue),,drop=F],col='blue')
    # points(combinedObjectiveValue[removedPoint,,drop=F],col='red')
  }
  hv_log <- append(hv_log,GetHypervolume(t(populationObjective),
                                         reference = refPoint,
                                         method = hypervolumeMethod))
  rescaled_pop <-t (t(population)*scale_multip + scale_shift)
  save(cdlog,mdlog,file='adapt_step.Rdata')
  return(list(x=rescaled_pop,y=populationObjective,hv_log=hv_log))
}
