CC_SMSEMOA <- function(CVs,fun,nObjective,control=list()){
  if(is.vector(CVs))
    CVs <- matrix(CVs,nrow=1)
  nVar <- ncol(CVs) # row major!

  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v))
      return (default)
    else
      return (v)
  }

  ## Parameters:
  group <-  controlParam("group", NULL)
  importance <- controlParam('importance',NULL)
  groupmethod  <- controlParam("groupmethod", 'DG2')
  pm <- controlParam("mutationrate", 1)
  pc <- controlParam("crossoverrate", 1)
  hypervolumeMethod <- controlParam('hvmethod','exact')
  hypervolumeMethodParam <- controlParam('hvmethodparam',list())
  refPoint <- controlParam('refpoint',NULL)
  lbound <- controlParam('lower',rep(0,nVar))
  ubound <- controlParam('upper',rep(1,nVar))
  md <- controlParam( 'mutationdistribution',20)
  cd <- controlParam( 'crossoverdistribution',30)

  gctrl <- controlParam('groupingcontrol',list(r=20))

  gctrl$lbound <- lbound
  gctrl$ubound <- ubound

  # grouping
  if(is.null(group)){
    if(groupmethod=='DG2'){

      group <- DG2(length(group[[i]]),
                   fun,
                   control=gctrl,...)
    }
    if(groupmethod='morris'){
      a <- sensitivity::morris(model=fun,
                               factors=nVar,
                               r = gctrl$r,
                               design = list(type='oat',levels=8,grid.jump=4),
                               binf=lbound,
                               bsup=ubound,
                               scale=F,...)

      bestPopIndex <- which.min(a$y)
      bestPop <- a$X[bestPopIndex,]
      bestObj <- min(a$y)

      contextVector <- bestPop
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

  for(groupIndex in 1:nGroup)
  MaOEA::SMSEMOA()

}
