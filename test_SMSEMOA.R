rm(list=ls())
library(MaOEA)
qq <- 9
#setwd("~/Desktop/MOLSGO/ccoevolution")
files.sources = list.files(path = "R",full.names = T,pattern = ".R")
sapply(files.sources, source)
source('local_mutation_gaussSMSEMOA.R')
set.seed(100)
nIterInner <- c(1,250,500,1000)
chosennIter <- (qq-1)%%4+1
nIterInner <- nIterInner[chosennIter]
func <- MaOEA::WFG4
adapt_step_log <- list()

m <- 2
evalInterval <- 1e5
vargroup <- vector(mode = 'list',length = m)
a <- vector(mode = 'list',length = m)
r <- 4
nVar <- 1000
lbound <- rep(0,nVar)
ubound <- (1:nVar)*2
nEval <- 0
k <- 250
vars <- 1:nVar

transposedFun <- function(x,mainfunc,...){
  return(t(mainfunc(t(x),...)))# transpose needed because MaOEA is set to work with x in column major order while morris and ccoevolution in row major
}
morrisMultOut <- function (model = NULL, factors, r = r,
                           design = list(type = "oat",
                                         levels = 5, grid.jump = 3), binf = 0, bsup = 1, scale = TRUE,
                           ...) {
  M = sensitivity::morris(model = NULL, factors = factors, r = r, design = design,
                          binf = binf, bsup = bsup, scale = scale, ...)
  print("Test")
  if (!is.null(model)) {
    Y = model(M$X,...)
    M = sensitivity:::.morrisMultOut(t(Y), M)
    M$objs <- Y
  }
  class(M) = c("morrisMultOut", "morris")
  return(M)
}
timelog <- NULL
abc_log <- vector('list')

for (testID in (qq-1)%/%4+1)
{
  testIDinUse <- testID+1
  
  nLevel <- c(1,50,10,4,50,10,4,50,10,4,50,10,4,50,10,4,50,10,4)
  nLevel <- nLevel[testIDinUse]
  budget <- 10000
  popSize <- 100
  
  print('SAs')
  # consider an SVD of the output
  SA <- morrisMultOut(model=transposedFun,mainfunc=func,k=k,nObj=m,
                      factors=nVar,
                      r = r,
                      design = list(type='oat',levels=8,grid.jump=4),
                      binf=lbound,
                      bsup=ubound,
                      scale=T)
  nEval <- nEval + r*(nVar+1)
  leftBudget <- budget - nEval
  nRepeat <- floor(nEval/evalInterval)
  if(nRepeat>0){
    for(i in 1:nRepeat){
      convergence_history <- append(convergence_history,bestObF)
    }
  }
  
  
  samples <- sample(1:nrow(SA$X),popSize)
  x <- SA$X[samples,]
  objs <- SA$objs[samples,]
  
  reg_coefficient <- vector('list')
  group <- NULL
  groupWeight <- NULL
  mu.star <- apply(SA$ee, 2, function(x) mean(abs(x)))
  
  ranks <- order(mu.star,decreasing = T)
  
  reg_coefficient <- matrix(mu.star)
  totalCoef <- sum(reg_coefficient)
  
  rankedVars <- vars[ranks]
  nLevelVar <- floor(nVar/nLevel)
  for(groupingIndex in 1:(nLevel-1)){
    start <- (groupingIndex-1)*nLevelVar+1
    end <- groupingIndex*nLevelVar
    unorderedCurrentGroup <- rankedVars[start:end]
    orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
    group <- append(group,list(orderedGroup))
    groupWeight <- append(groupWeight,sum(reg_coefficient[orderedGroup]))
  }
  start <- (nLevel-1)*nLevelVar+1
  end <- nVar
  unorderedCurrentGroup <- rankedVars[start:end]
  orderedGroup <- unorderedCurrentGroup[order(unorderedCurrentGroup)]
  groupWeight <- append(groupWeight,sum(reg_coefficient[orderedGroup]))
  group <- append(group,list(orderedGroup))
  print('grouping done')
  CV <- x
  CVobjs <- objs
  
  base <- F
  if(testIDinUse==1)
    base <- T
  
  for(i in 1:1){
    set.seed(i*100)
    start_1 <- Sys.time()
    gs <- c("",rep(c("roundrobin","roundrobin","roundrobin",
                     "adaptive","adaptive","adaptive"),3))
    
    gs <- gs[testIDinUse]
    
    # group <- list((1:20),(21:40),(41:60),(61:80),(81:100))
    
    population <- NULL
    for(groupingIndex in 1:length(group)){
      thisnVar <- length(group[[groupingIndex]])
      member <- group[[groupingIndex]]
      population_group <- t(InitializePopulationSobol(ncol = 98,
                                                    nrow = thisnVar,
                                                    minVal = lbound[member],
                                                    maxVal = ubound[member]))
      # population_group <- rbind(population_group,(lbound[]))
      newpop <- lbound[member]
      newpop[1:10] <- ubound[member[1:10]]
      population_group <- rbind(population_group,newpop)

      
      newpop <- lbound[member]
      newpop[11:20] <- ubound[member[11:20]]
      population_group <- rbind(population_group,newpop)
      
      population <- cbind(population,population_group)
    }
    for(i in 1:100)population[i,21:100] <- (21:100)*2*0.35
    # group <- list((1:20))
    
    # population <- MaOEA::InitializePopulationLHS(100,nVar,lbound,ubound)
    # population <- matrix(runif(100*nVar),nrow=100)*2*(1:100)
    # for(i in 1:100)population[21:100,i] <- (21:100)*2*0.35
    # population = t(population)
    # 
    # population <- NULL
    
    ctrl <- list(group = group,importance=reg_coefficient,
                 lower = lbound,upper = ubound,
                 budget=budget-nEval,refmultiplier=2,
                 stepadapt=T,
                 crossoverrate=0.5,
                 mutationrate=0.05,
                 stepsize=0.3,
                 alpha=1.2,
                 crossoveroperator="uniform",
                 mutationoperator="truncgauss",
                 localmutation=T,
                 vectorized=T, scaleinput=F,loghv=F,refpoint=c(2.2,4.4))
    if(!base){
      abc <- CC_SMSEMOA( population = population, nVar = nVar, populationSize  = 100,
                         fun = transposedFun, mainfunc=func, k=k, nObj =2,growPop=F,
                         maxPopulationSize=100,nIterInner=nIterInner,
                         control=ctrl,group_selection=gs)
    }else{
      population <- MaOEA::InitializePopulationLHS(100,nVar,lbound,ubound)
      for(i in 1:budget){
        print( paste("progress:",100*i/budget) )
        start_2 <- Sys.time()
        if(i==budget){
          final <- MaOEA::SMSEMOA(population,fun = func,nObjective = 2,
                                  nObj =2, k=k,
                                  control = list(lower=lbound,
                                                 upper=ubound,
                                                 ref_multiplier=2))
          abc <- list()
          abc$x <- t(final$population)
          abc$y <- t(final$objective)
        }else{
          population <- MaOEA::SMSEMOA(population,fun = func,nObjective = 2,
                                       nObj=2, k=k,
                                       control = list(lower=lbound,
                                                      upper=ubound,
                                                      ref_multiplier=2))$population
        }
        print(difftime(Sys.time(),start_2,units="secs"))
      }
    }
    end_1 <- Sys.time()
    #abc_log <- append(abc_log,abc)
    timelog <- append(timelog,as.numeric(end_1-start_1))
    load("adapt_step.Rdata")
    adapt_step_log <- append(adapt_step_log,list(cdlog=cdlog,
                                                 mdlog=mdlog,
                                                 sec_steplog=altstep_log,
                                                 stepsizelog=stepsize_log))
    save(list=ls(),file=paste0("checkgauss",i,".Rdata"))
  }
}

# nPoint <- 1000
# # paretoSetSample <- matrix((randtoolbox::sobol(n = nPoint,dim = k,scrambling = 3))*(1:k)*2,nrow=k)
# # paretoSetSample <- t(matrix((randtoolbox::sobol(k*nPoint))*(1:k)*2,ncol=k))
# # paretoSetSample <- matrix((runif(k*nPoint)^50)*(1:k)*2,nrow=k)
# #paretoSetSample <- matrix((runif(k*nPoint))*(1:k)*2,nrow=k)
# paretoSetSample <- matrix(1:nPoint * (1/nPoint),nrow=1)*(1:k)*2 # wfg2
# aa <- pracma::repmat(rep(0,nPoint),k-1,1)
# paretoSetSample <- rbind(paretoSetSample,aa)
# bb <- pracma::repmat(matrix(0.7*((k+1):nVar),ncol=1),1,nPoint)
# paretoSetSample <- rbind(paretoSetSample,bb)
# paretoFrontSample <- func(individual = paretoSetSample,nObj = 2,k = 20)
# xx <- (0:1000)/1000
# yy <- (1-xx*xx)^0.5
# xx <- paretoFrontSample[1,order(paretoFrontSample[1,])]
# yy <- paretoFrontSample[2,order(paretoFrontSample[1,])]
# plot(xx,yy,type="l",xlim= c(0,2),ylim=c(0,2))
# points(t(t(abc$y)/c(2,4)),col='red')
