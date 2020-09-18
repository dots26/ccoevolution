rm(list=ls())
qq <- 1#as.integer(Sys.getenv("PBS_ARRAYID"))
testID <- ((qq-1)%/%10)+1 # 1-10, 11-20
instance <- (qq-1)%%10+1

library(MaOEA)
# setwd("~/Desktop/MOLSGO/ccoevolution")
files.sources = list.files(path = "R",full.names = T,pattern = ".R")
sapply(files.sources, source)
source('WOF_SMSEMOA.R')
set.seed(100)
func <- MaOEA::WFG2
adapt_step_log <- list()
# formals(func) <- alist(individual=nObj=m,k= nObj-1)
m <- 2
evalInterval <- 1e5
vargroup <- vector(mode = 'list',length = m)
a <- vector(mode = 'list',length = m)
r <- 1
nVar <- 100
lbound <- rep(0,nVar)
ubound <- (1:nVar)*2
nEval <- 0
k <- 20
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

# for (testID in c(4))
{
  print(paste("TestID",testID))
  if(testID>13){
    testIDinUse <- testID-13
  }else{
    testIDinUse <- testID
  }
  nLevel <- c(1,50,10,4,50,10,4,50,10,4,50,10,4)
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
  rnk <- nsga2R::fastNonDominatedSorting(objs)
  #print('non_domsort done')
  rnkIndex <- integer(ncol(objs))
  sortingIndex <- 1;
  while (sortingIndex <= length(rnk)) {
    rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
    sortingIndex <- sortingIndex + 1;
  }
  print('sorting done')
  # fill new population
  
  print('population built')
  
  reg_coefficient <- vector('list')
  group <- NULL
  groupWeight <- NULL
  mu.star <- apply(SA$ee, 2, function(x) mean(abs(x)))
  
  ranks <- order(mu.star,decreasing = T)
  
  reg_coefficient <- matrix(mu.star)
  # for(i in 1:m){
  totalCoef <- sum(reg_coefficient)
  
  group <- NULL
  group <- append(group,list(1:20))
  group <- append(group,list(21:40))
  group <- append(group,list(41:60))
  group <- append(group,list(61:80))
  group <- append(group,list(81:100))
  print('grouping done')
  CV <- x
  CVobjs <- objs
  
  ctrl <- list(group = group,importance=reg_coefficient,
               lower = lbound,upper = ubound,
               budget=budget-nEval,refmultiplier=4,
               stepadapt=F,
               # refpoint=c(2.2,4.4),
               #refpointmethod="static",
               #budget = budget,
               vectorized=T, scaleinput=F,loghv=F)
  
  base <- F
  for(i in 1:10){
    set.seed(i*100)
    start_1 <- Sys.time()
    gs <- c("",rep(c("roundrobin","roundrobin","roundrobin",
                     "adaptive","adaptive","adaptive"),2))
    
    gs <- gs[testIDinUse]
    abc <- WOF_SMSEMOA( population = NULL, nVar = nVar, populationSize  = 80,
                        fun = func, k=k, nObj =2,growPop=T,
                        maxPopulationSize=80,nIterInner=10,
                        control=ctrl,group_selection=gs)
    
    end_1 <- Sys.time()
    abc_log <- append(abc_log,abc)
    timelog <- append(timelog,as.numeric(end_1-start_1))
    save(list=ls(),file=paste0('results/niter_sensitivity/WFG2_100iter',i,'.Rdata'))
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
# xx <- paretoFrontSample[1,order(paretoFrontSample[1,])]/2
# yy <- paretoFrontSample[2,order(paretoFrontSample[1,])]/4
# plot(xx,yy,type="l",xlim= c(0,2),ylim=c(0,2))
# points(t(t(abc$y)/c(2,4)),col='red')
