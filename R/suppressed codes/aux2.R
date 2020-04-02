library(ccoevolution2)
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
func <- f5cec
res <- NULL
resVal <- NULL
for(qqIndex in 1:length(qq_mat)){
  qq <- qq_mat[qqIndex]
  methods <- 3
  seed <- ((qq-1)%/%8)*1000
  groupStart <- 1
  group <- NULL
  groupEnd <- 0
  permutation <- NULL
  set.seed(seed)
  rotation_matrix <- matrix(numeric(1000*1000),nrow=1000)
  weight <- ccoevolution:::cec_weight(7)
  weight <- append(weight,1)
  for(i in 1:8){
    group <- append(group,list(groupEnd+(1:group_size[i])))
    permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
    groupEnd <- groupEnd + group_size[i]
  }
  for(i in 1:7){
    rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
  }
  rotation_matrix[group[[8]],group[[8]]] <- diag(1,group_size[8])
  print(paste0('algo',methods))
  if(methods==1){
    print('CCDG')
    a <- cc(nVar = nVar,fun=func,budget=budget,grouping_control=ctrl,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)
  }
  if(methods==2){
    print('cc_2')
    a <- cc_2(nVar = nVar,fun=func,budget=budget,grouping_control=ctrl,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)
  }
  if(methods==3){
    print('SACC_mustar')
    a <- SACC(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight,SA_method='morris_mu')
  }
  if(methods==4){
    print('TSCC_mustar')
    a <- TSCC(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight,SA_method='morris_mu')
  }
  if(methods==5){
    print('SACC_sobol')
    a <- SACC(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight,SA_method='sobol')
  }
  if(methods==6){
    print('TSCC_sobol')
    a <- TSCC(nVar = nVar,fun=func,budget=budget,nLevel=4,lbound=lbound,ubound=ubound,o=optimum,rotation_matrix=rotation_matrix,permutation=permutation,weight=weight,SA_method='sobol')
  }
  res <- cbind(res,a)
  resVal <- cbind(resVal,a$y)
  save(list=ls(),file=paste('algo',methods,'_fun',func_index+3,'.Rdata',sep=''))
}
