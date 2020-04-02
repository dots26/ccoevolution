rm(list=ls())
qq_mat <- 1:25
qq_mat <- (qq_mat)*4

library('ccoevolution')
nVar <- 100
optimum <- rep(13,nVar)
func_set <- c(f4cec_s,f5cec_s,f6cec_s,f7cec_s)
nfunc <- length(func_set)
ctrl <- list(lbound=rep(-100,nVar),ubound=rep(100,nVar),delta=rep(20,nVar))
budget <- 500000
lbound <- rep(-100,nVar)
ubound <- rep(100,nVar)
delta <- rep(20,nVar)

group_size <- c(25,10,10,55)
for(func_index in 1:nfunc){
  func <- func_set[[func_index]]
  res <- NULL
  resVal <- NULL
  for(qqIndex in 1:length(qq_mat)){
    qq <- qq_mat[qqIndex]
    methods <- ((qq-1)%%4)+1
    seed <- ((qq-1)%/%4)*1000
    set.seed(seed)

    rotation_matrix <- matrix(numeric(100*100),nrow=100)
    weight <- ccoevolution:::cec_weight(3)
    weight <- append(weight,1)
    group <- NULL
    groupEnd <- 0
    permutation <- NULL
    for(i in 1:4){
      group <- append(group,list(groupEnd+(1:group_size[i])))
      permutation <- append(permutation,sample.int(group_size[i])+groupEnd)
      groupEnd <- groupEnd + group_size[i]
    }
    for(i in 1:3){
      rotation_matrix[group[[i]],group[[i]]] <- soobench:::random_rotation_matrix(group_size[i])
    }
    rotation_matrix[group[[4]],group[[4]]] <- diag(1,group_size[4])


    print('CCDG')
    a <- cc_2(nVar = nVar,fun=func,budget=budget,grouping_control=ctrl,lbound=lbound,ubound=ubound,o=optimum,
              rotation_matrix=rotation_matrix,permutation=permutation,weight=weight)

    res <- cbind(res,a)
    resVal <- cbind(resVal,a$y)
    if(qqIndex==25)
      save(list=ls(),file=paste('CC2_fun',func_index+3,'.Rdata',sep=''))
  }
}
