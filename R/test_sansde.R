
# Optimize the subcomponent using SaNSDE
# The SaNSDE algorithm can be found in:
# Zhenyu Yang, Ke Tang and Xin Yao, "Self-adaptive Differential Evolution with
# Neightborhood Search", in Proceedings of the 2008 IEEE Congress on
# Evolutionary Computation (CEC2008), Hongkong, China, 2008, pp. 1110-1116.

## need to consider bestmem as parent?
sansde <- function(fname, pop, bestmem=NULL, bestval=Inf, Lbound, Ubound,control =list(),...){
  # each row is an individual
  con <- list(ccm=0.5,
              itermax=100)
  con[names(control)] <- control
  control <- con
  ccm <- control$ccm
  itermax <- control$itermax

  popsize <- dim(pop)[1]
  D <- dim(pop)[2]
  NP <- popsize
  tracerst <- NULL

  F_mat <- pracma::zeros(NP,1)

  linkp <- 0.5
  l1 <- 1
  l2 <- 1
  nl1 <- 1
  nl2 <- 1

  fp <- 0.5
  ns1 <- 1
  nf1 <- 1
  ns2 <- 1
  nf2 <- 1

  pm1 <- pracma::zeros(NP,D)# initialize population matrix 1
  pm2 <- pracma::zeros(NP,D)# initialize population matrix 2
  pm3 <- pracma::zeros(NP,D)# initialize population matrix 3
  pm4 <- pracma::zeros(NP,D)# initialize population matrix 4
  pm5 <- pracma::zeros(NP,D)# initialize population matrix 5
  bm  <- pracma::zeros(NP,D)# initialize DE_gbestber  matrix
  ui  <- pracma::zeros(NP,D)# intermediate population of perturbed vectors
  mui <- pracma::zeros(NP,D)# mask for intermediate population
  mpo <- pracma::zeros(NP,D)# mask for old population
  rot <- (0:(NP-1))# rotating index array (size NP)
  rotd <-  (0:(D-1))# rotating index array (size D)
  rt  <- pracma::zeros(NP)# another rotating index array
  rtd <- pracma::zeros(D)# rotating index array for exponential crossover
  a1  <- pracma::zeros(NP)# index array
  a2  <- pracma::zeros(NP)# index array
  a3  <- pracma::zeros(NP)# index array
  a4  <- pracma::zeros(NP)# index array
  a5  <- pracma::zeros(NP)# index array
  ind <- pracma::zeros(4)

  cc_rec <- NULL
  f_rec <- NULL

  gpop <- pop

  val <- fname(gpop,...)
  # val <- fname(oldpop)
  # val <- fname(gpop,rotation_matrix = rotation_matrix,permutation = permutation,o=rep(0,1000),weight=weight)

  ibest <- which.min(val)
  best <- min(val)
  subbestmem <- pop[ibest,]
  if(is.null(bestmem)){
    bestmem <- subbestmem
  }
  used_FEs <- 0
  # if (gcount > 1){
  #   used_FEs <- NP
  # }

  if (best < bestval){
    bestval <- best
    bestmem <- gpop[ibest, ]
  }

  iter <- 0
  while (iter < itermax){
    popold <- pop# save the old population

    ind <- pracma::randperm(4)# index pointer array

    a1  <- pracma::randperm(NP)# shuffle locations of vectors
    rt <- pracma::rem(rot+ind[1],NP)# rotate indices by ind(1) positions
    a2  <- a1[rt+1]# rotate vector locations
    rt <- pracma::rem(rot+ind[2],NP)
    a3  <- a2[rt+1]
    rt <- pracma::rem(rot+ind[3],NP)
    a4  <- a3[rt+1]
    rt <- pracma::rem(rot+ind[4],NP)
    a5  <- a4[rt+1]

    pm1 <- popold[a1,,drop=F]# shuffled population 1
    pm2 <- popold[a2,,drop=F]# shuffled population 2
    pm3 <- popold[a3,,drop=F]# shuffled population 3
    pm4 <- popold[a4,,drop=F]# shuffled population 4
    pm5 <- popold[a5,]# shuffled population 5

    bm <- pracma::ones(NP, 1) %*% bestmem #changed from subbestmem

    if (pracma::rem(iter,24)==0){
      if ((iter!=0) && (!is.null(cc_rec))){
        ccm <- sum(f_rec*cc_rec)/sum(f_rec)
      }
      cc_rec <- NULL
      f_rec <- NULL
    }
    if (pracma::rem(iter,5)==0){
      cc_sansde <- matrix(rnorm(NP*3,ccm, 0.1))
      index <- which((cc_sansde < 1) & (cc_sansde > 0))
      cc_sansde <- cc_sansde[index[1:NP],,drop=F]
    }

    fst1 <- (runif(NP) <= fp)
    fst2 <- 1-fst1

    fst1_index <- which(fst1 != 0)
    fst2_index <- which(fst1 == 0)

    tmp <- matrix(rnorm(NP,0.5, 0.3))
    F_mat[fst1_index] <- tmp[fst1_index]

    tmp <- matrix(rnorm(NP,0, 1)) / matrix(rnorm(NP,0, 1)) # cauchy
    F_mat[fst2_index] <- tmp[fst2_index]

    F_mat <- matrix(abs(F_mat))

    # all random numbers < CR are 1, 0 otherwise

    aa <- matrix(runif(NP*D),nrow=NP) < pracma::repmat(cc_sansde,1,D)
    index <- which(rowSums((aa)) == 0)
    tmpsize <- length(index)

    if(tmpsize>0){
      for (k in 1:tmpsize){
        bb <- ceiling(D*runif(1))
        aa[index[k], bb] <- 1
      }
    }

    mui <- aa
    mpo <- mui < 0.5# inverse mask to mui

    aaa <- (runif(NP) <= linkp)
    aindex=which(aaa == 0)
    bindex=which(aaa != 0)

    save(pm1,pm2,pm3,pm4,file='pm1.Rdata')
    if (!pracma::isempty(bindex)){
      # mutation
      ui[bindex,] <- popold[bindex,]+ pracma::repmat(F_mat[bindex,,drop=F],1,D) * (bm[bindex,,drop=F]-popold[bindex,,drop=F]) + pracma::repmat(F_mat[bindex,,drop=F],1,D) * (pm1[bindex,,drop=F] - pm2[bindex,,drop=F] + pm3[bindex,,drop=F] - pm4[bindex,,drop=F])
      # ui[bindex,] <- popold[bindex,] + pracma::repmat(F_mat[bindex,,drop=F],1,D) * (pm1[bindex,] - pm2[bindex,] + pm3[bindex,] - pm4[bindex,])
      # crossover
      ui[bindex,] <- popold[bindex,,drop=F]*mpo[bindex,,drop=F] + ui[bindex,,drop=F]*mui[bindex,,drop=F]
    }
    if (!pracma::isempty(aindex)){
      ui[aindex,] <- pm3[aindex,,drop=F] + pracma::repmat(F_mat[aindex,,drop=F],1,D)*(pm1[aindex,,drop=F] - pm2[aindex,,drop=F])
      ui[aindex,] <- popold[aindex,,drop=F]*mpo[aindex,,drop=F] + ui[aindex,,drop=F]*mui[aindex,,drop=F]
    }

    bbb <- 1-aaa

    # -----Select which vectors are allowed to enter the new population-------
    for(varIndex in 1:D){
      index <- which(ui[,varIndex] > Ubound[varIndex])
      ui[index,varIndex] <- Ubound[varIndex] - pracma::mod((ui[index,varIndex]-Ubound[varIndex]),(Ubound[varIndex]-Lbound[varIndex]))
      index <- which(ui[,varIndex] < Lbound[varIndex])
      ui[index,varIndex] <- Lbound[varIndex] + pracma::mod((Lbound[varIndex]-ui[index,varIndex]),(Ubound[varIndex]-Lbound[varIndex]))
    }

    gpop <- ui
    tempval <- fname(gpop,...)
    # tempval <- fname(gpop)
    # tempval <- fname(gpop,rotation_matrix = rotation_matrix,permutation = permutation,o=rep(0,1000),weight=weight)

    used_FEs <- used_FEs + NP

    for (i in 1:NP){
      if (tempval[i] <= val[i]){
        if (tempval[i] < val[i]){
          cc_rec <- cbind(cc_rec, cc_sansde[i,1])
          f_rec <- cbind(f_rec,(val[i] - tempval[i]))
        }

        pop[i,] <- ui[i,]
        val[i]   <- tempval[i]

        l1 <- l1 + aaa[i]
        l2 <- l2 + bbb[i]

        ns1 <- ns1 + fst1[i]
        ns2 <- ns2 + fst2[i]
      } else {
        nl1 <- nl1 + aaa[i]
        nl2 <- nl2 + bbb[i]

        nf1 <- nf1 + fst1[i]
        nf2 <- nf2 + fst2[i]
      }
    }

    if ((pracma::rem(iter,24) == 0) && (iter!=0)){
      linkp <- (l1/(l1+nl1))/(l1/(l1+nl1)+l2/(l2+nl2))
      l1 <- 1
      l2 <- 1
      nl1 <- 1
      nl2 <- 1
      fp <- (ns1 * (ns2 + nf2))/(ns2 * (ns1 + nf1) + ns1 * (ns2 + nf2))
      ns1 <- 1
      nf1 <- 1
      ns2 <- 1
      nf2 <- 1
    }

    # val <- fname (gpop,...)
    ibest <- which.min(val)
    bestval <- min(val)
    subbestmem <- pop[ibest,]

    if (best < bestval){
      bestval <- best
      bestmem <- pop[ibest, ]
    }

    #print(paste(iter,bestval))

    tracerst <- rbind(tracerst, bestval)
    iter <- iter + 1
  }

  popnew <- pop
  bestmemnew <- bestmem
  bestvalnew <- bestval

  return (list(popnew= popnew,
               par=bestmemnew,
               value=bestvalnew,
               tracerst=tracerst,
               ccm=ccm,
               used_FEs=used_FEs))
}



