#' Create groups of variable using extended differential grouping (XDG). Tolerance is set following DG2 scheme.
#'
#' @title Extended Differential Grouping
#' @param nVar Number of variables to be grouped
#' @param fun The objective function to be solved
#' @param control List of control parameter for DG
#' \code{lbound} A vector of lower bound values for each variable
#' \code{ubound} A vector of upper bound values for each variable
#' \code{delta} The shift from \code{lbound} to compute the gradient
#' @param ... Further arguments to be passed to fun
#' @examples
#' control <- list(lbound=rep(-100,100),ubound=rep(100,100),delta=0.1)
#' groups <- XDG(nVar=100,f1cec,control,o=rep(13,100))
#' @references <doi:10.1145/2739480.2754666>
#' @export
XDG <- function(nVar,fun,control=list(),...){
  con <- list(lbound=rep(0,nVar),
              ubound=rep(1,nVar),
              delta=NULL)

  con[names(control)] <- control
  center <- con$lbound + (con$ubound-con$lbound)/2

  nEval <- 0
  #doParallel::registerDoParallel(cores=4)
  group <- vector(mode='list')
  separable <- NULL
  ungrouped <- 1:nVar

  if((length(con$lbound) != nVar) || (length(con$ubound) != nVar))  stop('lbound and ubound must be a vector of size nVar')

  if(any(con$lbound > con$ubound)) stop('lbound > ubound')

  if(!is.null(con$delta)){ # if delta is specified, check if lbound + delta will cross ubound
    if(any((con$lbound+con$delta)>con$ubound)) stop('lbound+delta > ubound, consider smaller delta or removing it to use the midpoint.')
    center <- con$lbound+con$delta
  }

  p1 <- con$lbound
  fun1 <- fun(p1,...)


  p2 <- t(matrix(rep(p1,length(ungrouped)),ncol=length(ungrouped)))
  diag(p2) <- center
  fun_repeat <- fun(p2,...)

  groupID <- 1
  varGroup <- rep(0,length(ungrouped))
  for(i in 1:(length(ungrouped)-1)){
    this.group <- i
    a <- NULL

    p2 <- t(matrix(rep(p1,length(ungrouped)),ncol=length(ungrouped)))
    p2[,i] <- center[i]
    diag( p2) <- center
    funj <- fun(p2[(i+1):length(ungrouped),],...)

    for(j in (i+1):length(ungrouped)) {
      maxfun <- max(c(fun_repeat[i], fun1,funj[j-i],fun_repeat[j]))
      maxsum <- max(c(fun_repeat[i]+fun_repeat[j],fun1+funj[j-i]) )

      eInf <- gammaFunc(2) * maxsum
      eSup <- gammaFunc(nVar) * maxfun


      tolerance <- eSup

      delta1 <- fun_repeat[i] - fun1 #up lo lo - lo lo lo
      delta2 <- funj[j-i] - fun_repeat[j] #up up lo - lo up lo
      if(abs(delta2-delta1)>tolerance){
        additional_member <- j
      }else{
        additional_member <- NULL
      }
      a <- append(a,additional_member)
    }
    this.group <- append(this.group,unlist(a))
    previouslyUngrouped <- (max(varGroup[this.group]) == 0 )

    if(previouslyUngrouped){
      varGroup[this.group] <- groupID
      groupID <- groupID +1
    }else{
      varGroup[this.group] <- min(varGroup[this.group])
    }
  }


  for(i in 1:groupID){
    currentGroupLength <- length(which(varGroup==i))
    if(currentGroupLength>0)
      group <- append(group,list(which(varGroup==i)))
  }
  nEval <- 1+(nVar*(nVar+1)/2)
  return(list(group=group,nEval = nEval))
}

#' Create groups of variable using Differential Grouping. The code uses the supremum tolerance.
#' Infimum tolerance is calculated, but given zero weight.
#'
#' @title Differential Grouping
#' @param nVar Number of variables to be grouped
#' @param fun The objective function to be solved
#' @param control List of control parameter for DG
#' \code{lbound} A vector of lower bound values for each variable
#' \code{ubound} A vector of upper bound values for each variable
#' \code{delta} The shift from \code{lbound} to compute the gradient
#' @param ... Further arguments to be passed to fun
#' @examples
#' control <- list(lbound=rep(-100,100),ubound=rep(100,100),delta=0.1)
#' groups <- DG(nVar=100,f1cec,control,o=rep(13,100))
#' @references <doi:10.1109/TEVC.2017.2694221>

DG <- function(nVar,fun,control=list(),...){
  if(nVar==1){
    separable <- 1
    group <- vector(mode='list')
    nEval <- 0
  }else{
    con <- list(lbound=rep(0,nVar),
                ubound=rep(1,nVar),
                delta=NULL)

    con[names(control)] <- control
    center <- con$lbound + (con$ubound-con$lbound)/2

    nEval <- 0

    group <- vector(mode='list')
    separable <- NULL
    ungrouped <- 1:nVar

    if((length(con$lbound) != nVar) || (length(con$ubound) != nVar))  stop('lbound and ubound must be a vector of size nVar')

    if(any(con$lbound > con$ubound)) stop('lbound > ubound')

    if(!is.null(con$delta)){ # if delta is specified, check if lbound + delta will cross ubound
      if(any((con$lbound+con$delta)>con$ubound)) stop('lbound+delta > ubound, consider smaller delta or removing it to use the midpoint.')
      center <- con$lbound+con$delta
    }

    p1 <- con$lbound
    fun1 <- fun(p1,...)

    p2_a <- t(matrix(rep(p1,length(ungrouped)),ncol=length(ungrouped)))
    diag(p2_a) <- center
    fun_repeat <- fun(p2_a,...)

    nEval <- 1+nVar
    separable <- NULL
    while(length(ungrouped)>1){
      i <- ungrouped[1]
      this.group <- i
      a <- NULL

      p2 <- t(matrix(rep(p1,length(ungrouped)-1),ncol=length(ungrouped)-1))
      p2[,ungrouped[1]] <- center[ungrouped[1]]
      for(ix in 1:nrow(p2)){
        p2[ix,ungrouped[ix+1]] <- center[ungrouped[ix+1]]
      }

      funj <- fun(p2,...)
      nEval <- nEval + nrow(p2)

      for(j in (1:nrow(p2)) ) {
        maxfun <- max(c(fun_repeat[i], fun1,funj[j],fun_repeat[j]))
        maxsum <- max(c(fun_repeat[i]+fun_repeat[j],fun1+funj[j]) )

        eInf <- gammaFunc(2) * maxsum
        eSup <- gammaFunc(nVar) * maxfun
        tolerance <- 0.5*eInf + 0.5*eSup
        delta1 <- fun_repeat[i] - fun1 #up lo lo - lo lo lo
        delta2 <- funj[j] - fun_repeat[ungrouped[j+1]] #up up lo - lo up lo



        if(abs(delta2-delta1)>tolerance){
          additional_member <- ungrouped[j+1]
        }else{
          additional_member <- NULL
        }
        a <- append(a,additional_member)
      }
      this.group <- append(this.group,unlist(a))
      if(length(this.group)==1){
        separable <- append(separable,this.group)
      }else{
        group <- append(group,list(this.group))
      }

      ungrouped <- ungrouped[-(match(this.group,ungrouped))]
    }
    if(length(ungrouped)==1)
      separable <- append(separable,ungrouped)
    # group <- append(group,list(separable))
  }
  return(list(group=group,separable=separable,nEval = nEval))
}

#' Create groups of variable using Differential Grouping. The code uses the supremum tolerance.
#' Infimum tolerance is calculated, but given zero weight.
#'
#' @title Differential Grouping version 2
#' @param nVar Number of variables to be grouped
#' @param fun The objective function to be solved
#' @param control List of control parameter for DG
#' \code{lbound} A vector of lower bound values for each variable
#' \code{ubound} A vector of upper bound values for each variable
#' \code{delta} The shift from \code{lbound} to compute the gradient
#' @param ... Further arguments to be passed to fun
#' @return A list
#'         \code{group} Lists of the groups found by DG2
#'         \code{separable} Vector of the separable variables found by DG2
#'         \code{DSM} The full design structure matrix indicating the interaction graph.
#'         \code{nEval} Number of evaluation used in grouping
#'         \code{x} Points evaluated
#'         \code{y} Objective value wrt x, i.e. f(x). If multiobjective, each column is an objective, each row is for an individual
#' @examples
#' control <- list(lbound=rep(-100,100),ubound=rep(100,100),delta=0.1)
#' groups <- DG2(nVar=100,f1cec,control,o=rep(13,100))
#' @references <doi:10.1109/TEVC.2017.2694221>
#' @export
DG2 <- function(nVar,fun,control=list(),...){
  if(nVar==1){
    separable <- 1
    group <- vector(mode='list')
    nEval <- 0
    DSM <- as.matrix(1)
  }else{
    x <- NULL
    y <- NULL
    con <- list(lbound=rep(0,nVar),
                ubound=rep(1,nVar),
                delta=NULL)

    con[names(control)] <- control
    center <- con$lbound + (con$ubound-con$lbound)/2

    nEval <- 0

    group <- vector(mode='list')
    separable <- NULL
    ungrouped <- 1:nVar

    if((length(con$lbound) != nVar) || (length(con$ubound) != nVar))  stop('lbound and ubound must be a vector of size nVar')

    if(any(con$lbound > con$ubound)) stop('lbound > ubound')

    if(!is.null(con$delta)){ # if delta is specified, check if lbound + delta will cross ubound
      if(any((con$lbound+con$delta)>con$ubound)) stop('lbound+delta > ubound, consider smaller delta or removing it to use the midpoint.')
      center <- con$lbound+con$delta
    }

    DSM <- pracma::zeros(nVar)*NA
    diag(DSM) <- 1
    ISM <- pracma::zeros(nVar)
    p1 <- con$lbound
    fun1 <- fun(p1,...)

    x <- rbind(x,p1)
    y <- rbind(y,fun1)
    p2_a <- t(matrix(rep(p1,nVar),ncol=nVar))
    diag(p2_a) <- center
    # vectorized, not working on ubuntu
    # fun_repeat <- fun(p2_a,...)
    x <- rbind(x,p2_a)

    fun_repeat <- vector()
    for(i in 1:nrow(p2_a))
      fun_repeat[i] <- fun(p2_a[i,],...)

    y <- rbind(y,fun_repeat)


    nEval <- 1+nVar
    funj <- list()
    for(i in 1:(nVar-1)){
      p2 <- t(matrix(rep(p1,nVar-i),ncol=nVar-i))
      p2[,i] <- center[i]
      for(ix in 1:nrow(p2)){
        p2[ix,ix+i] <- center[ix+i]
      }
      x <- rbind(x,p2)
      funj_this <- vector()
      for(i_this in 1:nrow(p2))
        funj_this[i_this] <- fun(p2[i_this,],...)
      y <- rbind(y,funj_this)

      # vectorized, not working on ubuntu
      # funj <- append(funj,list(fun(p2,...)))
      funj <- append(funj,list(funj_this))
      nEval <- nEval + nrow(p2)

      for(j in (1:nrow(p2)) ) {
       # print(c(fun_repeat[i],fun1,funj[[i]][j],fun_repeat[j+1]))
        delta1 <- fun_repeat[i] - fun1 #up lo lo - lo lo lo
        delta2 <- funj[[i]][j] - fun_repeat[j+i] #up up lo - lo up lo

        ISM[i,j+i] <- abs(delta2-delta1)
      }
    }

    eta0 <- 0
    eta1 <- 0
    for(i in 1:(nVar-1)){
      for(j in (i+1):nVar){
        maxfun <- max(c(fun_repeat[i], fun1,funj[[i]][j-i],fun_repeat[j]))
        maxsum <- max(c(fun_repeat[i]+fun_repeat[j],fun1+funj[[i]][j-i]) )
        eInf <- gammaFunc(2) * maxsum
        eSup <- gammaFunc(nVar^0.5) * maxfun

        # if(!is.nan(ISM[i,j])){
        # print(c(ISM[i,j],eInf,eSup))
          if(ISM[i,j]<eInf){
            DSM[i,j] <- 0
            DSM[j,i] <- 0
            eta0 <- eta0+1
          }else if(ISM[i,j]>eSup){
            DSM[i,j] <- 1
            DSM[j,i] <- 1
            eta1 <- eta1+1
          }
      }
    }


    for(i in 1:(nVar-1)){
      for(j in (i+1):nVar){
        maxfun <- max(c(fun_repeat[i], fun1,funj[[i]][j-i],fun_repeat[j]))
        maxsum <- max(c(fun_repeat[i]+fun_repeat[j],fun1+funj[[i]][j-i]) )

        eInf <- gammaFunc(2) * maxsum
        eSup <- gammaFunc(nVar^0.5) * maxfun

        # print(DSM[i,j])
        if(is.na(DSM[i,j])){
          if((eta0+eta1)>0){
            tolerance <- eta0/(eta0+eta1)*eInf + eta1/(eta0+eta1)*eSup
          }else{
            tolerance <- eSup
          }
          # print(c(ISM[i,j],tolerance))
          if(ISM[i,j]>tolerance){
            # print('suc')
            DSM[i,j] <- 1
            DSM[j,i] <- 1
          }else{
            # print('fail')
            DSM[i,j] <- 0
            DSM[j,i] <- 0
          }
        }
      }
    }
# save(DSM,ISM,eta0,eta1,fun_repeat,fun1,funj,file='DSM.Rdata')
    separable <- NULL
    group <- NULL
    groupID <- 1
    assignedGroup <- integer(nVar)

    for(i in 1:(nVar)){
      if(sum(DSM[i,])==1){
        separable <- append(separable,i)
      }else{
        this.group <- which(DSM[i,]>0)
        for(j in this.group){
          if(max(assignedGroup[this.group]!=0)){
            assignedGroup[this.group] <- max(assignedGroup[this.group])
          }else{
            assignedGroup[this.group] <- groupID
            groupID <- groupID+1
          }
        }
      }
    }
    for(i in 1:groupID){
      groupI <- which(assignedGroup==i)
      if(length(groupI)>0){
        group <- append(group,list(groupI))
      }
    }
  }
  return(list(group=group,separable=separable,DSM=DSM,nEval = nEval,x=x,y=y))
}



gammaFunc <- function(d){
  machinePrecision <- .Machine$double.eps/2
  return((d * machinePrecision)/(1 - (d * machinePrecision)))
}
