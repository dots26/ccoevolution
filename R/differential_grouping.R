#' Create groups of variable using differential grouping
#'
#' @title Differential Grouping
#' @param nVar Number of variables to be grouped
#' @param fun The objective function to be solved
#' @param control List of control parameter for DG
#' \code{lbound} A vector of lower bound values for each variable
#' \code{ubound} A vector of upper bound values for each variable
#' \code{delta} The shift from \code{lbound} to compute the gradient
#' \code{tolerance} Tolerance indicates the threshold on the objective function to determine whether two variables are separable or not. Smaller tolerance will make it easier for variables to be grouped together.
#' @param ... Further arguments to be passed to fun
#' @examples
#' control <- list(lbound=rep(-100,100),ubound=rep(100,100),delta=0.1)
#' groups <- differential_grouping(nVar=100,f1cec,control,o=rep(13,100))
#' @references Omidvar, et al.
#' @export
differential_grouping <- function(nVar,fun,control=list(),...){
  con <- list(lbound=rep(0,nVar),
              ubound=rep(1,nVar),
              delta=NULL,
              tolerance=1e-3)
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
  nEval <- nEval + 1

  p2 <- t(matrix(rep(p1,length(ungrouped)),ncol=length(ungrouped)))
  fun_repeat <- fun(p2,...)

  nEval <- nEval + length(ungrouped)
  groupID <- 1
  varGroup <- rep(0,length(ungrouped))
  for(i in 1:(length(ungrouped)-1)){
    this.group <- i
    a <- NULL

    p2 <- t(matrix(rep(p1,length(ungrouped)),ncol=length(ungrouped)))
    p2[,i] <- center[i]
    for(j in (i+1):length(ungrouped)) {
      p2[,j] <- center[j]
    }
    funj <- fun(p2[(i+1):length(ungrouped),],...)

    for(j in (i+1):length(ungrouped)) {
      delta1 <- fun_repeat[i] - fun1
      delta2 <- funj[j-i] - fun_repeat[j]

      if(abs(delta2-delta1)>con$tolerance){
        additional_member <- j
      }else{
        additional_member <- NULL
      }
      a <- append(a,additional_member)
    }
    nEval <- nEval + length(ungrouped)*(length(ungrouped)-1)
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

  return(list(group=group,nEval = nEval))
}


differential_grouping_old <- function(nVar,fun,control=list(),...){
  con <- list(lbound=rep(0,nVar),
              ubound=rep(1,nVar),
              delta=NULL,
              tolerance=1e-3)
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
  while(length(ungrouped)>1){
    varIndex <- ungrouped[1]
    this.group <- varIndex

    a <- foreach::foreach(secondaryVar = ungrouped[2:length(ungrouped)],.combine = 'cbind') %dopar% {
      p1 <- con$lbound
      p2 <- con$lbound
      p2[varIndex] <- con$ubound[varIndex]

      delta1 <- fun(p1,...)- fun(p2,...)
      nEval <- nEval + 2
      p1[secondaryVar] <- center[secondaryVar]
      p2[secondaryVar] <- center[secondaryVar]
      delta2 <- fun(p1,...)- fun(p2,...)
      nEval <- nEval + 2

      if(abs(delta2-delta1)>con$tolerance){
        additional_member <- secondaryVar
      }else{
        additional_member <- NULL
      }
      additional_member
    }
    this.group <- append(this.group,unlist(a))
    if(length(this.group)==1){ # the variable is separable
      separable <- append(separable,this.group)
    }else{
      group <- append(group,list(this.group))
    }
    ungrouped <- ungrouped[-(match(this.group,ungrouped))] # remove the grouped variables from the ungrouped list
  }
  if(length(ungrouped)==1){
    separable <- append(separable,ungrouped)
  }
  group <- append(group,(separable)) # add the separable vars as the last group
  return(list(group=group,nEval = nEval))
}
