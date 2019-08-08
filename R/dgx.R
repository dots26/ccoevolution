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
#' control <- list(lbound=-100,ubound=100,delta=0.1)
#' groups <- differential_grouping_ex(100,f1cec,control)
#' @references Omidvar, et al.
#' @export
differential_grouping_ex <- function(nVar,fun,control=list(),...){
  con <- list(lbound=rep(0,nVar),
              ubound=rep(1,nVar),
              delta=NULL,
              tolerance=1e-3,
              maxLowerLevel = round(nVar/10))
  con[names(control)] <- control
  center <- con$lbound + (con$ubound-con$lbound)/2

  registerDoParallel(cores=4)
  print(getDoParWorkers())



  group <- vector(mode='list')
  separable <- NULL
  ungrouped <- 1:nVar

  bestPopIndex <- which.min(LASSOtarget)
  bestPop <- LASSOpop[,bestPopIndex]
  bestVal <- LASSOtarget[,bestPopIndex]

  if((length(con$lbound) != nVar) || (length(con$ubound) != nVar))  stop('lbound and ubound must be a vector of size nVar')

  if(any(con$lbound > con$ubound)) stop('lbound > ubound')

  if(!is.null(con$delta)){ # if delta is specified, check if lbound + delta will cross ubound
    if(any((con$lbound+con$delta)>con$ubound)) stop('lbound+delta > ubound, consider smaller delta or removing it to use the midpoint.')
    center <- con$lbound+con$delta
  }
  while(length(ungrouped)>1){
    #  print(ungrouped)
    varIndex <- ungrouped[1]
    this.group <- varIndex

    a <- foreach(secondaryVar = ungrouped[2:length(ungrouped)],.combine = 'cbind')%dopar%{
      p1 <- con$lbound
      p2 <- p1
      p2[varIndex] <- con$ubound[varIndex]
      delta1 <- fun(p1,...)- fun(p2,...)

      p1[secondaryVar] <- center[secondaryVar]
      p2[secondaryVar] <- center[secondaryVar]
      delta2 <- fun(p1,...)- fun(p2,...)

      if(abs(delta2-delta1)>con$tolerance){
        # this.group <- c(this.group,secondaryVar)
        additional_member <- secondaryVar
      }else{
        additional_member <- NULL
      }
      additional_member
    }
    this.group <- append(this.group,unlist(a))
    if(length(this.group==1)){ # the variable is separable
      separable <- append(separable,this.group)
    }else{
      group <- append(group,this.group)
    }
    #  print('avx')
    #  print(a)
    #  print(this.group)
    #  print(ungrouped)
    #  print(match(this.group,ungrouped))
    ungrouped <- ungrouped[-(match(this.group,ungrouped))] # remove the grouped variables from the ungrouped list
    #  print('new gro')
    print(ungrouped)
  }
  if(length(ungrouped)==1){
    separable <- append(separable,ungrouped)
  }
  group <- append(group,(separable)) # add the separable vars as the last group
}
