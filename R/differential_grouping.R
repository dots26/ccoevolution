#' Create groups of variable using differential grouping
#'
#' @title Differential Grouping
#' @param nVar Number of variables to be grouped
#' @param fun The objective function to be solved
#' @param lbound A vector of lower bound values for each variable
#' @param ubound A vector of upper bound values for each variable
#' @param delta The shift from \code{lbound} to compute the gradient
#' @param tolerance Tolerance \epsilon indicates the threshold on the objective function to determine whether two variables are separable or not. Smaller tolerance will make it easier for variables to be grouped together.
#' @param ... Further arguments to be passed to fun
#' @references Omidvar, et al.
#' @export
differential_grouping <- function(nVar,fun,lbound=rep(0,nVar),ubound=rep(1,nVar),delta=NULL,tolerance=1e-3,...){
  center <- lbound + (ubound-lbound)/2
  group <- vector(mode='list')
  separable <- NULL
  ungrouped <- 1:nVar

  if((length(lbound) != nVar) || (length(ubound) != nVar))  stop('lbound and ubound must be a vector of size nVar')

  if(any(lbound > ubound)) stop('lbound > ubound')

  if(!is.null(delta)){ # if delta is specified, check if lbound + delta will cross ubound
    if(any((lbound+delta)>ubound)) stop('lbound+delta > ubound, consider smaller delta or removing it to use the midpoint.')
    center <- lbound+delta
  }
  for(varIndex in ungrouped){
    this.group <- varIndex
    for(secondaryVar in (varIndex+1):nVar){
      p1 <- lbound
      p2 <- p1
      p2[varIndex] <- ubound
      delta1 <- fun(p1,...)- fun(p2,...)
      p1[secondaryVar] <- center[secondaryVar]
      p2[secondaryVar] <- center[secondaryVar]
      delta2 <- fun(p1,...)- fun(p2,...)

      if(abs(delta2-delta1)>tolerance){
        this.group <- c(this.group,secondaryVar)
      }
    }
    if(length(this.group==1)){ # the variable is separable
      separable <- append(separable,this.group)
    }else{
      group <- append(group,this.group)
    }
    ungrouped <- ungrouped[-(match(this.group,ungrouped))] # remove the grouped variables from the ungrouped list
  }
  group <- append(group,separable) # add the separable vars as the last group
}
