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
#' groups <- RDG3(nVar=100,f1cec,control,o=rep(13,100))
#' @references <doi:10.1109/TEVC.2017.2694221>
#' @export
RDG3 <- function(fun,lbound,ubound,max_group_size=50,...){
  nVar <- (length(lbound))
  if(nVar==1){
    separable <- 1
    group <- vector(mode='list')
    nEval <- 0
    DSM <- as.matrix(1)
  }else{
    x <- NULL
    y <- NULL

    if(any(is.infinite(lbound))){
      x <- rep(0,nVar)
    }else{
      x <- t(matrix(lbound))
    }
    fun_repeat <- fun(x,...)
    y <- fun_repeat
    nEval <- 1
    group <- vector('list')
    seps <- NULL
    vars <- 1:nVar

    thisvar_group <- vars[1]
    vars <- vars[-1]
    while(length(vars)>0){
      #print('main length')
      #print(paste(length(vars),thisvar_group))

      cur_length <- length(thisvar_group)
      interact <- check_interaction(fun,thisvar_group,vars,lbound,ubound,fun_repeat,...)
      thisvar_group <- interact$interactingVar
      #print('NCI')
      #print(thisvar_group)
      new_length <- length(thisvar_group)

      if(new_length==cur_length || new_length>=max_group_size){
        if(length(thisvar_group)>1){
          group <- append(group,list(thisvar_group))
          nEval <- nEval + interact$nEval
          x <- rbind(x,interact$x)
          y <- append(y,interact$y)
        }else{
          seps <- append(seps,thisvar_group)
          nEval <- nEval + interact$nEval
          x <- rbind(x,interact$x)
          y <- append(y,interact$y)
        }
        vars <- setdiff(vars,thisvar_group)
        thisvar_group <- vars[1]
      }
      # use setdiff
      vars <- setdiff(vars,thisvar_group)
    }

  }
  return(list(group=group,separable=seps,nEval = nEval,x=x,y=y))
}

check_interaction <- function(fun,first_group,second_group,lbound,ubound,f_lbound,...){
  first_group_upper <- matrix(lbound)
  first_group_upper[first_group] <- (lbound[first_group]+ubound[first_group])/2
  f1_upper <- fun(drop(first_group_upper),...)

  second_group_upper <- matrix(lbound)
  second_group_upper[second_group] <- (lbound[second_group] +ubound[second_group] )/2
  f2_upper <- fun(drop(second_group_upper),...)

  both_group_upper <- matrix(lbound)
  both_group_upper[c(first_group,second_group)] <- (lbound[c(first_group,second_group)]+
                                                      ubound[c(first_group,second_group)])/2
  both_upper <- fun(drop(both_group_upper),...)
  nEval <- 3

  x <- cbind(first_group_upper,second_group_upper)
  x <- cbind(x,both_group_upper)
  x <- t(x)

  y <- append(f1_upper,f2_upper)
  y <- append(y, both_upper)

  delta1 <- f1_upper - f_lbound
  delta2 <- both_upper - f2_upper

  difference <- abs(delta2-delta1)
  maxfun <- max(c(f_lbound,f1_upper,f2_upper,both_upper))
  maxsum <- max(c(f1_upper+f2_upper,both_upper+f_lbound))

  summedfun <- abs(f_lbound)+abs(f1_upper)+abs(f2_upper)+abs(both_upper)
  eInf1 <- gammaFunc(nVar^0.5+2) * summedfun

  eInf <- gammaFunc(2) * maxsum
  eSup <- gammaFunc(nVar^0.5) * maxfun

  tolerance <- eInf1 # original RDG3
  # tolerance <- (eInf+eSup)/2
  interactingVar <- first_group
  if(difference>tolerance){
    second_group_length <- length(second_group)
    if(second_group_length==1){
      interactingVar <- union(first_group,second_group)

    }else{
      G1 <- second_group[1:round(second_group_length/2)]
      G2 <- second_group[-(1:round(second_group_length/2))]
      X1_1 <- check_interaction(fun,first_group,G1,lbound,ubound,f_lbound,...)
      X1_2 <- check_interaction(fun,first_group,G2,lbound,ubound,f_lbound,...)
      interactingVar <- union(X1_1$interactingVar,X1_2$interactingVar)
      nEval <- nEval+X1_2$nEval+X1_1$nEval
      x <- rbind(x,X1_1$x)
      x <- rbind(x,X1_2$x)

      y <- append(y,X1_1$y)
      y <- append(y,X1_2$y)
    }
  }#else{
   #  print(paste(difference,tolerance))
   #}
  return(list(interactingVar=interactingVar,nEval=nEval,x=x,y=y))
}

