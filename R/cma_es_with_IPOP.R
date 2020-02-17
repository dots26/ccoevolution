##
## sep-cmaes.R - covariance matrix adapting evolutionary strategy for separable problem
##
##' Global optimization procedure using a covariance matrix adapting
##' evolutionary strategy for separable functions.
##'
##' @source The code is a minor modification from cmaes::cma_es(), adding cutoff on the box constraint. (line 128)
##'
##' @seealso \code{\link[cmaes]{cma_es}}
##'
##' @title Covariance matrix adapting evolutionary strategy for separable problem
##' @export
cma_es <- function(par, fn, ..., lower, upper, control=list(), logFeasible=F,limit_sigma=F) {
  norm <- function(x)
    drop(sqrt(crossprod(x)))

  # set.seed(100)
  # print(control)
  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v))
      return (default)
    else
      return (v)
  }

  ## Initial solution:
  xmean <- drop(par)
  N <- length(xmean)
  ## Box constraints:
  if (missing(lower))
    lower <- rep(-Inf, N)
  else if (length(lower) == 1)
    lower <- rep(lower, N)

  if (missing(upper))
    upper <- rep(Inf, N)
  else if (length(upper) == 1)
    upper <- rep(upper, N)

  range <- max(upper-lower)

  ## Parameters:
  trace       <- controlParam("trace", FALSE)
  fnscale     <- controlParam("fnscale", 1)
  stopfitness <- controlParam("stopfitness", -Inf)
  maxiter     <- controlParam("maxit", 100 * N^2)
  sigma       <- controlParam("sigma", 0.5*max(upper-lower))
  sc_tolx     <- controlParam("stop.tolx", 1e-12 * sigma) ## Undocumented stop criterion
  keep.best   <- controlParam("keep.best", TRUE)
  vectorized  <- controlParam("vectorized", FALSE)

  ## Logging options:
  log.all    <- controlParam("diag", FALSE)
  log.sigma  <- controlParam("diag.sigma", log.all)
  log.eigen  <- controlParam("diag.eigen", log.all)
  log.value  <- controlParam("diag.value", log.all)
  log.pop    <- controlParam("diag.pop", log.all)

  ## Strategy parameter setting (defaults as recommended by Nicolas Hansen):
  lambda      <- controlParam("lambda", 4+floor(3*log(N)))
  mu          <- controlParam("mu", floor(lambda/2))
  weights     <- controlParam("weights", log(mu+1) - log(1:mu))
  weights     <- weights/sum(weights)
  mueff       <- controlParam("mueff", sum(weights)^2/sum(weights^2))
  cc          <- controlParam("ccum", 4/(N+4))
  cs          <- controlParam("cs", (mueff+2)/(N+mueff+3))
  mucov       <- controlParam("ccov.mu", mueff)
  ccov        <- controlParam("ccov.1",
                              (1/mucov) * 2/(N+1.4)^2
                              + (1-1/mucov) * min(1,((2*mucov-1)/((N+2)^2+mucov))))
  damps       <- controlParam("damps",
                              1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs)
  C           <- controlParam("cov",NULL)
  max_no_improve <- controlParam("noImprove", 10+ceiling(30*N/lambda))
  TolUpSigma  <- controlParam("tolUpSigma", 1e5)

  term_code   <- 0 # no error
  no_improve_count <- 0
  ## Safety checks:
  stopifnot(length(upper) == N)
  stopifnot(length(lower) == N)
  stopifnot(all(lower < upper))
  stopifnot(length(sigma) == 1)

  ## Bookkeeping variables for the best solution found so far:
  # best.fit <- Inf
  # best.par <- NULL
  best.fit <- fn(par, ...) * fnscale
  best.par <- par

  best_arfit <- Inf

  best.fit_cut <- fn(par, ...) * fnscale
  best.par_cut <- par

  starting_sigma <- sigma
  ## Preallocate logging structures:
  if (log.sigma)
    sigma.log <- numeric(maxiter)
  if (log.eigen)
    eigen.log <- matrix(0, nrow=maxiter, ncol=N)
  if (log.value)
    value.log <- matrix(0, nrow=maxiter, ncol=mu)
  if (log.pop)
    pop.log <- array(0, c(N, mu, maxiter))

  ## Initialize dynamic (internal) strategy parameters and constants
  # pc <- rep(0.0, N)
  # ps <- rep(0.0, N)
  pc <- controlParam("pc", rep(0.0, N))
  ps <- controlParam("ps", rep(0.0, N))
  chiN <- sqrt(N) * (1-1/(4*N)+1/(21*N^2))
  if(is.null(C)){
    B <- diag(N)
    D <- diag(N)
    BD <- B %*% D
    C <- BD %*% t(BD)
  }else{
    e <- eigen(C, symmetric = TRUE)
    if (log.eigen)
      eigen.log[iter, ] <- rev(sort(e$values))
    if (!all(e$values >= sqrt(.Machine$double.eps) * abs(e$values[1]))) {
      msg <- "Covariance matrix 'C' is numerically not positive definite."
    }
    sigma <- sigma * exp((norm(ps)/chiN - 1)*cs/damps)
    B <- e$vectors
    D <- diag(sqrt(e$values), length(e$values))
    BD <- B %*% D
  }


  iter <- 0L      ## Number of iterations
  counteval <- 0L ## Number of function evaluations
  cviol <- 0L     ## Number of constraint violations
  msg <- NULL     ## Reason for terminating
  nm <- names(par) ## Names of parameters

  ## Preallocate work arrays:
  arx <- matrix(0.0, nrow=N, ncol=lambda)
  arfitness <- numeric(lambda)
  while (iter < maxiter) {
    iter <- iter + 1L

    if (!keep.best) {
      best.fit <- Inf
      best.par <- NULL
      best.fit_cut <- Inf
      best.par_cut <- NULL
    }


    if (log.sigma)
      sigma.log[iter] <- sigma


    ## Generate new population:

    arz <- matrix(rnorm(N*lambda), ncol=lambda)
    arx <- xmean + sigma * (BD %*% arz)

    vx <- ifelse(arx > lower, ifelse(arx < upper, arx, upper), lower)

    # cutoff to force feasibility # not written in Olaf Mersmann's version
    # arx <- vx

    if (!is.null(nm))
      rownames(vx) <- nm

    pen <- 1 + colSums((arx - vx)^2)
    pen[!is.finite(pen)] <- .Machine$double.xmax / 2
    cviol <- cviol + sum(pen > 1)

    if (vectorized) {
      y <- fn(vx, ...) * fnscale
    } else {
      y <- apply(vx, 2, function(x) fn(x, ...) * fnscale)
    }
    counteval <- counteval + lambda

    arfitness <- y * pen
    valid <- pen <= 1

    if (any(valid)) {
      wb <- which.min(y[valid])
      if (y[valid][wb] < best.fit) {
        best.fit <- y[valid][wb]
        best.par <- arx[,valid,drop=FALSE][,wb]
      }
    }
    wb_cut <- which.min(y)
    if (y[wb_cut] < best.fit_cut) {
      best.fit_cut <- y[wb_cut]
      best.par_cut <- vx[,wb_cut]
    }

    ## Order fitness:
    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]

    aripop <- arindex[1:mu]
    selx <- arx[,aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[,aripop]
    zmean <- drop(selz %*% weights)

    # check for arfit improvement
    if( arfitness[1] < best_arfit){
      best_arfit <- arfitness[1]
      no_improve_count <- 0
    }else{
      no_improve_count <- no_improve_count + 1
    }

    ## Save selected x value:
    if(!logFeasible){
      if (log.pop) pop.log[,,iter] <- selx
      if (log.value) value.log[iter,] <- arfitness[aripop]
    }else{
      if (log.pop) pop.log[,,iter] <- vx[,aripop]
      if (log.value) value.log[iter,] <- y[aripop]
    }

    ## Cumulation: Update evolutionary paths
    ps <- (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B %*% zmean)
    hsig <- drop((norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN) < (1.4 + 2/(N+1)))
    pc <- (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * drop(BD %*% zmean)

    ## Adapt Covariance Matrix:
    BDz <- BD %*% selz
    C <- (1-ccov) * C + ccov * (1/mucov) *
      (pc %o% pc + (1-hsig) * cc*(2-cc) * C) +
      ccov * (1-1/mucov) * BDz %*% diag(weights) %*% t(BDz)

    ## Adapt step size sigma:
    e <- eigen(C, symmetric = TRUE)
    if (log.eigen)
      eigen.log[iter, ] <- rev(sort(e$values))
    if (!all(e$values >= sqrt(.Machine$double.eps) * abs(e$values[1]))) {
      msg <- "Covariance matrix 'C' is numerically not positive definite."
      break
    }

    sigma <- sigma * exp((norm(ps)/chiN - 1)*cs/damps)
    # limiting sigma
    if(limit_sigma)
      if(sigma>=0.1*range) sigma <- 0.1*range
    # D <- diag(N)*C
    # D <- sqrt(D)
    B <- e$vectors
    D <- diag(sqrt(e$values), length(e$values))
    BD <- B %*% D


    ## break if fit:
    if (arfitness[1] <= stopfitness * fnscale) {
      msg <- "Stop fitness reached."
      break
    }

    ## Check stop conditions:

    ## Condition 1 (sd < tolx in all directions):
    if (all(D < sc_tolx) && all(sigma * pc < sc_tolx)) {
      msg <- "All standard deviations smaller than tolerance."
      break
    }

    ## Escape from flat-land:
    if (arfitness[1] == arfitness[min(1+floor(lambda/2), 2+ceiling(lambda/4))]) {
      sigma <- sigma * exp(0.2+cs/damps);
      if (trace)
        message("Flat fitness function. Increasing sigma.")
    }
    if (trace)
      message(sprintf("Iteration %i of %i: current fitness %f",
                      iter, maxiter, arfitness[1] * fnscale))

    if(sigma>=starting_sigma*TolUpSigma){
      message(sprintf("Sigma divergence detected! Terminating."))
      term_code <- 1
      break
    }

    if(no_improve_count >= max_no_improve){
      message(sprintf("No fitness improvement for %i generation. Terminating.",max_no_improve))
      term_code <- 2
      break
    }

    if(checkNoEffectAxis(xmean,BD,sigma)){
      message(sprintf("Addition of 0.1 times sigma does not change mean value."))
      term_code <- 3
      break
    }

    if(checkNoEffectCoord(xmean,sigma)){
      message(sprintf("Addition of 0.2 times sigma on any variable does not change mean value."))
      term_code <- 4
      break
    }

    if(checkCondNumber(C)){
      message(sprintf("Cov matrix condition number exceed 1e14."))
      term_code <- 5
      break
    }
  }
  cnt <- c(`function`=as.integer(counteval), gradient=NA)

  log <- list()
  ## Subset lognostic data to only include those iterations which
  ## where actually performed.
  if (log.value) log$value <- value.log[1:iter,]
  if (log.sigma) log$sigma <- sigma.log[1:iter]
  if (log.eigen) log$eigen <- eigen.log[1:iter,]
  if (log.pop)   log$pop   <- pop.log[,,1:iter]

  ## Drop names from value object
  names(best.fit) <- NULL
  if(logFeasible){
    best.fit <- best.fit_cut
    best.par <- best.par_cut
  }

  res <- list(par=best.par,
              value=best.fit / fnscale,
              counts=cnt,
              convergence=ifelse(iter >= maxiter, 1L, 0L),
              message=msg,
              constr.violations=cviol,
              diagnostic=log,
              cov=C,
              sigma=sigma,
              pc=pc,
              ps=ps,
              termination_code=term_code)
  # termination codes:
  # 0: no error
  # 1: sigma divergence
  # 2: no improvement for 20 generations
  # 3: sigma is too small to have effect on covariance matrix
  # 4: sigma is too small to have effect on mean x
  # 5: condition number too large
  class(res) <- "cma_es.result"
  return(res)
}

checkNoEffectAxis <-function(xmean,BD,sigma){
 m <- xmean
 return(sum((m - (m + 0.1 * sigma * diag(C) ))^2) <
          .Machine$double.eps)
}

checkNoEffectCoord <-function(xmean,sigma){
  m <- xmean
  return(sum((m - (m + 0.2 * sigma))^2) < .Machine$double.eps)
}

checkCondNumber <- function(C){
  return(kappa(C)>1e14)
}

