#' Main function for the cooperative coevolution procedure.
#'
#' @title Cooperative coevolution.
#' @param population Initial population
#' @param fun The objective function to be solved
#' @param group Vector of list. Each list contains a group of non-separable variables.
#' @export
cc <- function(population,fun,group=NULL){
  # grouping
  if(is.null(group)){ # if no group is supplied, then the group is determined by differential grouping based on the first individual
    group <- differential_grouping(population[,1],fun)
  }

  # error checking on groups
  if(!is.list(group)) stop('group is of wrong mode, it should be a list.')
  if(any(summary(group)[,3]!='vector')) stop('Sublist of group is of wrong mode, all of them should also be a vector')

  # Optimize for each group
  for(groupIndex in 1:length(group)){

  }
}

pygmo <- NULL
rndGen <- NULL

.onLoad <- function(libname, pkgname){
  library(reticulate)

  pygmo <<- reticulate::import("pygmo", delay_load = TRUE)
  rndGen <<- reticulate::import("numpy", delay_load = TRUE)

  have_pygmo <- py_module_available("pygmo")
  if (!have_pygmo)
    print("PyGMO not available, install dependencies using MaOEA::install_python_dependencies()")
}

#' Install the required python package: PyGMO
#' @title Install PyGMO python package
#' @param method Default: auto
#' @param conda Default: auto
#' @export
install_python_dependencies <- function(method = "auto", conda = "auto") {
  reticulate::py_install("pygmo", method = method, conda = conda)
}
