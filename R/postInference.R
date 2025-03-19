#' Extract samples from pigeons results
#'
#' This function returns the matrix of samples from pigeons results.
#'
#' @param pt pigeons results from sampling function
#' @return a matrix of samples
#' @export
pigeons_samples<-function(pt){
  julia_assign("result",pt)
  julia_library("MCMCChains")
  julia_command("samples = Chains(result)")
  samplesMatrix = julia_eval("Array(samples)")
  names = julia_eval("sample_names(pt)")
  colnames(samplesMatrix) = names
  return(samplesMatrix)
}


#' Trace and density plot from pigeons results
#'
#' This function returns the plots from pigeons results.
#'
#' @param pt pigeons results from sampling function
#' @export
pigeons_plots<-function(pt){
  julia_library("MCMCChains")
  julia_assign("result",pt)
  julia_command("samples = Chains(result)")
  julia_library("MCMCChains")
  julia_library("StatsPlots")
  julia_command("StatsPlots.plot(samples)")
}


