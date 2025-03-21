#' Initial setup function.
#'
#' Set up Julia environment, and install Julia if necessary.
#' Install required Julia packages.
#' @param installJulia boolean if TRUE, it will install Julia
#' @param JULIA_HOME (Optional) string indicating the location of specific Julia binary
#' @export
setup<-function(installJulia = FALSE, JULIA_HOME = NULL){
  julia_setup(installJulia = installJulia, JULIA_HOME = JULIA_HOME)
  julia_library("Pkg")
  julia_do.call("Pkg.add", list(url="https://github.com/Julia-Tempering/Pigeons.jl"))
  julia_do.call("Pkg.add", list(url="https://github.com/TuringLang/JuliaBUGS.jl"))
  julia_do.call("Pkg.add", list(url = "https://github.com/JuliaIO/FFMPEG.jl"))
  julia_install_package_if_needed("MCMCChains")
  julia_install_package_if_needed("StatsPlots")
}

#' Initialize Julia session
#'
#' It needs to be required to run before any functions in this package
#' The number of threads in Julia can't be changed once session starts.
#' To restart Julia with new number of threads, restart of R session and run this function.
#' @param threads The number of threads used in Julia session
#' @export
init_julia_session_multi<-function(threads){
  Sys.setenv("JULIA_NUM_THREADS" = threads)
  if(julia_eval("Threads.nthreads()") != threads){
    warnings("threads options are not properly set, please restart R session.")
  }
}
