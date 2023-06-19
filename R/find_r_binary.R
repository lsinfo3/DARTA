

#' Compute a fitting autocorrelation structure for the base process such that
#' it approximates the provided target process, by use of a binary search method.
#' @details A map of base process autocorrelation values to target process
#'   autocorrelation values is stored in the \code{.binary_caches} file in the
#'   current working directory. Marginal distributions that assign significant
#'   probability mass to a large amount of their value range tend to require
#'   longer computations in order to find a fitting base autocorrelation
#'   structure, and caching removes the need for this computation on repeated
#'   generation attempts.
#'
#' @param cdf Marginal CDF of target process.
#' @param mean Mean of marginal CDF of target process.
#' @param var Variance of marginal CDF of target process.
#' @param rho Target autocorrelation Structure of target process, provided as a
#'   vector of numerics. Position on the vector describes the lag for which the
#'   numeric autocorrelation value is supposed to occur.
#' @param cdf_name_parameterized Name of the target Distribution, together with
#'   its parameters, separated by dashes. Used for naming the cache-file, which
#'   stores a map of values from the base process autocorrelation space to the
#'   target process autocorrelation space.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process. Defaults to
#'   0.001.
#' @param gamma Threshold for computation of the expected value of the
#'   product of two random variables of the target process.
#' @param print_progress Flag indicating whether updates on the generation
#'   process should be printed. Used for debugging, or estimating computation
#'   time.
#'
#' @return A vector containing a fitting autocorrelation structure for the base
#'   process, or \code{NULL} if the process was not successful.
#'
find_r_binary <- function(cdf, mean, var, rho, cdf_name_parameterized, epsilon, gamma, use_caching = T, print_progress=FALSE){

  r <- vector(mode="double", length = length(rho))
  approximation <- vector(mode="double", length = length(rho))

  if(use_caching){
    cache_dir <- ".binary_caches"
    cache_path <- file.path(cache_dir,paste(cdf_name_parameterized, "gamma",gamma, sep = "_"))
    boundary_cache_path <-file.path(cache_dir, paste("correlation_bound",cdf_name_parameterized,gamma, sep = "_"))
    dir.create(cache_dir, showWarnings = F)
    cache <-if(file.exists(cache_path)) readRDS(cache_path) else r2r::hashmap(default=-2)
    search_boundaries<- if(file.exists(boundary_cache_path)) readRDS(boundary_cache_path) else get_correlation_bound(cdf=cdf, mean=mean, var=var, gamma = gamma)
  }
  else{
    search_boundaries <- get_correlation_bound(cdf=cdf, mean=mean, var=var, gamma = gamma)
  }

  if(use_caching & is.numeric(search_boundaries)){
    saveRDS(search_boundaries, file = boundary_cache_path)
  }

  for(i in 1:length(rho)){
    if(rho[i]<search_boundaries[1] | rho[i]>search_boundaries[2]){
      stop(
        paste(
          cdf_name_parameterized,
          ": correlation boundary violated by value ",
          rho[i],
          " for gamma ",
          gamma,
          ", resulting bounds: [",
          search_boundaries[1],
          ",",
          search_boundaries[2],
          "]",
          sep = " "
        )
      )
    }
  }
  for(i in 1:length(rho)){
    interval <- c(-1,1)
    r_found = FALSE
    search_iteration = 0
    if(print_progress){
      print(paste0("Now searching base correlation for rho[",i,"] = ",rho[i]))
    }
    if(rho[i] == 0){
      r[i] <- 0
      approximation[i] <- 0
      r_found <-  TRUE
    }
    while(!r_found & search_iteration < 30){
      # search was not successful yet
      search_iteration <- search_iteration+1
      search_r <- (interval[1]+interval[2])/2
      if(use_caching){
        current <-cache[[search_r]]
        if(current==-2){ #cache miss
          current <- get_target_correlation(expected_target_product(cdf = cdf, r = search_r, gamma = gamma), mean=mean,  var=var)
          if(!is.null(current)){ # valid result
            cache[[search_r]]<-current
            saveRDS(cache, file = cache_path)
          }else{
            stop(paste("An Error occurred during computation of the target process correlation."))
          }
        }
      }else{
        current <- get_target_correlation(expected_target_product(cdf = cdf, r = search_r, gamma = gamma), mean=mean,  var=var)
      }
      if(print_progress){
        print(paste0("r[",i,"]:", search_r, " || approximation for rho: ", current))
      }
      r[i] <- search_r
      approximation[i] <- current
      if(abs(current - rho[i])<epsilon){
        # approximation successful
        r_found = TRUE
      }else{
        interval <- if(rho[i]<current) c(interval[1], search_r) else c(search_r, interval[2])
      }
    }
  }
  return(data.frame(r, rho, approximation))
}
