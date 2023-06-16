#' Compute a fitting autocorrelation structure for the base process such that
#' it approximates the provided target process, by use of a interpolation method.
#' @details A map of base process autocorrelation values to target process
#'   autocorrelation values is stored in the \code{.interpol_caches} file in the
#'   current working directory, if the 'use_caching' parameter is set to TRUE.
#'   Marginal distributions that assign significant
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
#' @param gamma Threshold for computation of the expected value of the
#'   product of two random variables of the target process.
#' @param n_interpol Controls the number of equidistant points in the base autocorrelation space for which resulting autocorrelations in the target autocorrelation space are computed, in order to fit a polynomial through them.
#' @param poly_deg Sets the degree of the polynomial fit through 'n_interpol' many points.
#' @param use_caching Should cached values be used, and saved after computation, or not.
#' @param cdf_name_parameterized Name of the target Distribution, together with
#'   its parameters, separated by dashes. Used for naming the cache-file, which
#'   stores a map of values from the base process autocorrelation space to the
#'   target process autocorrelation space.
#'
#' @return A vector containing a fitting autocorrelation structure for the base
#'   process.
#'
find_r_interpol <- function(cdf, mean, var, rho, gamma, n_interpol, poly_deg = 9, use_caching = T, cdf_name_parameterized){
  lower_bound_fitting = if(min(rho)<0) -1 else 0
  upper_bound_fitting = if(max(rho)>0) 1 else 0

  # check if autocorrelation values are feasible
  if(use_caching){
    full_name_parameterized = paste(cdf_name_parameterized, gamma, sep = "_")
    cache_file <- file.path(".interpol_caches", full_name_parameterized)
    bound_cache_file <- file.path(".interpol_caches", paste("correlation_bound", full_name_parameterized, sep = "_"))
    search_boundaries<-if(file.exists(bound_cache_file)) readRDS(bound_cache_file) else  get_correlation_bound(cdf=cdf, mean=mean, var=var, gamma = gamma)
    dir.create(".interpol_caches", showWarnings = F)
    saveRDS(object = search_boundaries, file = bound_cache_file)
  }else{
    search_boundaries <- get_correlation_bound(cdf=cdf, mean=mean, var=var, gamma = gamma)
  }
  if(any(rho<search_boundaries[1]) | any(rho>search_boundaries[2])){
    stop(paste("correlation boundary violated for gamma ",gamma, ", resulting bounds: [", search_boundaries[1],", ", search_boundaries[2],"]", sep = ""))
  }

  # calculate values for fitting polynomial to correlation of target series
  fitting_base <- seq(lower_bound_fitting,upper_bound_fitting, length.out = n_interpol)
  if(use_caching){
    poly_target <-
      get_poly_target_cached(
        cache_file = cache_file,
        fitting_base = fitting_base,
        cdf = cdf,
        gamma  = gamma,
        mean = mean,
        var = var
      )
  }else{
    poly_target <-
      get_poly_target(
        fitting_base = fitting_base,
        cdf = cdf,
        gamma  = gamma,
        mean = mean,
        var = var
      )
  }
  # fit polynomial
  interpol_base=seq(lower_bound_fitting, upper_bound_fitting, 0.0002)
  poly_function = pracma::polyfit(x = fitting_base, y = as.numeric(poly_target), n = poly_deg)
  # generate values from which to interpolate autocorrelation
  interpol_target=pracma::polyval(poly_function ,interpol_base)
  interpol_values=cbind(interpol_base,interpol_target)
  # interpolate autocorrelation of target series
  r=approx(interpol_values[,2],interpol_values[,1], rho, method = "linear", rule = 2)$y
  return(r)
}

#' Compute values to which to fit the polynomial, without caching
#'
get_poly_target <- function(fitting_base, cdf, gamma, mean, var){
  fitting_intermed_res<- sapply(fitting_base,FUN=expected_target_product, cdf = cdf, gamma  = gamma)
  target_vals <- sapply(fitting_intermed_res, FUN=get_target_correlation,mean=mean,  var=var)
  return(target_vals)
}

#' Compute values to which to fit the polynomial, with caching enabled
#'
get_poly_target_cached <- function(cache_file, fitting_base, cdf, gamma, mean, var){

  cache <- if(file.exists(cache_file)) readRDS(cache_file) else r2r::hashmap(default = -2)
  # either compute target vals to which to fit polynomial function, or load them from the cache
  target_vals <- sapply(fitting_base, FUN = get_target_val_cached, cache = cache, cdf = cdf, gamma = gamma, mean=mean, var=var)
  cache[fitting_base] = target_vals
  if(!file.exists(".interpol_caches")){
    dir.create(".interpol_caches", showWarnings = F)
  }
  saveRDS(cache, file = cache_file, compress = T)
  return(target_vals)
}

#' Compute single autocorrelation value in target autocorrelation space from value in base autocorrelation space
get_target_val_cached <-
  function(base_val, cache, cdf, gamma, mean, var) {
    cached_val <- cache[base_val]
    target_val <-
      if (cached_val != -2) { # cache hit
        cached_val
      } else{
        # cache miss
        get_target_correlation(
          expected_target_product = expected_target_product(
            cdf = cdf,
            r = base_val,
            gamma = gamma
          ),
          mean = mean,
          var = var
        )
      }
  }
