#'@import VGAM
library(VGAM)
#'@import r2r
library(r2r)
#'@import polynom
library(polynom)
#'@import R.filesets
library(R.filesets)
#'@import mvtnorm
library(mvtnorm)
#'@import purrr
library(purrr)
#'@import pracma
library(pracma)

#'Computes the value of the integral for a plateau on the integrand function
#'surface. Plateau is specified by the values of the target distribution value
#'range that define its edges.
#'
#'@param j Element of target distribution's value range.
#'@param k Element of target distribution's value range.
#'@param cdf Target distribution's Cumulative Distribution Function. All of its
#'  parameters should be set, so the argument is unambiguous.
#'@param r Autocorrelation between two variables of the base process. Should be
#'  selected such that it induces the targeted autocorrelation in the target
#'  process.
#' @seealso \code{\link{VGAM::pbinorm}}
integrand_plateau <- function(j, k, cdf, r) {
  val <- k * j * (pbinorm(qnorm(cdf(k)), qnorm(cdf(j)), cov12 = r) -
                    pbinorm(qnorm(cdf(k)), qnorm(cdf(j - 1)), cov12 = r) -
                    pbinorm(qnorm(cdf(k - 1)), qnorm(cdf(j)), cov12 = r) +
                    pbinorm(qnorm(cdf(k - 1)), qnorm(cdf(j - 1)), cov12 = r))
  val[val < 0] <- 0
  return(val)
}

#'Computes the value of the integral for a plateau on the integrand function
#'surface. Plateaus edges have to be defined by the same value from the target
#'distribution's value range.
#'
#'@param k Element of target distribution's value range.
#'@param cdf Target distribution's Cumulative Distribution Function. All of its
#'  parameters need to be set, so that an element from the value range is the
#'  only argument to provided.
#'@param r Autocorrelation between two variables of the base process. Should be
#'  selected such that it induces the targeted autocorrelation in the target
#'  process.
#'
integrand_plateau_symmetric<- function(k, cdf, r) {
  (k^2) * (pbinorm(qnorm(cdf(k)), qnorm(cdf(k)), cov12 = r) -
             2 * pbinorm(qnorm(cdf(k)), qnorm(cdf(k-1)), cov12 = r) +
             pbinorm(qnorm(cdf(k-1)), qnorm(cdf(k-1)), cov12 = r))
}

#'Computes an approximation of the integral over the integral function. This is
#'done by summing up the exact values of integrals over plateaus, until a
#'threshold is reached.
#'@param cdf Target distribution's Cumulative Distribution Function. All of its
#'  parameters need to be set, so that an element from the value range is the
#'  only argument to provided.
#'@param r Autocorrelation between two variables of the base process. Should be
#'  selected such that it induces the targeted autocorrelation in the target
#'  process.
#'@param precision Controls the quality of the integral approximation.
#'  Computation stops once the next subsum of integral plateaus being computed
#'  has a relative size compared to the current total sum smaller then the
#'  precision.
#'
expected_target_product <- function(cdf,r, precision){
  k <- 1
  total <- 0
  s <- 0
  s_max <- 0
  while(s >= precision*total | k < 50){
    s <- 2*sum(integrand_plateau(seq(1, k-1), k, cdf, r))
    s <- s + integrand_plateau_symmetric(k, cdf, r)
    # Avoid numeric errors by not allowing negative results, which should not be possible
    if(s < 0){
      s <- 0
    }
    total <- total + s
    k <- k + 1
    if(cdf(k-1) == cdf(k)){
      break
    }
  }
  return(total)
}

#'Helper function for computing the autocorrelation between to random variables
#'of the target process based on the properties of the corresponding random
#'variables of the currently evaluated base process, i.e., the selected base
#'process autocorrelation and target process CDF.
#'@param expected_target_product Approximation (or exact value, if possible) of
#'  the integral function defining the expected value of two random variables of
#'  the target process. Usually expects the return value of
#'  \link{expected_target_product}.
#'@param mean Mean of the target distribution.
#'@param var Variance of the target distribution.
get_target_correlation <- function(expected_target_product,mean, var){
  (expected_target_product - mean^2)/var
}


#'Computes the minimum and maximum possible values for autocorrelation between
#variables of the target process. If any target autocorrelation present in the
#target autocorrelation structures lies outside of the computed bounds, the
#process cannot produce valid results. '@param cdf Target distribution's
#Cumulative Distribution Function. All of its '  parameters need to be set, so
#that an element from the value range is the '  only argument to provided.
#'@param mean Mean of the target distribution. '@param var Variance of the
#target distribution. '@param precision Defines a threshhold for the relative
#size of the subsummand in the computation of function
#\code{\link{expected_target_product}}

get_correlation_bound<- function(cdf, mean, var, precision){
  return(c(get_target_correlation(expected_target_product = expected_target_product(cdf = cdf,r = -1, precision = precision), mean=mean, var=var), get_target_correlation(expected_target_product(cdf = cdf,r = 1, precision = precision), mean=mean, var=var)))
}


#' Computes a fitting autocorrelation structure for the base process such that
#' it approximates the provided target process.
#' @details A map of base process autocorrelation values to target process
#'   autocorrelation values is stored in the \code{.DARTA_caches} file in the
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
#' @param precision Threshold for computation of the expected value of the
#'   product of two random variables of the target process.
#' @param print_progress Flag indicating whether updates on the generation
#'   process should be printed. Used for debugging, or estimating computation
#'   time.
#'
#' @return A vector containing a fitting autocorrelation structure for the base
#'   process, or \code{NULL} if the process was not successful.
#'
find_r <- function(cdf, mean, var, rho, cdf_name_parameterized, epsilon, precision, print_progress=FALSE){

  r <- vector(mode="double", length = length(rho))
  approximation <- vector(mode="double", length = length(rho))
  
  cache_dir <- ".DARTA_caches"
  cache_path <- file.path(cache_dir,paste(cdf_name_parameterized, "precision",precision, sep = "_"))
  boundary_cache_path <-file.path(cache_dir, paste("correlation_bound",cdf_name_parameterized,precision, sep = "_"))
  if(file.exists(cache_path)){
    cache <- loadRDS(cache_path)
  } else{
    dir.create(cache_dir, showWarnings = F)
    cache <- hashmap(default=-2)
  }
  if(file.exists(boundary_cache_path)){
    search_boundaries<- loadRDS(boundary_cache_path)
  }else{
    search_boundaries <- get_correlation_bound(cdf=cdf, mean=mean, var=var, precision = precision)
    dir.create(cache_dir, showWarnings = F)
    saveRDS(search_boundaries, file = boundary_cache_path)
  }
  for(i in 1:length(rho)){
    if(rho[i]<search_boundaries[1] | rho[i]>search_boundaries[2]){
      print(paste(cdf_name_parameterized,":correlation boundary violated by value ",rho[i]," for precision ",precision, ", resulting bounds: ", search_boundaries[1],"|", search_boundaries[2], sep = " "))
      return(NULL)
    }
  }
  for(i in 1:length(rho)){
    interval <- c(-1,1)
    r_found = FALSE
    search_limiter = 0
    if(print_progress){
      print(paste0("rho[",i,"]:",rho[i]))
    }
    if(rho[i] == 0){
      r[i] <- 0
      approximation[i] <- 0
      r_found <-  TRUE
    }
    while(!r_found & search_limiter < 20){
      search_limiter <- search_limiter+1
      search_r <- (interval[1]+interval[2])/2
      current <-cache[[search_r]]
      if(current==-2){ #cache miss
        current <- get_target_correlation(expected_target_product(cdf = cdf, r = search_r, precision = precision), mean=mean,  var=var)
        if(!is.null(current)){
          cache[[search_r]]<-current
          saveRDS(cache, file = cache_path)
        }else{
          print(paste("An Error occurred during computation of the target process correlation."))
          return(NULL)
        }
      }
      if(print_progress){
        print(paste0("r[",i,"]:", search_r, " || approximation for rho: ", current))
      }
      if(abs(current - rho[i])<epsilon){
        r_found = TRUE
        r[i] <- search_r
        approximation[i] <- current
      }else{
        if(rho[i]<current){
          interval <- c(interval[1], search_r)
        }
        else{
          interval <- c(search_r, interval[2])
        }
      }
    }
  }
  return(data.frame(r, rho, approximation))
}

#' Determines Gamma matrix for solving Yule-Walker-Equations.
#' @param r Vector of \code{numerics} describing the autocorrelation structure
#'   of the base process.
get_gamma <- function(r){
  return(toeplitz(c(1,r))[1:length(r),1:length(r)])
}

#'Check if distribution resulting from alpha is stationary.
#'@param alpha Vector of \code{numerics} defining the coefficients in a
#'  recursively defined stochastic process.
#'@return Logical indicating if the resulting distribution is stationary or not.
is_stationary <- function(alpha){
  roots <- solve(polynomial(c(1, -alpha)))
  roots <- abs(roots)
  is_stationary <- TRUE
  if(sum(roots < 1)){
    is_stationary <- FALSE
  }
  return(is_stationary)
}

#' Generate autocorrelated time-series with Negative Binomial marginal distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param size Size parameter of the Negative Binomial distribution.
#' @param prob Probability parameter of the Negative Binomial distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process. Defaults to
#'   0.001.
#' @param precision Threshold for computation of the expected value of the
#'   product of two random variables of the target process. defaults to 10^-5.
#'
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Negative Binomial marginal distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_DARTA}}
generate_nbinomial <- function(n, size, prob, rho, epsilon = 0.001, precision = 0.00001, method = "interpol" ){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "nbinomial", param1 = size, param2=prob,epsilon = epsilon, precision = precision, method = method))
}

#' Generate autocorrelated time-series with Binomial marginal distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param size Size parameter of the Binomial distribution.
#' @param prob Probability parameter of the Binomial distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process. Defaults to
#'   0.001.
#' @param precision Threshold for computation of the expected value of the
#'   product of two random variables of the target process. defaults to 10^-5.
#'
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Binomial marginal distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_DARTA}}
generate_binomial <- function(n, size, prob, rho, epsilon = 0.001, precision = 0.00001, method = "interpol" ){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "binomial", param1 = size, param2=prob, epsilon = epsilon, precision = precision, method = method ))
}

#' Generate autocorrelated time-series with marginal poisson-distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param lambda Lambda parameter of the Poisson distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process. Defaults to
#'   0.001.
#' @param precision Threshold for computation of the expected value of the
#'   product of two random variables of the target process. defaults to 10^-5.
#'
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal Poisson-distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_DARTA}}
generate_poisson <- function(n, lambda, rho, epsilon = 0.001, precision = 0.00001, method = "interpol" ){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "poisson", param1 = lambda,epsilon = epsilon, precision = precision, method = method ))
}

#' Generate time-series of the provided marginal distribution and
#' autocorrelation structure.
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   numerics.
#' @param distribution_name Name of the Distribution to be generated. Currently
#'   supported are "nbinomial", "binomial", "poisson", and "uniform".
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process. Defaults to
#'   0.001.
#' @param precision Threshold for computation of the expected value of the
#'   product of two random variables of the target process. defaults to 10^-5.
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal distribution defined by distribution_name, or NULL if base
#'   process is not stationary.
#' @examples
#' generate_DARTA(n = 10000, rho = c(0.7,0.5,0.3), distribution_name = "nbinomial", param1 = 10, param2 = 0.6)
#' @seealso \code{\link{gen_DARTA}}
#' @export
generate_DARTA <- function(n, rho, distribution_name = c("nbinomial", "binomial","poisson", "uniform"), epsilon = 0.001, n_interpol = 20, param1 = NULL, param2=NULL, precision = 0.00001, method = c("interpol", "binary"), use_caching = T){
  if(distribution_name == "nbinomial"){
    mean = param1*(1-param2)/param2
    var = param1*(1-param2)/(param2**2)
    cdf_name_parameterized <- paste("nbinomial",as.character(param1),as.character(param2),sep = "-")
    return(gen_DARTA(n = n, cdf = partial(pnbinom, size = param1, prob = param2),inv = partial(qnbinom, size = param1, prob = param2), mean = mean, var = var, cdf_name_parameterized = cdf_name_parameterized, rho = rho, epsilon = epsilon,n_interpol = n_interpol, method = method, precision = precision, use_caching = use_caching))
  }
  if(distribution_name == "binomial"){
    mean = param1*param2
    var = param1*param2*(1-param2)
    cdf_name_parameterized <-paste("binomial",as.character(param1),as.character(param2),sep = "-")
    return(gen_DARTA(n = n, cdf = partial(pbinom, size = param1, prob = param2),inv = partial(qbinom, size = param1, prob = param2), cdf_name_parameterized = cdf_name_parameterized,mean = mean, var = var, rho = rho, epsilon = epsilon,n_interpol = n_interpol, method = method,  precision = precision, use_caching = use_caching))
  }
  if(distribution_name == "poisson"){
    mean = param1
    var = param1
    cdf_name_parameterized <- paste("poisson",as.character(param1),sep = "-")
    return(gen_DARTA(n = n, cdf = partial(ppois, lambda = param1),inv = partial(qpois, lambda = param1), cdf_name_parameterized = cdf_name_parameterized, mean = mean, var = var, rho = rho, epsilon = epsilon,n_interpol = n_interpol, method = method, precision = precision, use_caching = use_caching))
  }
  if(distribution_name == "uniform"){
    mean = (param2+param1)/2
    var = ((param2- param1)**2 -1)/12
    cdf_name_parameterized <- paste("uniform",as.character(param1),as.character(param2),sep = "-")
    return(gen_DARTA(n = n, cdf = partial(pdunif, min = param1, max = param2), inv = partial(qdunif, min = param1,max = param2), cdf_name_parameterized = cdf_name_parameterized, mean = mean, var = var, rho = rho, epsilon = epsilon,n_interpol = n_interpol, method = method, precision = precision, use_caching = use_caching))
  }
}


#' Creates time-series of target distribution and autocorrelation structure, if
#' possible.
#' @details It is generally encouraged to use \code{\link{generate_DARTA}}
#'   instead, which calls upon this function.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param cdf Marginal CDF of target process.
#' @param inv Inverse of marginal CDF of target process, a.k.a. quantile
#'   function.
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
#' @param precision Threshold for computation of the expected value of the
#'   product of two random variables of the target process.
#' @examples
#' size = 6
#' prob = 0.4
#' # Compute mean and variance for the negative binomial distribution
#' mean = size*(1-prob)/prob
#' # 9
#' var = size*(1-prob)/(prob**2)
#' # 22.5
#' gen_DARTA(n = 10, cdf = partial(pnbinom, size = 6, prob = 0.4), inv = partial(qnbinom, size = 6, prob = 0.4), mean = 9, var = 22.5, rho = c(0.6,0.2,0.1), cdf_name_parameterized = "nbinomial-6-0.4", precision = 0.000001, epsilon = 0.0001)
#' # [1] 12 19 19 11 16 15 23 21 21 14
#' @seealso \code{\link{generate_DARTA}}
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal distribution defined by distribution_name, or NULL if base
#'   process is not stationary.
#' @export


gen_DARTA <- function(n, cdf, inv, mean, var, rho, cdf_name_parameterized, epsilon,n_interpol, method, precision, use_caching =T){
  # check if arguments are valid
  if(!(method %in% c("interpol", "binary"))){
    stop(paste("Method should be either 'interpol' or 'binary', but is '", method,"'", sep = ""))
  }
  if(length(n_interpol)>1){ 
    stop(paste("n_interpol should have length 1"))
  }
  if(length(epsilon)>1){
    stop(paste("epsilon should have length 1"))
  }
  if(equals(method, "interpol") & (!is.numeric(n_interpol) | n_interpol <= 0 )){
    stop(paste("When method 'interpol' is selected, n_interpol should be a positive natural number, but is ",n_interpol, sep = ""))
  }
  if(equals(method, "binary") & (epsilon <= 0 | epsilon >= 1)){
    stop(paste("When method 'binary' is selected, epsilon should be a positive number much smaller than 1, but is ", epsilon, sep =""))
  }
  
  r <- switch(method, 
              "interpol" = find_r_interpol(cdf = cdf,cdf_name_parameterized=cdf_name_parameterized, mean = mean, var = var, precision = precision,rho = rho, n_interpol = n_interpol),
              "binary" = find_r(cdf = cdf, mean = mean, var = var, rho = rho, cdf_name_parameterized = cdf_name_parameterized, epsilon = epsilon, precision = precision)[,1])
  if(is.null(r)){
    print("No suitable base process found, return NULL")
    return(NULL)
  }
  Y <- vector(mode = "double", length = n)
  p <- length(r)
  Z <- vector(mode = "numeric", length = p)
  gamma <- get_gamma(r)
  alpha <- solve(gamma, r)
  sigma <- sqrt(1 - alpha %*% r)
  alpha <- rev(alpha)
  if(p == 1){
    init_fun <- rnorm
  }else{
    init_fun <- rmvnorm
  }
  if(is_stationary(alpha)){
    Z <- rev(as.vector(init_fun(1, rep(0, p),gamma)))
    for(i in 1:n){
      Z[p+1] <-alpha %*% Z+ rnorm(1, mean=0,sd= sigma)
      Z <- Z[-1]
      Y[i] <- inv(pnorm(Z[p]))
    }
    return(Y)
  }else{
    print("Base process is not stationary")
  }
}


find_r_interpol <- function(cdf, cdf_name_parameterized, mean, var, precision, rho, n_interpol, poly_deg = 9, use_caching = T){
  lower_bound_fitting = if(min(rho)<0) -1 else 0
  upper_bound_fitting = if(max(rho)>0) 1 else 0
  
  # check if autocorrelation values are feasible
  if(use_caching){
    full_name_parameterized = paste(cdf_name_parameterized, precision, sep = "_")
    cache_file <- file.path(".interpol_caches", full_name_parameterized)
    bound_cache_file <- file.path(".interpol_caches", paste("correlation_bound", full_name_parameterized, sep = "_"))
    search_boundaries<-if(file.exists(bound_cache_file)) loadRDS(bound_cache_file) else  get_correlation_bound(cdf=cdf, mean=mean, var=var, precision = precision) 
    dir.create(".interpol_caches", showWarnings = F)
    saveRDS(object = search_boundaries, file = bound_cache_file)
  }else{
    search_boundaries <- get_correlation_bound(cdf=cdf, mean=mean, var=var, precision = precision)
  }
  if(any(rho<search_boundaries[1]) | any(rho>search_boundaries[2])){
    stop(paste("correlation boundary violated for precision ",precision, ", resulting bounds: [", search_boundaries[1],", ", search_boundaries[2],"]", sep = ""))
    return(NULL)
  }
  
  # calculate values for interpolating correlation of target series
  fitting_base <- seq(lower_bound_fitting,upper_bound_fitting, length.out = n_interpol)
  if(use_caching){
    interpol_target <- get_interpol_target_cached(cache_file= cache_file,fitting_base = fitting_base, cdf = cdf, precision  = precision, mean = mean, var = var)
  }else{
    interpol_target <- get_interpol_target(fitting_base = fitting_base, cdf = cdf, precision  = precision, mean = mean, var = var)
  }
  
  # interpolate autocorrelation of target series
  interpol_base=seq(lower_bound_fitting, upper_bound_fitting, 0.0002)
  interpol_function = pracma::polyfit(x = fitting_base, y = as.numeric(interpol_target), n = poly_deg)
  interpol_target=pracma::polyval(interpol_function ,interpol_base)
  interpol_values=cbind(interpol_base,interpol_target)
  r=approx(interpol_values[,2],interpol_values[,1], rho, method = "linear", rule = 2)$y
  return(r)
}

get_interpol_target <- function(fitting_base, cdf, precision, mean, var){
  fitting_intermed_res<- sapply(fitting_base,FUN=expected_target_product, cdf = cdf, precision  = precision)
  target_vals <- sapply(fitting_intermed_res, FUN=get_target_correlation,mean=mean,  var=var)
  return(target_vals)
}

get_interpol_target_cached <- function(cache_file, fitting_base, cdf, precision, mean, var){
  cache <- if(file.exists(cache_file)) loadRDS(cache_file) else hashmap(default = -2)
  target_vals <- sapply(fitting_base, FUN = get_target_val_cached, cache = cache, cdf = cdf, precision = precision, mean=mean, var=var)
  cache[fitting_base] = target_vals
  if(!file.exists(".interpol_caches")){
    dir.create(".interpol_caches", showWarnings = F)
  }
  saveRDS(cache, file = cache_file, compress = T)
  return(target_vals)
}

get_target_val_cached <-
  function(base_val, cache, cdf, precision, mean, var) {
    cached_val <- cache[base_val]
    target_val <-
      if (cached_val != -2) {
        cached_val
      } else{
        get_target_correlation(
          expected_target_product = expected_target_product(
            cdf = cdf,
            r = base_val,
            precision = precision
          ),
          mean = mean,
          var = var
        )
      }
  }


