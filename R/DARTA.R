
#'Compute the minimum and maximum possible values for autocorrelation between
#'variables of the target process. If any target autocorrelation present in the
#'target autocorrelation structures lies outside of the computed bounds, the
#'process cannot produce valid results.
#'@param cdf Target distribution's Cumulative Distribution Function. All of its
#'  parameters need to be set, so that an element from the value range is the
#'  only argument to provided.
#'@param mean Mean of the target distribution.
#'@param var Variance of the target distribution.
#'@param gamma Defines a threshhold for the relative size of the subsummand in
#'  the computation of function \code{\link{expected_target_product}}

get_correlation_bound<- function(cdf, mean, var, gamma){
  return(c(get_target_correlation(expected_target_product = expected_target_product(cdf = cdf,r = -1, gamma = gamma), mean=mean, var=var), get_target_correlation(expected_target_product(cdf = cdf,r = 1, gamma = gamma), mean=mean, var=var)))
}

#' Build Gamma matrix for solving Yule-Walker-Equations.
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
  roots <- solve(polynom::polynomial(c(1, -alpha)))
  roots <- abs(roots)
  is_stationary <- TRUE
  if(sum(roots < 1)){
    is_stationary <- FALSE
  }
  return(is_stationary)
}


#' Create time-series of target distribution and autocorrelation structure, if
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
#'   target process autocorrelation space. Could theoretically be chosen arbitrarily, but should be unique
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process. Defaults to
#'   0.001.
#' @param gamma Threshold for computation of the expected value of the
#'   product of two random variables of the target process.
#' @seealso \code{\link{generate_distribution}}
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal distribution defined by distribution_name, or NULL if base
#'   process is not stationary.
#' @export
generate_DARTA <- function(n, cdf, inv, mean, var, rho, cdf_name_parameterized, epsilon, n_interpol, method, gamma, use_caching =T){
  # check if method arguments are valid
  match.arg(arg = method, choices =  c("interpol", "binary"))
  if(identical(method, "interpol")){
    if(length(n_interpol)>1){
      stop(paste("n_interpol should have length 1"))
    }
    if(!is.numeric(n_interpol) | n_interpol <= 0 ){
      stop(paste("When method 'interpol' is selected, n_interpol should be a positive natural number, but is ",n_interpol, sep = ""))
    }
  }
  if(identical(method, "binary")){
    if(length(epsilon)>1){
      stop(paste("epsilon should have length 1"))
    }
    if(epsilon <= 0 | epsilon >= 1){
      stop(paste("When method 'binary' is selected, epsilon should be a positive number much smaller than 1, but is ", epsilon, sep =""))
    }
  }
  # find fitting base process autocorrelation structure by specified method
  r <- switch(method,
              "interpol" = find_r_interpol(
                cdf = cdf,
                cdf_name_parameterized = cdf_name_parameterized,
                mean = mean,
                var = var,
                gamma = gamma,
                rho = rho,
                n_interpol = n_interpol,
                use_caching = use_caching
              ),
              "binary" = find_r_binary(
                cdf = cdf,
                mean = mean,
                var = var,
                rho = rho,
                cdf_name_parameterized = cdf_name_parameterized,
                epsilon = epsilon,
                gamma = gamma,
                use_caching = use_caching
              )[, 1])
  if(is.null(r)){
    stop("No suitable base process found")
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
    init_fun <- mvtnorm::rmvnorm
  }
  if(is_stationary(alpha)){
    # generate initial sample vector
    Z <- rev(as.vector(init_fun(1, rep(0, p),gamma)))
    # generate sample by inverse transform
    for(i in 1:n){
      Z[p+1] <-alpha %*% Z+ rnorm(1, mean=0,sd= sigma)
      Z <- Z[-1]
      Y[i] <- inv(pnorm(Z[p]))
    }
    return(Y)
  }else{
    stop("Base process is not stationary")
  }
}


