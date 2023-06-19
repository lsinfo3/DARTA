# This file contains functions representing the recommended way of using DARTA.
# Each method corresponds to a single marginal distribution, and creates a
# time-series whose values follow this marginal distribution, as well as the
# specified autocorrelation structure rho

#' Generate autocorrelated time-series with marginal negative binomial distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param size Size parameter of the Negative Binomial distribution.
#' @param prob Probability parameter of the Negative Binomial distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process when the 'binary' method is used. Defaults to
#'   0.001.
#' @param gamma Controls accuracy of DARTA's autocorrelation approximation mechanism, lower being more accurate. Defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @example man/examples/generate_nbinomial_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Negative Binomial marginal distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_nbinomial <- function(n, size, prob, rho, gamma = 0.00001, method = "interpol", n_interpol=20, epsilon = 0.001, use_caching = T){
  return(generate_distribution(n= n, rho = rho, distribution_name = "nbinomial", size = size, prob=prob, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon, use_caching = use_caching))
}

#' Generate autocorrelated time-series with marginal binomial distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param size Size parameter of the Binomial distribution.
#' @param prob Probability parameter of the Binomial distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process when the 'binary' method is used. Defaults to
#'   0.001.
#' @param gamma Controls accuracy of DARTA's autocorrelation approximation mechanism, lower being more accurate. Defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#'
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Binomial marginal distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_binomial <- function(n, size, prob, gamma = 0.00001, method = "interpol", n_interpol=20, rho, epsilon = 0.001, use_caching = T){
  return(generate_distribution(n= n, rho = rho, distribution_name = "binomial", size = size, prob=prob, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon, use_caching = use_caching))
}

#' Generate autocorrelated time-series with marginal poisson distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param lambda Lambda parameter of the Poisson distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process when the 'binary' method is used. Defaults to
#'   0.001.
#' @param gamma Controls accuracy of DARTA's autocorrelation approximation mechanism, lower being more accurate. Defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#'
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal Poisson-distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_poisson <- function(n, lambda, rho, gamma = 0.00001, method = "interpol", n_interpol=20, epsilon = 0.001, use_caching = T){
  return(generate_distribution(n= n, rho = rho, distribution_name = "poisson", lambda = lambda, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon, use_caching = use_caching))
}

#' Generate autocorrelated time-series with marginal uniform distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param min min parameter of the Binomial distribution.
#' @param max max parameter of the Binomial distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param epsilon Controls acceptable error within which the target
#'   autocorrelation is to be approximated by the base process when the 'binary' method is used. Defaults to
#'   0.001.
#' @param gamma Controls accuracy of DARTA's autocorrelation approximation mechanism, lower being more accurate. Defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#'
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Binomial marginal distribution.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_uniform <- function(n, min, max, rho, gamma = 0.00001, method = "interpol", n_interpol=20, epsilon = 0.001, use_caching = T){
  return(generate_distribution(n= n, rho = rho, distribution_name = "uniform", min = min, max=max, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon, use_caching = use_caching))
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
#' @param gamma Threshold for computation of the expected value of the
#'   product of two random variables of the target process. defaults to 10^-5.
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal distribution defined by distribution_name, or NULL if base
#'   process is not stationary.
#' @examples
#' generate_DARTA(n = 10000, rho = c(0.7,0.5,0.3), distribution_name = "nbinomial", param1 = 10, param2 = 0.6)
#' @seealso \code{\link{generate_DARTA}}
#' @export
generate_distribution <- function(n, rho, distribution_name = c("nbinomial", "binomial","poisson", "uniform"), gamma = 0.00001, method = c("interpol", "binary"), epsilon = 0.001, n_interpol = 20, use_caching = T, ...){
  distribution_name <- match.arg(distribution_name)
  params <- list(...)
  if(distribution_name == "nbinomial"){
    stopifnot("Please provide 'size' and 'prob' arguments for nbinomial distribution as additional function arguments" = ("size" %in% names(params) & "prob" %in% names(params)))
    size <- params$size
    prob <- params$prob
    stopifnot("Please provide valid size for nbinomial distribution" = is.numeric(size))
    stopifnot("Please provide valid prob for nbinomial distribution" = is.numeric(prob))
    mean = size*(1-prob)/prob
    var = size*(1-prob)/(prob**2)
    cdf_name_parameterized <- paste("nbinomial",as.character(size),as.character(prob),sep = "-")
    return(
      generate_DARTA(
        n = n,
        cdf = purrr::partial(pnbinom, size = size, prob = prob),
        inv = purrr::partial(qnbinom, size = size, prob = prob),
        mean = mean,
        var = var,
        cdf_name_parameterized = cdf_name_parameterized,
        rho = rho,
        epsilon = epsilon,
        n_interpol = n_interpol,
        method = method,
        gamma = gamma,
        use_caching = use_caching
      )
    )
  }
  if(distribution_name == "binomial"){
    stopifnot("Please provide 'size' and 'prob' arguments for binomial distribution as additional function arguments" = ("size" %in% names(params) & "prob" %in% names(params)))
    size <- params$size
    prob <- params$prob
    stopifnot("Please provide valid size for binomial distribution" = is.numeric(size))
    stopifnot("Please provide valid prob for binomial distribution" = is.numeric(prob))
    mean = size*prob
    var = size*prob*(1-prob)
    cdf_name_parameterized <-paste("binomial",as.character(size),as.character(prob),sep = "-")
    return(
      generate_DARTA(
        n = n,
        cdf = purrr::partial(pbinom, size = size, prob = prob),
        inv = purrr::partial(qbinom, size = size, prob = prob),
        cdf_name_parameterized = cdf_name_parameterized,
        mean = mean,
        var = var,
        rho = rho,
        epsilon = epsilon,
        n_interpol = n_interpol,
        method = method,
        gamma = gamma,
        use_caching = use_caching
      )
    )
  }
  if(distribution_name == "poisson"){
    stopifnot("Please provide 'lambda' argument for nbinomial distribution as additional function arguments" = "lambda" %in% names(params))
    lambda <- params$lambda
    stopifnot("please provide valid lambda for poisson distribution" = is.numeric(lambda))
    mean = lambda
    var = lambda
    cdf_name_parameterized <- paste("poisson",as.character(lambda),sep = "-")
    return(
      generate_DARTA(
        n = n,
        cdf = purrr::partial(ppois, lambda = lambda),
        inv = purrr::partial(qpois, lambda = lambda),
        cdf_name_parameterized = cdf_name_parameterized,
        mean = mean,
        var = var,
        rho = rho,
        epsilon = epsilon,
        n_interpol = n_interpol,
        method = method,
        gamma = gamma,
        use_caching = use_caching
      )
    )
  }
  if(distribution_name == "uniform"){
    stopifnot("Please provide 'min' and 'max' arguments for nbinomial distribution as additional function arguments" = ("min" %in% names(params) & "max" %in% names(params)))
    min <- params$min
    max <- params$max
    stopifnot("please provide valid min for uniform distribution" = is.numeric(min))
    stopifnot("please provide valid max for uniform distribution" = is.numeric(max))
    mean = (max+min)/2
    var = ((max- min)**2 -1)/12
    cdf_name_parameterized <- paste("uniform",as.character(min),as.character(max),sep = "-")
    return(
      generate_DARTA(
        n = n,
        cdf = purrr::partial(extraDistr::pdunif, min = min, max = max),
        inv = purrr::partial(extraDistr::qdunif, min = min, max = max),
        cdf_name_parameterized = cdf_name_parameterized,
        mean = mean,
        var = var,
        rho = rho,
        epsilon = epsilon,
        n_interpol = n_interpol,
        method = method,
        gamma = gamma,
        use_caching = use_caching
      )
    )
  }
}
