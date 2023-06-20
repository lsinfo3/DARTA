# This file contains functions representing the recommended way of using DARTA.
# Each method corresponds to a single marginal distribution, and creates a
# time-series whose values follow this marginal distribution, as well as the
# specified autocorrelation structure rho

#' Generate autocorrelated time-series with marginal negative binomial
#' distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param size Size parameter of the Negative Binomial distribution.
#' @param prob Probability parameter of the Negative Binomial distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param gamma Numeric, controls accuracy of DARTA's autocorrelation
#'   approximation mechanism, lower being more accurate. Defaults to 10^-5.
#'   product of two random variables of the target process. defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @param n_interpol Natural number, decides how many equidistant points in the
#'   base process autocorrelation space should be sampled in order to fit a
#'   polynomial function. Defaults to 30.
#' @param epsilon Numeric between 0 and 1, controls acceptable error within
#'   which the target autocorrelation is to be approximated by the base process
#'   when the 'binary' method is used. Defaults to 0.001.
#' @param use_caching Logical, deciding wether results should be cached, and
#'   cached results used in the computation.
#' @example man/examples/generate_nbinomial_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and negative binomial marginal distribution, or NULL if base process is not
#'   stationary.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_nbinomial <- function(n, size, prob, rho, gamma = 0.00001, method = "interpol", n_interpol=30, epsilon = 0.001, use_caching = T){
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
#' @param gamma Numeric, controls accuracy of DARTA's autocorrelation
#'   approximation mechanism, lower being more accurate. Defaults to 10^-5.
#'   product of two random variables of the target process. defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @param n_interpol Natural number, decides how many equidistant points in the
#'   base process autocorrelation space should be sampled in order to fit a
#'   polynomial function. Defaults to 30.
#' @param epsilon Numeric between 0 and 1, controls acceptable error within
#'   which the target autocorrelation is to be approximated by the base process
#'   when the 'binary' method is used. Defaults to 0.001.
#' @param use_caching Logical, deciding wether results should be cached, and
#'   cached results used in the computation. Defaults to TRUE.
#' @example man/examples/generate_binomial_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Binomial marginal distribution, or NULL if base process is not
#'   stationary.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_binomial <- function(n, size, prob, rho, gamma = 0.00001, method = "interpol", n_interpol=30, epsilon = 0.001, use_caching = T){
  return(generate_distribution(n= n, rho = rho, distribution_name = "binomial", size = size, prob=prob, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon, use_caching = use_caching))
}

#' Generate autocorrelated time-series with marginal poisson distribution.
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param lambda Lambda parameter of the Poisson distribution.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param gamma Numeric, controls accuracy of DARTA's autocorrelation
#'   approximation mechanism, lower being more accurate. Defaults to 10^-5.
#'   product of two random variables of the target process. defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @param n_interpol Natural number, decides how many equidistant points in the
#'   base process autocorrelation space should be sampled in order to fit a
#'   polynomial function. Defaults to 30.
#' @param epsilon Numeric between 0 and 1, controls acceptable error within
#'   which the target autocorrelation is to be approximated by the base process
#'   when the 'binary' method is used. Defaults to 0.001.
#' @param use_caching Logical, deciding wether results should be cached, and
#'   cached results used in the computation. Defaults to TRUE.
#' @example man/examples/generate_poisson_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal poisson distribution, or NULL if base process is not
#'   stationary.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_poisson <- function(n, lambda, rho, gamma = 0.00001, method = "interpol", n_interpol=30, epsilon = 0.001, use_caching = T){
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
#' @param gamma Numeric, controls accuracy of DARTA's autocorrelation
#'   approximation mechanism, lower being more accurate. Defaults to 10^-5.
#'   product of two random variables of the target process. defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @param n_interpol Natural number, decides how many equidistant points in the
#'   base process autocorrelation space should be sampled in order to fit a
#'   polynomial function. Defaults to 30.
#' @param epsilon Numeric between 0 and 1, controls acceptable error within
#'   which the target autocorrelation is to be approximated by the base process
#'   when the 'binary' method is used. Defaults to 0.001.
#' @param use_caching Logical, deciding wether results should be cached, and
#'   cached results used in the computation. Defaults to TRUE.
#' @example man/examples/generate_uniform_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and uniform marginal distribution.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_uniform <- function(n, min, max, rho, gamma = 0.00001, method = "interpol", n_interpol=30, epsilon = 0.001, use_caching = T){
  return(generate_distribution(n= n, rho = rho, distribution_name = "uniform", min = min, max=max, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon, use_caching = use_caching))
}

#' Generate time-series of a custom distribution
#'
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param name Name of the custom distribution, used for caching. Should be chosen uniquely for each cdf, but different parameterizations can be cached for the same name.
#' @param cdf Cumulative distribution function of the custom distribion
#' @param inv Inverse cumulative distribution function, i.e. qunatile function, of cumstom distribution
#' @param mean Mean of Custom distribution
#' @param var Variance of custom distribution
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   \code{numerics}.
#' @param gamma Numeric, controls accuracy of DARTA's autocorrelation
#'   approximation mechanism, lower being more accurate. Defaults to 10^-5.
#'   product of two random variables of the target process. defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated fdevrom equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @param n_interpol Natural number, decides how many equidistant points in the
#'   base process autocorrelation space should be sampled in order to fit a
#'   polynomial function. Defaults to 30.
#' @param epsilon Numeric between 0 and 1, controls acceptable error within
#'   which the target autocorrelation is to be approximated by the base process
#'   when the 'binary' method is used. Defaults to 0.001.
#' @param use_caching Logical, deciding whether results should be cached, and
#'   cached results used in the computation. Defaults to TRUE.
#' @example man/examples/generate_uniform_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and cumstom marginal distribution.
#' @export
#' @seealso \code{\link{generate_distribution}}
generate_custom_distribution <-
  function(n, cdf, inv, mean, name,var, rho, gamma = 0.00001, method = "interpol", n_interpol = 30, epsilon = 0.001, use_caching = T,  ...) {
    return(
      generate_distribution(
        n = n,
        rho = rho,
        distribution_name = "custom",
        cdf = cdf,
        inv = inv,
        name = name,
        mean = mean,
        var = var,
        gamma = gamma,
        method = method,
        n_interpol = n_interpol,
        epsilon = epsilon,
        use_caching = use_caching,
        ... = ...
      )
    )
  }

#' Generate time-series of the provided marginal distribution and
#' autocorrelation structure.
#' @param n length of time-series to be generated, i.e., total size of the
#'   vector of random numbers that is returned.
#' @param rho Autocorrelation structure to be approximated, as a vector of
#'   numerics.
#' @param distribution_name Name of the Distribution to be generated. Currently
#'   supported are "nbinomial", "binomial", "poisson", and "uniform".
#' @param gamma Numeric, controls accuracy of DARTA's autocorrelation
#'   approximation mechanism, lower being more accurate. Defaults to 10^-5.
#'   product of two random variables of the target process. defaults to 10^-5.
#' @param method Defines method for autocorrelation structure approximation.
#'   'binary' employs a binary search, ending when autocorrelation in target
#'   process is within error margin 'epsilon' of desired autocorrelation.
#'   'interpol' fits a function through 'n_interpol' autocorrelation values
#'   calculated from equidistant points in the base process autocorrelation
#'   space, then interpolates a fitting base process autocorrelation from these
#'   values.
#' @param n_interpol Natural number, decides how many equidistant points in the
#'   base process autocorrelation space should be sampled in order to fit a
#'   polynomial function. Defaults to 30.
#' @param epsilon Numeric between 0 and 1, controls acceptable error within
#'   which the target autocorrelation is to be approximated by the base process
#'   when the 'binary' method is used. Defaults to 0.001.
#' @param use_caching Logical, deciding wether results should be cached, and
#'   cached results used in the computation.
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal distribution defined by distribution_name, or NULL if base
#'   process is not stationary.
#' @examples
#' generate_DARTA(n = 10000, rho = c(0.5,0.4,0.3), distribution_name = "nbinomial", size = 10, prob = 0.6, method = "interpol")
#' @seealso \code{\link{generate_DARTA}}
#' @export
generate_distribution <- function(n, rho, distribution_name = c("nbinomial", "binomial","poisson", "uniform", "custom"), gamma = 0.00001, method = c("interpol", "binary"), epsilon = 0.001, n_interpol = 30, use_caching = T, ...){
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
    stopifnot("Please provide valid lambda for poisson distribution" = is.numeric(lambda))
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
    stopifnot("Please provide 'min' and 'max' arguments for nbinomial distribution as additional function arguments" = ("cdf" %in% names(params) & "cdf" %in% names(params)))
    min <- params$min
    max <- params$max
    stopifnot("Please provide valid min for uniform distribution" = is.numeric(min))
    stopifnot("Please provide valid max for uniform distribution" = is.numeric(max))
    stopifnot("Please make sure min < max" = min < max)
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
  if(distribution_name == "custom"){
    if(use_caching){
      warning("Be careful to provide distinct custom names for different cdf/inv when using caching.")
      stopifnot("Please provide name, or disable caching" = !("name" %in% params))
      name <- params$name
    }
    stopifnot("Please provide cdf" = ("cdf" %in% names(params)))
    stopifnot("Please provide inv" = ("inv" %in% names(params)))
    stopifnot("Please provide mean" = ("mean" %in% names(params)))
    stopifnot("Please provide var" = ("var" %in% names(params)))
    mean <- params$mean
    var <- params$var
    cdf <- params$cdf
    inv <- params$inv

    # filter non distribution params
    dstr_params <- params[names(params) %in% c("cdf", "inv", "mean", "var", "name") == FALSE]
    cdf_name_parameterized <- if(use_caching) paste(name, paste(dstr_params, collapse = "-"), sep ="_") else "custom"

    if(length(dstr_params)>0){
      # set params for distribution functions
      cdf <- do.call(partial,c(list(.f = cdf), dstr_params))
      inv <- do.call(partial,c(list(.f = inv), dstr_params))
    }
    return(
      generate_DARTA(
        n = n,
        cdf = cdf,
        inv = inv,
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
