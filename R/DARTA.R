

#'Computes the value of the integral for a plateau on the integrand function
#'surface. Plateau is specified by the values of the target distribution value
#'range that define its edges.
#'
#'@param j Element of target distribution's value range.
#'@param k Element of target distribution's value range.
#'@param cdf Target distribution's Cumulative Distribution Function. Needs to be fully parameterized.
#'@param r Autocorrelation between two variables of the base process. Should be
#'  selected such that it induces the targeted autocorrelation in the target
#'  process.
#' @seealso \code{\link{VGAM::pbinorm}}
integrand_plateau <- function(j, k, cdf, r) {
  val <- k * j * (VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(j)), cov12 = r) -
                    VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(j - 1)), cov12 = r) -
                    VGAM::pbinorm(qnorm(cdf(k - 1)), qnorm(cdf(j)), cov12 = r) +
                    VGAM::pbinorm(qnorm(cdf(k - 1)), qnorm(cdf(j - 1)), cov12 = r))
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
  (k^2) * (VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(k)), cov12 = r) -
             2 * VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(k-1)), cov12 = r) +
             VGAM::pbinorm(qnorm(cdf(k-1)), qnorm(cdf(k-1)), cov12 = r))
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
#'@param gamma Controls the quality of the integral approximation.
#'  Computation stops once the next subsum of integral plateaus being computed
#'  has a relative size compared to the current total sum smaller then the
#'  gamma.
#'
expected_target_product <- function(cdf,r, gamma){
  k <- 1
  total <- 0
  s <- 0
  s_max <- 0
  while(s >= gamma*total | k < 50){
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
      if(print_progress){
        print(paste0("r[",i,"]:", search_r, " || approximation for rho: ", current))
      }
      if(abs(current - rho[i])<epsilon){
        # approximation successful
        r_found = TRUE
        r[i] <- search_r
        approximation[i] <- current
      }else{
        interval <- if(rho[i]<current) c(interval[1], search_r) else c(search_r, interval[2])
      }
    }
  }
  return(data.frame(r, rho, approximation))
}

#' Determine Gamma matrix for solving Yule-Walker-Equations.
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
#' @example man/examples/'generate_nbinomial'_example.R
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and Negative Binomial marginal distribution, or NULL if base
#'   process is not stationary.
#' @export
#' @seealso \code{\link{generate_DARTA}}
generate_nbinomial <- function(n, size, prob, rho, gamma = 0.00001, method = "interpol", n_interpol=20, epsilon = 0.001){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "nbinomial", size = size, prob=prob, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon))
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
#' @seealso \code{\link{generate_DARTA}}
generate_binomial <- function(n, size, prob, gamma = 0.00001, method = "interpol", n_interpol=20, rho, epsilon = 0.001){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "binomial", size = size, prob=prob, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon))
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
#' @seealso \code{\link{generate_DARTA}}
generate_poisson <- function(n, lambda, rho, gamma = 0.00001, method = "interpol", n_interpol=20, epsilon = 0.001){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "poisson", lambda = lambda, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon))
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
#' @seealso \code{\link{generate_DARTA}}
generate_uniform <- function(n, min, max, rho, gamma = 0.00001, method = "interpol", n_interpol=20, epsilon = 0.001){
  return(generate_DARTA(n= n, rho = rho, distribution_name = "uniform", min = min, max=max, gamma = gamma, method = method, n_interpol=n_interpol, epsilon = epsilon))
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
#' @seealso \code{\link{gen_DARTA}}
#' @export
generate_DARTA <- function(n, rho, distribution_name = c("nbinomial", "binomial","poisson", "uniform"), gamma = 0.00001, method = c("interpol", "binary"), epsilon = 0.001, n_interpol = 20, use_caching = T, ...){
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
      gen_DARTA(
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
      gen_DARTA(
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
      gen_DARTA(
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
      gen_DARTA(
        n = n,
        cdf = purrr::partial(pdunif, min = min, max = max),
        inv = purrr::partial(qdunif, min = min, max = max),
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
#' @seealso \code{\link{generate_DARTA}}
#' @return Vector containing a time-series of length n with autocorrelation rho
#'   and marginal distribution defined by distribution_name, or NULL if base
#'   process is not stationary.
#' @export
gen_DARTA <- function(n, cdf, inv, mean, var, rho, cdf_name_parameterized, epsilon, n_interpol, method, gamma, use_caching =T){
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
  # find fitting base process autocorrelation structure
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
                gamma = gamma
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
