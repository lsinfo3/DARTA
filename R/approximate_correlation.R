# This File contains functions pertaining to the computation of the
# autocorrelation in the target process, based on the autocorrelation in the base process.
# Specifically, if two random variables in the base process have autocorrelation r, and
# we aim for a target process with CDF cdf, then inserting the return value of
# the expected_target_product() function into the get_target_correlation()
# function yields the the autocorrelation function of the corresponding two
# random variables in the target process.

#'Compute the value of the integral for a plateau on the integrand function
#'surface.
#'
#'@param j Element of target distribution's value range.
#'@param k Element of target distribution's value range.
#'@param cdf Target distribution's Cumulative Distribution Function.
#'  Needs to be fully parameterized.
#'@param r Autocorrelation between two variables of the base process.
#'@seealso \code{\link{VGAM::pbinorm}}
integrand_plateau <- function(j, k, cdf, r) {
  val <- k * j * (VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(j)), cov12 = r) -
                    VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(j - 1)), cov12 = r) -
                    VGAM::pbinorm(qnorm(cdf(k - 1)), qnorm(cdf(j)), cov12 = r) +
                    VGAM::pbinorm(qnorm(cdf(k - 1)), qnorm(cdf(j - 1)), cov12 = r))
  val[val < 0] <- 0
  return(val)
}

#'Compute the value of the integral for a plateau on the integrand function
#'surface. Plateaus edges have to be defined by the same value from the target
#'distribution's value range.
#'
#'@param k Element of target distribution's value range.
#'@param cdf Target distribution's Cumulative Distribution Function.
#'  Needs to be fully parameterized.
#'@param r Autocorrelation between two variables of the base process.
#'@seealso \code{\link{VGAM::pbinorm}}
integrand_plateau_symmetric<- function(k, cdf, r) {
  val <-(k^2) * (VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(k)), cov12 = r) -
             2 * VGAM::pbinorm(qnorm(cdf(k)), qnorm(cdf(k-1)), cov12 = r) +
             VGAM::pbinorm(qnorm(cdf(k-1)), qnorm(cdf(k-1)), cov12 = r))
  val <- if(val < 0) 0 else val
  return(val)
}

#'Compute intermediate result for computation of autocorrelation between two
#'random variables in the target process.
#'@details
#'This function approximates the expected value of the product of two random
#'variables in the target process. By making use of the fact that the function
#'being approximated is piecewise steady, it sums up  the exact values of
#'integrals over plateaus, where a solution can be found analytically, until a
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
#'@returns Expected value of product of two random variables in the target
#'process, which are the product of the inverse transform of two random
#'variables in the base process with autocorrelation r
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
  }
  return(total)
}

#'Compute the autocorrelation between to random variables
#'of the target process based on the properties of the corresponding random
#'variables of the currently evaluated base process, i.e., the selected base
#'process autocorrelation and target process CDF.
#'@param expected_target_product Approximation (or exact value, if possible) of
#'  the integral function defining the expected value of two random variables of
#'  the target process. Usually return value of
#'  \link{expected_target_product}.
#'@param mean Mean of the target distribution.
#'@param var Variance of the target distribution.
#'@returns Autocorrelation between two variables of the target process.
get_target_correlation <- function(expected_target_product,mean, var){
  (expected_target_product - mean^2)/var
}
