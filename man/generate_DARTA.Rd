% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DARTA.R
\name{generate_DARTA}
\alias{generate_DARTA}
\title{Create time-series of target distribution and autocorrelation structure, if
possible.}
\usage{
generate_DARTA(
  n,
  cdf,
  inv,
  mean,
  var,
  rho,
  cdf_name_parameterized,
  epsilon,
  n_interpol,
  method,
  gamma,
  use_caching = T
)
}
\arguments{
\item{n}{length of time-series to be generated, i.e., total size of the
vector of random numbers that is returned.}

\item{cdf}{Marginal CDF of target process.}

\item{inv}{Inverse of marginal CDF of target process, a.k.a. quantile
function.}

\item{mean}{Mean of marginal CDF of target process.}

\item{var}{Variance of marginal CDF of target process.}

\item{rho}{Target autocorrelation Structure of target process, provided as a
vector of numerics. Position on the vector describes the lag for which the
numeric autocorrelation value is supposed to occur.}

\item{cdf_name_parameterized}{Name of the target Distribution, together with
its parameters, separated by dashes. Used for naming the cache-file, which
stores a map of values from the base process autocorrelation space to the
target process autocorrelation space. Could theoretically be chosen arbitrarily, but should be unique}

\item{epsilon}{Controls acceptable error within which the target
autocorrelation is to be approximated by the base process. Defaults to
0.001.}

\item{gamma}{Threshold for computation of the expected value of the
product of two random variables of the target process.}
}
\value{
Vector containing a time-series of length n with autocorrelation rho
and marginal distribution defined by distribution_name, or NULL if base
process is not stationary.
}
\description{
Create time-series of target distribution and autocorrelation structure, if
possible.
}
\details{
It is generally encouraged to use \code{\link{generate_DARTA}}
instead, which calls upon this function.
}
\seealso{
\code{\link{generate_distribution}}
}
