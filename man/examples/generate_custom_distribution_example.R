library(DARTA)
library(purrr)

# In this example, we are going to create autocorrelated values from a hypergeometric distribution.
# To do this, we need a fitting cumulative distribution function (cdf), quantile function (inv),
# and the distributions mean and variance

# Define Function parameters
n <-30
m <- 30
k <- 30

# Here, we pre-set the n parameter early, as it clashes with the n parameter in the generate_custom_distribution function
cdf <- partial(phyper, n = n)
inv <- partial(qhyper, n = n)

# Calculate mean and variance for the Hypergeometric distribution, see also https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html
mean <- k*m/(m+n)
var <- k*m/(m+n)*(1-m/(m+n))*(m+n-k)/(m+n-1)

# Atutocorrelation structure
rho <- c(0.5,0.3,0.2,0.15)

# generate sample.
# We pass the parameters 'm' and 'k' as function arguments. We could have also set them as we did with 'n' above, but when they are passed into the function, caching will be able to discern between different parameterizations of the same distribution
# The 'name' argument is only used when caching is enabled. Make sure to use different names for different custom distributions, or the wrong cache may be used.
vec <-
  generate_custom_distribution(
    n = 1e6,
    name = paste("hypergeometric", n, sep = "-"),
    cdf = cdf,
    inv = inv,
    mean = mean,
    var = var,
    m = m,
    k = k,
    rho = rho,
    method = "interpol",
    use_caching = T
  )

# plot autocorrelation
acf(vec)

# plot empirical CDF
plot(ecdf(vec))

# compare to binomial distribution
plot(ecdf(rhyper(1e5,n =n, m =m, k = k)))

