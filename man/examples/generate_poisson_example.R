library(DARTA)
# define distribution parameter
lambda <- 1.5

# define autocorrelation structure as a vector, starting with lag 1:
rho <- c(0.9)

# generate sample
vec <- generate_poisson(n = 1e5, lambda = lambda, rho=rho, method = "interpol", use_caching = F)

# plot autocorrelation
acf(vec)

# plot empirical CDF
plot(ecdf(vec))

# compare to poisson distribution
plot(ecdf(rpois(1e5,lambda)))
