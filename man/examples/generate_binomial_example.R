library(DARTA)
# define distribution parameters
size <- 30
prob <- 0.7

# define autocorrelation structure as a vector, starting with lag 1:
rho <- c(0.7,0.4,0.1)

# generate sample
vec <- generate_binomial(n = 1e6, size =size, prob=prob, rho=rho, method = "interpol", use_caching = F)

# plot autocorrelation
acf(vec)

# plot empirical CDF
plot(ecdf(vec))

# compare to binomial distribution
plot(ecdf(rbinom(1e5,size,prob)))

