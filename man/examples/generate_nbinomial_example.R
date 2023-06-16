library(DARTA)
# define distribution parameters
size <- 10
prob <- 0.9

# define autocorrelation structure as a vector, starting with lag 1:
rho <- c(0.5,0.4,0.3)

# generate sample
vec <- generate_nbinomial(n = 1e5, size =size, prob=prob, rho=rho, method = "interpol", use_caching = F)

# plot autocorrelation
acf(vec)

# plot empirical CDF
plot(ecdf(vec))

# compare to negative bin distribution
plot(ecdf(rnbinom(1e5,size,prob)))

