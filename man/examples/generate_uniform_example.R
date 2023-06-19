library(DARTA)
# define distribution parameters
min <- 50
max <- 150

# define autocorrelation structure as a vector, starting with lag 1:
rho <- c(0.3,0.2,0.3,0.2,0.3,0.2)

# generate sample
vec <- generate_uniform(n = 1e5, min = min, max = max, rho=rho, method = "interpol", use_caching = F)

# plot autocorrelation
acf(vec)

# plot empirical CDF
plot(ecdf(vec))

# compare to uniform distribution
plot(ecdf(runif(1e5,min, max)))

