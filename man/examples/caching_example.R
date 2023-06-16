library(DARTA)

## Caching example

# define distribution parameters, autocorrelation structure
size <- 100
prob <- 0.25
rho <- c(0.5,0.4,0.3)

# generate sample, should take a few seconds

system.time(vec <- generate_nbinomial(n = 1e5, size =size, prob=prob, rho=rho, method = "interpol", use_caching = T))

# since values were cached, repeated execution should be quicker

system.time(vec <- generate_nbinomial(n = 1e5, size =size, prob=prob, rho=rho, method = "interpol", use_caching = T))

# A subdirectory named '.interpol_caches' with the corresponding cache can be found in the current working directory.
# There is also a separate '.binary_caches' directory containing caches for the 'binary' method
