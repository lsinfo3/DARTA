source("./DARTA.R")
#generate sample
vec <- generate_nbinomial(n = 1e5, size =20, prob=0.9, rho=c(0.5,0.4,0.3), epsilon = 0.001, precision = 0.00001)
# plot autocorrelation
acf(vec)

#plot empirical CDF
plot(ecdf(vec))

# compare to dist
plot(ecdf(rnbinom(1e5,size,prob)))
