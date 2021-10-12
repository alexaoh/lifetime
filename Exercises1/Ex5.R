# Problem 5 d) Estimate the parameters of the Weibull distribution.

library(fitdistrplus)
df <- data.frame(left = c(11, 12, 15, 33, 45, 28, 16, 17, 19, 30),
                right = c(11, 12, 15, NA, 45, NA, NA, NA, NA, NA))
obj <- fitdistcens(censdata = df, distr = "weibull")
k <- obj$estimate[1]
k
rho <- 1/obj$estimate[2]
rho
plot(obj)
