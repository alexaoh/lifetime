# a) Survival during 40 and 60 days.

1-pweibull(40, shape = 2, scale = 1/0.03)
1-pweibull(60, shape = 2, scale = 1/0.03)

# b) 80\% quantile. 
qweibull(0.8, shape = 2, scale = 1/0.03)
