# Exercise 2 (kidney dataframe from KMSurv)
library(dplyr)
library(KMsurv)
library(survival)
library(bshazard)
library(km.ci)
data(kidney)
?kidney
summary(kidney)
str(kidney)
head(kidney)

# a) Draw the estimated survival curves in both groups and estimate the results. 
surv1 <- with(kidney, Surv(time, delta)) # Make survival object. 
summary(surv1)
sfit1 <- survfit(surv1 ~ type, data = kidney) # Fit the survival curves (estimation with KM).
sfit1
summary(sfit1)

plot(sfit1, col = 3:4, xlab = "Time to renal infection [Months]",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to catheter placement")
legend("bottomleft", legend = levels(as.factor(kidney$type)), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)

# b) Estimated median times. 
# sfit1 gives estimated median time for type 1 but not for type 2. Why? This is what the q is in task 1 also!

# Biased estimate of the median. Is this what they want here?
# DOES NOT WORK FOR MEDIAN!! FIND HOW I CAN ESTIMATE THE MEDIAN!
print(sfit1, rmedian = "common")
print(sfit1, rmedian = 30) 

# c) bshazard for drawing hazard functions. Compare with the survival functions from earlier. 
hazard1 <- bshazard(with(kidney[kidney$type == 1,], Surv(time, delta)) ~ 1, data = kidney[kidney$type == 1,])
plot(hazard1)
hazard2 <- bshazard(with(kidney[kidney$type == 2,], Surv(time, delta)) ~ 1, data = kidney[kidney$type == 2,])
plot(hazard2)

# Kan jeg beregne Hazardene med Nelson-Aalen for å se om disse gir mening!?

# d) Estimated probabilities of not suffering renal infection after 6, 12, 18, 24 and 30 months?
surv.group1 <- with(kidney[kidney$type == 1,], Surv(time, delta))
sfit.group1 <- survfit(surv.group1 ~ type, data = kidney[kidney$type == 1,])
plot(sfit.group1)
abline(v=c(6,12,18,24,30))
summary(sfit.group1) # Kanskje best å bare lese av fra summaryen!
(surv.times.df <- data.frame(sfit.group1$time, sfit.group1$surv, sfit.group1$lower, sfit.group1$upper))
surv.times.df[c(5.5, 11.5, 16.5, 23.5, 26.5), ] # 23.5 does not work for some reason, but I can see it.

surv.group2 <- with(kidney[kidney$type == 2,], Surv(time, delta))
sfit.group2 <- survfit(surv.group2 ~ type, data = kidney[kidney$type == 2,])
plot(sfit.group2)
abline(v=c(6,12,18,24,30))
summary(sfit.group2) # Kanskje best å bare lese av fra summaryen!
(surv.times.df.2 <- data.frame(sfit.group2$time, sfit.group2$surv, sfit.group2$lower, sfit.group2$upper))
surv.times.df.2[c(5.5, 11.5, 16.5, 24.5, 28.5, 15.5), ] # Are some NAs because some of them dont exist in summary I think. 

# e) Lower and upper limits of linear EP confidence bands with level 95% at the same levels as d).
CI.lin.group1 <- km.ci(sfit.group1, conf.level = 0.95, method = "linear")
summary(CI.lin.group1) # Can see that the confidence intervals are modified!
CI.lin.group2 <- km.ci(sfit.group2, conf.level = 0.95, method = "linear")
summary(CI.lin.group2) # Can see that the confidence intervals are modified!
