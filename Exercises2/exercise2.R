# Exercise 2 (kidney dataframe from KMSurv)
library(dplyr)
library(KMsurv)
library(survival)
library(bshazard)
library(km.ci)
data(kidney)
#?kidney
summary(kidney)
str(kidney)
head(kidney)

kidney$type <- factor(kidney$type, labels = c("Surgically", "Percutaneously"))

# a) Draw the estimated survival curves in both groups and estimate the results. 
surv1 <- with(kidney, Surv(time, delta)) # Make survival object. 
summary(surv1)
sfit1 <- survfit(surv1 ~ type, data = kidney) # Fit the survival curves (estimation with KM).
sfit1
summary(sfit1)

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
#par(font = 2, font.axis = 2, font.lab = 4, las = 1)
plot(sfit1, col = 3:4, xlab = "Time to renal infection [Months]",
     ylab = expression(bold(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to catheter placement")
legend("bottomleft", legend = levels(kidney$type), title = "Catheter Placement",
      bty = "n", col = 3:4, lty = 1:2, lwd = 3)
axis(1, at = seq(0, 30, 5))
axis(2, at = seq(0, 1, 0.1))


# b) Estimated median times. 
# sfit1 gives estimated median time for type 1 but not for type 2. This is because Type 2 never crosses the median quantile, 
# so the median cannot be estimated for the Percutaneosly placed catheter. 

# c) bshazard for drawing hazard functions. Compare with the survival functions from earlier. 
hazard1 <- bshazard(with(kidney[kidney$type == "Surgically",], Surv(time, delta)) ~ 1, 
                    data = kidney[kidney$type == "Surgically",])

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(hazard1, main = "Smoothed estimate of hazard function for surgical placement", ylab = expression(bold(hat(lambda)(t))), 
     xlab = "Time to renal infection [Months]")

hazard2 <- bshazard(with(kidney[kidney$type == "Percutaneously",], Surv(time, delta)) ~ 1, 
                    data = kidney[kidney$type == "Percutaneously",])

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(hazard2, main = "Smoothed estimate of hazard function for percutaneous placement", ylab = expression(bold(hat(lambda)(t))), 
     xlab = "Time to renal infection [Months]")

# Kan jeg beregne Hazardene med Nelson-Aalen for å se om disse gir mening!?

# d) Estimated probabilities of not suffering renal infection after 6, 12, 18, 24 and 30 months?
surv.group1 <- with(kidney[kidney$type == "Surgically",], Surv(time, delta))
sfit.group1 <- survfit(surv.group1 ~ type, data = kidney[kidney$type == "Surgically",])
plot(sfit.group1)
abline(v=c(6,12,18,24,30))
summary(sfit.group1) # Kanskje best å bare lese av fra summaryen!
(surv.times.df <- data.frame(sfit.group1$time, sfit.group1$surv, sfit.group1$lower, sfit.group1$upper))
surv.times.df[c(5.5, 11.5, 16.5, 23.5, 26.5), ] # 23.5 does not work for some reason, but I can see it.

surv.group2 <- with(kidney[kidney$type == "Percutaneously",], Surv(time, delta))
sfit.group2 <- survfit(surv.group2 ~ type, data = kidney[kidney$type == "Percutaneously",])
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
