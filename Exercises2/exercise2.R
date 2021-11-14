# Exercise 2 (kidney dataframe from KMSurv)
library(KMsurv)
library(survival)
library(bshazard)
data(kidney)
?kidney
summary(kidney)
str(kidney)
head(kidney)

# a) Draw the estimated survival curves in both groups and estimate the results. 
surv1 <- with(kidney, Surv(time, delta)) # Make survival object. 
summary(surv1)
sfit1 <- survfit(surv1 ~ type, data = kidney) # Fit the survival curves (estimation with KM).
summary(sfit1)

plot(sfit1, col = 3:4, xlab = "Time to renal infection [Months]",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to catheter placement")
legend("bottomleft", legend = levels(as.factor(kidney$type)), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)

# b) Estimated median times. 
# sfit1 gives estimated median time for type 1 but not for type 2. Why? This is what the q is in task 1 also!

# c) bshazard for drawing hazard functions. Compare with the survival functions from earlier. 

