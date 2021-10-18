## =============================================================
## Lifetime Data Analysis, Course 2021/22
## 14.10.2021: Nonparametric estimation of the survival function
## =============================================================
setwd("/home/ajo/gitRepos/lifetime/lab2")
load("LdaLab2.RData")
ls()


## The data set ---------------------------------------------------------------
## ----------------------------------------------------------------------------
View(remission)
#
# 42 entries, 5 variables
# Labels with info about each variable
# Clinical Trial in order to compare a new treatment for leucemia with the standard treatment
# 21 pair of patients: 21 new + 21 old
# Time given in weeks, until relapse
# 5th variable not important today
# 
str(remission)
summary(remission)


## Kaplan-Meier estimation of S(t) with R -------------------------------------
## ----------------------------------------------------------------------------
library(survival)

## The Surv object
## ---------------
srem <- with(remission, Surv(time, cens))
srem
str(srem)
summary(srem)

## The Kaplan-Meier estimate is obtained with function survfit
## -----------------------------------------------------------
survfit(srem ~ 1)

# Estimation of the median but not of the mean, why?
# The censored observations make the estimation
# of the mean difficult, but not of the median.

## More information is shown with function summary.survfit
## -------------------------------------------------------
svf <- survfit(srem ~ 1)
summary(svf)

# Several options of function summary.survfit:
summary(svf, censored = TRUE)
summary(svf, times = 1:9 * 4)   # same as times = seq(4, 36, 4)
summary(svf, scale = 4) # Time in months.

# Why does function survfit not provide an estimation
# of the mean survival time?
# See comment above. 
survfit(srem ~ 1)

# A (biased) estimation of mean(T) can be obtained as follows
# (Restricted Mean Survival Time, RMST):
print(svf, rmean = "common")
print(svf, rmean = 50)


# Which information is saved in the survfit and summary.survfit objects?
# ----------------------------------------------------------------------
names(svf)
str(svf)
names(summary(svf))
str(summary(svf))

# Differences between both objects
with(svf, cbind(time, n.event, surv))
with(summary(svf), cbind(time, n.event, surv))


## Treatment comparison with respect to S(t)
## -----------------------------------------
(svf2t <- survfit(srem ~ treat, remission))
summary(svf2t)

##### New treatment quiet effective. 


## Computation of several quantiles
## --------------------------------
quantile(svf2t)


## Graphical representation of the estimated survival functions ---------------
## ----------------------------------------------------------------------------
plot(svf)

# A little bit nicer
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(svf, lwd = 3, xlab = "Weeks", ylab = expression(bolditalic(hat(S)(t))))
title("Survival function of the time to leukemia relapse")

# Both treatment arms
# -------------------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(svf2t, col = 1:2, conf.int = TRUE, lwd = 3, xlab = "Weeks", bty = "n",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i")
legend("topright", levels(remission$treat), col = 1:2, lwd = 3, bty = "n")

# Maybe better without the confidence intervals
plot(svf2t, col = 3:4, xlab = "Time to relapse [Weeks]",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to treatment")
legend("bottomleft", levels(remission$treat), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)

# Representation of the distribution function
# -------------------------------------------
plot(svf2t, col = 3:4, xlab = "Time to relapse [Weeks]",
     ylab = expression(bolditalic(hat(F)(t))), xaxs = "i",
     lty = 1:2, fun = "event", lwd = 3, yaxs = "i", bty = "n")
title("Distribution functions according to treatment")
legend("topleft", levels(remission$treat), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)

# Representation of the cumulative hazard function
# ------------------------------------------------
plot(svf2t, col = 3:4, xlab = "Time to relapse [Weeks]",
     ylab = expression(bolditalic(hat(Lambda)(t))), xaxs = "i",
     lty = 1:2, fun = "cumhaz", lwd = 3, yaxs = "i", bty = "l")
title("Cumulative hazards functions according to treatment")
legend("topleft", levels(remission$treat), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)


# Close all figures (if you want)
graphics.off()

# How to find help on that plot function?
? plot.survfit


## The Nelson-Aalen estimator of the survival function ------------------------
## ----------------------------------------------------------------------------
# Kaplan-Meier estimation
svf2t
summary(svf2t)

# Nelson-Aalen estimation
survfit(srem ~ treat, remission, type = "fleming")
# The same with more recent versions of the survival package
? survfit.formula
(svf2tNA <- survfit(srem ~ treat, remission, stype = 2, ctype = 1))
summary(svf2tNA)

# Graphical comparison:
# ---------------------
par(font = 2, font.axis = 2, font.lab = 2, las = 1, mar = c(5, 5, 4, 2))
plot(svf2t, col = 3:4, xlab = "Time to relapse [Weeks]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i",
     lty = 1, lwd = 3, yaxs = "i", bty = "l")
lines(svf2tNA, col = 3:4, lty = 2, lwd = 3, yaxs = "i", bty = "l")
title("Survival functions according to treatment")
legend("bottomleft", paste(rep(c("NA estimator:", "KM estimator:"), each = 2),
                           levels(remission$treat)), bty = "n", col = 3:4,
                           lty = c(2, 2, 1, 1), lwd = 3)
