## =====================================================
## Lifetime Data Analysis, Course 2021/22
## 30.11.2021: Fit of the Accelerated Failure Time Model
## =====================================================
library(survival)
setwd("/home/ajo/gitRepos/lifetime/lab6")
load("LdaLab6.RData")

## The data set on larynx cancer
## =============================
str(larynx)
head(larynx, 10)
summary(larynx)

# The Surv object
# ---------------
larsurv <- with(larynx, Surv(time, cens))

## KM estimates of S(t) per disease stage
## --------------------------------------
windows(width = 10)
# quartz(width = 10)
par(font = 2, font.lab = 4, font.axis = 2, las = 1, mar = c(5, 5, 4, 2))
plot(survfit(larsurv ~ stagec, larynx), col = 1:4, xlab = "Years", lwd = 3,
     ylab = expression(bolditalic(hat(S)(t))), bty = "n")
legend("bottomleft", levels(larynx$stagec), col = 1:4, lwd = 3, bty = "n",
       title = "Disease stage")
title("Survival functions of patients with larynx cancer")


## Two Weibull regression models
## =============================
# Including disease stage
weimod2 <- survreg(larsurv ~ stagec, larynx)
summary(weimod2)
# Including disease stage and age
weimod3 <- update(weimod2, ~. + age)
summary(weimod3)


## Little exercise (Nov 23, 2021):
## · Compute the HRs associated with Stages 2, 3, and 4 with respect to Stage 1.
## · Compute the AFs associated with Stages 2, 3, and 4 with respect to Stage 1.
## · Which is the hazard ratio associated to an age difference of 5 years?
## -----------------------------------------------------------------------------
(HRs <- with(weimod3, exp(-coefficients[2:4] / scale)))
(AFs <- with(weimod3, exp(-coefficients[2:4])))
(HRage5 <- with(weimod3, exp(- 5 * coefficients[5] / scale)))


## Comparison with the log-logistic model
## --------------------------------------
loglomod3 <- update(weimod3, dist = "loglogistic")
summary(loglomod3)
# Which is the interpretation of exp(gamma/sigma) under this model?
with(loglomod3, round(exp(coefficients[2:4] / scale), 2))


## Check of the parametric model assumption
## ========================================
summary(weimod3)

# The residuals
# -------------
wm3pred <- predict(weimod3, type = "linear")
wm3pred
resids <- (log(larynx$time) - wm3pred) / weimod3$scale
resids

# The KM curve of the residuals
# -----------------------------
windows(width = 10)
par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(resids, larynx$cens) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the Weibull regression model")

# Graphical comparison with the theoretical survival function
survgumb <- function(x) {
  return(exp(-exp(x)))
}

curve(survgumb(x), from = min(resids), to = max(resids), col = 2, lwd = 3,
      add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Gumbel Distribution"),
       col = c(1, 1, 2), lty = c(1, 2, 1), lwd = 3, bty = "n")


## Exercise:
## Check the model fit with the lognormal distribution.
## Which parametric assumption seems to be more adequate?
## ------------------------------------------------------
lnomod3 <- update(weimod3, dist = "lognormal")
summary(lnomod3)

# The residuals
# -------------
lnopred <- predict(lnomod3, type = "linear")
residsLN <- (log(larynx$time) - lnopred) / lnomod3$scale
residsLN

# The KM curve of the residuals
# -----------------------------
windows(width = 10)
par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(residsLN, larynx$cens) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the lognormal regression model")

# Adding the theoretical survival function
curve(pnorm(x, lower.tail = FALSE), from = min(residsLN), to = max(residsLN),
      col = 2, lwd = 3, add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Normal Distribution"),
       col = c(1, 1, 2), lty = c(1, 2, 1), lwd = 3, bty = "n")


## Model-based prediction (using Model weimod3)
## ============================================
head(larynx)
predict(weimod3)[1:6]
## (a) Expected survival times
## ---------------------------
(newdat <- data.frame(age = rep(c(60, 65), each = 4),
                      stagec = paste("Stage", 1:4)))
predict(weimod3, newdata = newdat) * gamma(1 + weimod3$scale)


## (b) Quantile estimation
## -----------------------
predict(weimod3, type = "quantile", p = 1:3 / 4,
        newdata = newdat, se.fit = TRUE)
# A somewhat nicer presentation
quants <- cbind(newdat, predict(weimod3, type = "quantile", p = 1:3 / 4,
                newdata = newdat))
colnames(quants) <- c("Age", "Stage", "25%-Q.", "50%-Q.", "75%-Q.")
print(quants, digits = 2)


## (c) Estimation of survival probabilities
## ----------------------------------------
# install.packages("rms")
library(rms)
weimod3b <- psm(larsurv ~ stagec + age, larynx)
weimod3b
summary(weimod3b)
str(weimod3b)

survest(weimod3b, newdata = newdat, times = 2)
survest(weimod3, newdata = newdat, times = 2)     # Not very helpful ...
survest(weimod3b, newdat, times = 1:5)

# A somewhat nicer presentation
probs <- round(survest(weimod3b, newdat, times = 1:5), 3)
rownames(probs) <- paste(newdat$age, newdat$stagec, sep = " & ")
colnames(probs) <- c("1 Year", paste(2:5, "Years"))
{
  txt <- "Model-based prediction of survival probabilities"
  cat(txt, fill = TRUE)
  cat(rep("=", nchar(txt)), sep = "", fill = TRUE)
  probs
}

rm(txt, probs, quants, newdat)

## PD: In the case of left truncation, use function
##     aftreg (R package: eha)
