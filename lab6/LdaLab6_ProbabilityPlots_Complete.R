## ======================================
## Lifetime Data Analysis, Course 2021/22
## 30.11.2021: Probability plots
## ======================================
library(survival)
set.seed(3011)

## Generation of Weibull-distributed survival times
stimes <- round(rweibull(1000, 2, 4), 2)
summary(stimes)

## Generation of exponentially-distributed censoring times
ctimes  <- round(rexp(1000, 0.1), 2)
summary(ctimes)

## Observed survival times and event indicator
obs <- pmin(stimes, ctimes)
cens <- as.numeric(stimes <= ctimes)
table(cens)

## Checking
head(cbind(stimes, ctimes, obs, cens))


## Probability plot for the Weibull and log-logistic distribution
## --------------------------------------------------------------
## Nelson-Aalen estimate of the survival function
svf <- survfit(Surv(obs, cens) ~ 1, stype = 2, ctype = 1)
svf
summary(svf)

## Uncensored survival times
times <- summary(svf)$time

## Nelson-Aalen estimate of the cumulative hazard function
chaz <- -log(summary(svf)$surv)

## The probability plots
windows(width = 14)
## quartz(width = 14)
par(mfrow = c(1, 2), las = 1, font.lab = 4, font.axis = 2, pch = 16)
plot(log(chaz) ~ log(times), xlab = "Log(Time)",
     ylab = "Log(Cumulative hazard)", col = 4)
abline(lm(log(chaz) ~ log(times)), lwd = 3)
title("Check for Weibull distribution")

plot(log(exp(chaz) - 1) ~ log(times), xlab = "Log(Time)",
     ylab = "Log(Exp(Cumulative hazard) - 1)", cex = 1.3, col = 4)
abline(lm(log(exp(chaz) - 1) ~ log(times)), lwd = 3)
title("Check for Log-logistic distribution")


## Probability plots with the GofCens package
## ------------------------------------------
library(GofCens)
packageDescription("GofCens")
windows(width = 14)
## quartz(width = 14)
par(las = 1, font.lab = 4, font.axis = 2, pch = 16)
cumhazPlot(obs, cens, col = 4)

windows(width = 14)
## quartz(width = 14)
par(las = 1, font.lab = 4, font.axis = 2, pch = 16)
cumhazPlot(obs, cens, col = 4, distr = c("wei", "loglo"), font.lab = 4)

windows(width = 14)
## quartz(width = 14)
cumhazPlot(obs, cens, col = 4, distr = c("wei", "loglo"), ggplo = TRUE)

windows(width = 12, height = 10)
## quartz(width = 14)
probPlot(obs, cens, distr = "weibull", col = 4)

windows(width = 12, height = 10)
## quartz(width = 14)
probPlot(obs, cens, distr = "loglo", col = 4, ggplo = TRUE)
