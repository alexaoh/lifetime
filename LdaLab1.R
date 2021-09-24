## ===================================================
## Lifetime Data Analysis, Course 2021/22
## 23.9.2021: Introduction to Survival Analysis with R
## ===================================================

## A couple of survival functions
## ==============================
mu <- 3; si <- 2
al <- 0.5; la <- 55
windows(width = 15)
par(mfrow = 1:2, las = 1, font = 2, font.axis = 2, font.lab = 4,
    mar = c(5, 4, 1, 2), bty = "n", xaxs = "i", yaxs = "i")
curve(plnorm(x, mu, si, lower.tail = FALSE), from = 0, to = 100, lwd = 3,
      xlab = "Time", ylab = "S(t)", ylim = 0:1)
axis(1, seq(10, 90, 20))
axis(2, seq(0.1, 0.9, 0.1))
curve(pweibull(x, al, la, lower.tail = FALSE), from = 0, to = 100, lwd = 3,
      xlab = "Time", ylab = "S(t)", ylim = 0:1)
axis(1, seq(10, 90, 20))
axis(2, seq(0.1, 0.9, 0.1))

## What shape do the corresponding hazard functions have?
## ------------------------------------------------------
library(eha)
windows(width = 15)
par(mfrow = 1:2, las = 1, font = 2, font.axis = 2, font.lab = 4,
    mar = c(5, 4, 1, 2), bty = "n", xaxs = "i", yaxs = "i")
curve(hlnorm(x, mu, si), from = 0, to = 100, lwd = 3, ylim = c(0, 0.07),
      xlab = "Time", ylab = expression(lambda(t)))
axis(1, seq(10, 90, 20))
curve(hweibull(x, al, la), from = 0, to = 100, lwd = 3, ylim = c(0, 0.07),
      xlab = "Time", ylab = expression(lambda(t)))
axis(1, seq(10, 90, 20))

windows(width = 15)
par(mfrow = 1:2, las = 1, font = 2, font.axis = 2, font.lab = 4,
    mar = c(5, 4, 1, 2), bty = "n", xaxs = "i", yaxs = "i")
curve(hlnorm(x, mu, si), from = 0, to = 10, lwd = 3, ylim = c(0, 0.08),
      xlab = "Time", ylab = expression(lambda(t)))
axis(1, seq(1, 9, 2))
curve(hweibull(x, al, la), from = 0, to = 10, lwd = 3, ylim = c(0, 0.08),
      xlab = "Time", ylab = expression(lambda(t)))
axis(1, seq(1, 9, 2))


## =================================
## "Illustration DMBA" (Slide 45/89)
## =================================
ratcancer <- data.frame(group = factor(rep(c("Group 1", "Group 2"), c(19, 21))),
                        times = c(143, 164, 188, 188, 190, 192, 206, 209, 213,
                                  216, 220, 227, 230, 234, 246, 265, 304, 216,
                                  244, 142, 156, 163, 198, 205, 232, 232, 233,
                                  233, 233, 233, 239, 240, 261, 280, 280, 296,
                                  296, 323, 204, 344),
                        cens = rep(c(1, 0, 1, 0), c(17, 2, 19, 2)))
summary(ratcancer)

## Empirical distribution function (ignoring right censoring)
## ----------------------------------------------------------
edfG1 <- with(subset(ratcancer, group == "Group 1"), ecdf(times)) # ecdf: Empirical Cumulative Distribution Function. 
edfG1
edfG1(seq(100, 350, 10))
maxtim <- max(ratcancer$times)
(srvedfG1 <- round(1 - edfG1(0:maxtim), 3))
edfG2 <- with(subset(ratcancer, group == "Group 2"), ecdf(times))
srvedfG2 <- round(1 - edfG2(0:maxtim), 3)

## Graphical representation of the survival function
## -------------------------------------------------
par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l", mar = c(5, 5, 4, 2))
#par(mfrow = c(1,1))
plot(0:maxtim, srvedfG1, xlab = "Survival time [Days]", xaxs = "i", yaxs = "i",
     ylab = expression(bolditalic(hat(S)(t))), type = "s", lwd = 3,
     xlim = c(0, 350))
points(0:maxtim, srvedfG2, type = "s", lwd = 3, col = 2)
axis(1, seq(25, 325, 50))
legend("bottomleft", levels(ratcancer$group), col = 1:2, lwd = 3, bty = "n")
title("Estimated survival function (ignoring censoring)")


## Estimation of the survival function with right-censored data in R
## =================================================================
library(survival)

## The Surv object with right-censored data
## ----------------------------------------
with(ratcancer, Surv(times, cens))

## Estimation of the survival function in both groups
## --------------------------------------------------
svf <- survfit(Surv(times, cens) ~ group, ratcancer)
svf
summary(svf)

## Graphical representation
## ------------------------
par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l",
    mar = c(5, 5, 4, 2))
plot(svf, col = 1:2, lwd = 3, xlab = "Survival time [Days]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 350))
axis(1, seq(25, 325, 50))
legend("bottomleft", levels(ratcancer$group), lwd = 3, col = 1:2, bty = "n")
title("Survival functions of time to death due to vaginal cancer")


## We can superpose the estimates that ignored censoring
points(0:maxtim, srvedfG1, type = "s", lty = 2, lwd = 3, col = 1)
points(0:maxtim, srvedfG2, type = "s", lty = 2, lwd = 3, col = 2)
legend("left", paste("Group", 1:2, "(ignoring right censoring)"), lwd = 3,
       col = 1:2, bty = "n", lty = 2)


## Which parametric model might be the best choice?
## ================================================
library(GofCens)
with(subset(ratcancer, group == "Group 1"),
     cumhazPlot(times, cens, distr = c("weibull", "loglog", "lognormal")))

with(subset(ratcancer, group == "Group 2"),
     cumhazPlot(times, cens, distr = c("weibull", "loglog", "lognormal"),
                ggplo = TRUE))
