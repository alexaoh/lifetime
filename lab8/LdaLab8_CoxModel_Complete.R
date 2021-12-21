## =====================================================
## Lifetime Data Analysis, Course 2021/22
## 21.12.2021: Fit of the Cox Proportional Hazards model
## =====================================================
setwd("/home/ajo/gitRepos/lifetime/lab8")
load("LdaLab8.RData")
library(survival)


## A data set on patients with multiple myeloma
## ============================================
str(myeloma)
head(myeloma)
summary(myeloma)
Hmisc::hist.data.frame(myeloma, nclass = 10)

## The Surv object
## ---------------
myesurv <- with(myeloma, Surv(time, status))

## Two survival plots
## ------------------
windows(width = 14)
par(mfrow = c(1, 2), font = 2, font.lab = 4, font.axis = 2, las = 1,
    oma = c(0, 0, 1, 0), mar = c(5, 5, 4, 2))
plot(survfit(myesurv ~ sex, myeloma), xlab = "Months", col = 1:2,
     lwd = 3, ylab = expression(bolditalic(hat(S)(t))), bty = "l", xaxs = "i")
title("Survival functions according to gender")
legend("bottomleft", levels(myeloma$sex), col = 1:2, lwd = 3, bty = "n")
plot(survfit(myesurv ~ BJP, myeloma), xlab = "Months", col = 1:2,
     lwd = 3, ylab = expression(bolditalic(hat(S)(t))), bty = "l", xaxs = "i")
title("Survival functions according to Bence-Jones protein")
legend("bottomleft", levels(myeloma$BJP), col = 1:2, lwd = 3, bty = "n")


## The fit of a Cox model
## ======================
(cox1 <- coxph(myesurv ~ sex + BJP, myeloma))
summary(cox1)

## Including the interaction between both variables
## ------------------------------------------------
(cox1 <- update(cox1,  ~ . + sex:BJP))
summary(cox1)


## The model-based hazard ratios associated with the absence of the BJ protein
## ---------------------------------------------------------------------------
library(Epi)
ci.lin(cox1)
round(ci.lin(cox1, Exp = TRUE), 3)
(ctmat <- matrix(c(0, 1, 0, 0, 1, 1), byrow = TRUE, nr = 2))
round(ci.lin(cox1, ctr.mat = ctmat, Exp = TRUE), 3)

# A somewhat nicer presentation
HRmat <- round(ci.lin(cox1, ctr.mat = ctmat, Exp = TRUE), 3)[, c(1, 5:7)]
rownames(HRmat) <- c("HR| Women", "HR| Men")
colnames(HRmat) <- c("logHR", "HR", "Lower 95%", "Upper 95%")
HRmat


## The Cox model including, also, BUN, HB, and age
## ===============================================
cox2 <- update(cox1,  ~ . + BUN + HB + age - sex:BJP, myeloma)
summary(cox2)

## Which is the hazard ratio associated with an increase of 10 units in BUN?
## Which is the corresponding 95% confidence interval?
## -------------------------------------------------------------------------
round(ci.lin(cox2, ctr.mat = matrix(c(0, 0, 10, 0, 0), nr = 1), Exp = TRUE), 3)


## Estimation of the baseline hazard function
## ------------------------------------------
coxph.detail(cox2)
coxph.detail(cox2)$hazard


## Model-based prediction
## ======================
survfit(cox2)
plot(survfit(cox2), lwd = 2, col = 3)

(newdat <- data.frame(BUN = 25, sex = factor(rep(c("Female", "Male"), 2)),
                     HB = 10, age = 65,
                     BJP = factor(rep(c("No", "Yes"), each = 2))))

survfit(cox2, newdata = newdat)
summary(survfit(cox2, newdata = newdat))

## Graphical respresentation
## -------------------------
windows(width = 8)
par(font = 2, font.lab = 4, font.axis = 2, las = 1)
plot(survfit(cox2, newdata = newdat), col = 1:4, lwd = 3, xlab = "Time [Months]",
     bty = "n", xlim = c(0, 100))
legend("bottomleft", c("Woman with BJP", "Man with BJP", "Woman without BJP  ",
                     "Man without BJP"), col = 1:4, lwd = 2, bty = "n")
title("Model-based prediction of survival function (BUN = 25, HB = 10, Age = 65)")

graphics.off()


## Residuals and goodness-of-fit of the Cox model
## ==============================================
residuals(cox2)

## Checking the linear assumption of the continuos variables"age", "BUN",
## and HB.
## ----------------------------------------------------------------------
resids1 <- residuals(update(cox2,  ~ . - age))
resids2 <- residuals(update(cox2,  ~ . - BUN))
resids3 <- residuals(update(cox2,  ~ . - HB))

windows(width = 15, height = 5)
par(mfrow = c(1, 3), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
plot(resids1 ~ myeloma$age, xlab = "Age", ylab = "Residuals", pch = 19)
abline(h = 0, lwd = 2, lty = 2)
lines(lowess(myeloma$age, resids1), lwd = 3)
plot(resids2 ~ myeloma$BUN, xlab = "BUN", ylab = "Residuals", pch = 19)
abline(h = 0, lwd = 2, lty = 2)
lines(lowess(myeloma$BUN, resids2), lwd = 3)
plot(resids3 ~ myeloma$HB, xlab = "HB", ylab = "Residuals", pch = 19)
abline(h = 0, lwd = 2, lty = 2)
lines(lowess(myeloma$HB, resids3), lwd = 3)


## The dfbeta residuals: looking for influential observations.
## -----------------------------------------------------------
dfbet <- residuals(cox2, type = "dfbeta")
dim(dfbet)
windows(width = 12, height = 7)
par(mfrow = c(2, 3), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
for (i in 1:5) {
  plot(dfbet[, i], pch = 16, ylab = "")
  title(names(coef(cox2))[i])
  axis(1, at = seq(5, 45, 5))
}

## Closer look at the influential observations on the estimation of a
## couple of parameters.
## Variable BUN:
#windows(width = 8)
par(font = 2, font.lab = 4, font.axis = 2, las = 1, cex.lab = 1.3,
    cex.axis = 1.2)
plot(dfbet[, 3], pch = 16, ylab = "")
title(names(coef(cox2))[3])
axis(1, at = seq(5, 45, 5))
identify(dfbet[, 3],
         labels = paste0("Row: ", rownames(myeloma), "; Time: ", myeloma$time,
                         "; Cens: ", myeloma$status, "; BUN: ", myeloma$BUN))

# ID 38 is an influential observation:
myeloma[38, ]      # Medium BUN value, but VERY short survival time.
# Estimation of beta_3 excluding ID38
coxph(Surv(time, status) ~ sex + BJP + BUN + HB + age, data = myeloma[-38, ])

## Variable HB:
windows(width = 8)
par(font = 2, font.lab = 4, font.axis = 2, las = 1, cex.lab = 1.3,
    cex.axis = 1.2)
plot(dfbet[, 4], pch = 16, ylab = "")
title(names(coef(cox2))[4])
axis(1, at = seq(5, 45, 5))
identify(dfbet[, 4],
         labels = paste0("Row: ", rownames(myeloma), "; Time: ", myeloma$time,
                         "; Cens: ", myeloma$status, "; HB: ", myeloma$HB))

# ID 13 is an influential observation
myeloma[13, ]      # Low HB value, but large survival time.
# Estimation of beta_4 excluding ID13
coxph(Surv(time, status) ~ sex + BJP + BUN + HB + age, data = myeloma[-13, ])


## The Schoenfeld residuals:
## Can we assume the assumption of proportional hazards hold?
## ----------------------------------------------------------
residuals(cox2, "schoenfeld")
## (i) Use of function cox.zph
cox.zph(cox2)

## (ii) Use of function plot.cox.zph
windows(width = 12, height = 7)
par(mfrow = c(2, 3), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
plot(cox.zph(cox2), lwd = 2)

graphics.off()

## =================================================================
## A very special R package of a former MESIO student (Jose Barrera)
## =================================================================
library(christmas)
{
 xmastree(2022)
 Sys.sleep(3)
 xmastreewire(2022, "spanish")
 Sys.sleep(3)
 xmassnowman(2022, "catalan")
}
