## =====================================================
## Lifetime Data Analysis, Course 2021/22
## 14.12.2021: Fit of the Cox Proportional Hazards model
## =====================================================
setwd("/home/ajo/gitRepos/lifetime/lab7")
load("LdaLab7.RData")
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
#windows(width = 14)
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
# Here: exp(\beta_2) = 1.6 is the HR comparing NoBJP vs. BJP for only women. (since women is the reference category)
# exp(\beta_2 + \beta_3) = 1.12 is the HR comparing NoBJP vs. BJP for only men. 
# See some notes in my book from the lecture! (Lab7_CoxModel)
# exp(\beta_3) = 1.12 cannot be interpreted as a HR. 

## The model-based hazard ratios associated with the absence of the BJ protein
## ---------------------------------------------------------------------------
library(Epi)
ci.lin(cox1)
round(ci.lin(cox1, Exp = TRUE), 3)
(ctmat <- matrix(c(0,1,0,0,1,1), byrow = TRUE, nr = 2))
round(ci.lin(cox1, ctr.mat = ctmat, Exp = TRUE), 3)

# A somewhat nicer presentation
HRmat <- round(ci.lin(cox1, ctr.mat = ctmat, Exp = TRUE), 3)[, c(1, 5:7)]
rownames(HRmat) <- c("HR| Women", "HR| Men")
colnames(HRmat) <- c("logHR", "HR", "Lower 95%", "Upper 95%")
HRmat


## The Cox model including, also, BUN, HB, and age
## -----------------------------------------------
cox2 <- update(cox1,  ~ . + BUN + HB + age - sex:BJP, myeloma)
summary(cox2)

## Which is the hazard ratio associated with an increase of 10 units in BUN?
## Which is the corresponding 95% confidence interval?
## -------------------------------------------------------------------------


## Estimation of the baseline hazard function
## ------------------------------------------
coxph.detail(cox2)
coxph.detail(cox2)$hazard


## Model-based prediction
## ======================
cox2$mean
survfit(cox2)
#windows(width = 8)
par(font = 2, font.lab = 4, font.axis = 2, las = 1)
plot(survfit(cox2), lwd = 2, col = 3)

## Patient profiles for prediction
## -------------------------------
(newdat <- data.frame(BUN = 25, sex = factor(rep(c("Female", "Male"), 2)),
                      HB = 10, age = 65,
                      BJP = factor(rep(c("No", "Yes"), each = 2))))

survfit(cox2, newdata = newdat)
summary(survfit(cox2, newdata = newdat))

## Graphical respresentation
## -------------------------
#windows(width = 8)
par(font = 2, font.lab = 4, font.axis = 2, las = 1)
plot(survfit(cox2, newdata = newdat), col = 1:4, lwd = 3, xlab = "Time [Months]",
     bty = "n", xlim = c(0, 100))
legend("bottomleft", c("Woman without BJP", "Man without BJP",
                       "Woman with BJP", "Man with BJP"),
       col = 1:4, lwd = 2, bty = "n")
title("Model-based prediction of survival function (BUN = 25, HB = 10, Age = 65)")

graphics.off()
