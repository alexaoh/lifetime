## =====================================================
## Lifetime Data Analysis, Course 2021/22
## 23.11.2021: Fit of the Accelerated Failure Time Model
## =====================================================
setwd("/home/ajo/gitRepos/lifetime/lab5")
load("LdaLab5.RData")
library(survival)

# A data set on larynx cancer
# ===========================
str(larynx)
head(larynx, 10)
summary(larynx)

# Does the assumption of non-informative censoring hold in this case?
# (Censoring times independent of survival times).
# This is what we assume in (almost) ALL our model/theory in this course!
# Yes, we can make this assumption correctly here, since the end time of the study
# is the same for all. This means that the right-censoring occurs at the same times for all people.
# No left-truncation either, since all people are followed from starting the treatment. 
# The censoring is not dependent on the stage of the disease.

# Could be a good idea to google a bit about this assumption to make sure that I really understand it!!

table(larynx$stage)

# Kaplan-Meier estimates of the survival functions
# ------------------------------------------------
larsurv <- with(larynx, Surv(time, cens))

# Two KM plots
par(mfrow = c(1, 2), font = 2, font.lab = 4, font.axis = 2, las = 1,
    oma = c(0, 0, 1, 0), mar = c(5, 5, 4, 2))
plot(survfit(larsurv ~ stage, larynx), col = 1:4, xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))))
title("Survival functions according to disease stages")
legend("bottomleft", paste("Stage", 1:4), col = 1:4, lwd = 3, bty = "n")
plot(survfit(larsurv ~ yearc, larynx), col = 1:3, xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))))
title("Survival functions according to year of diagnosis")
legend("bottomleft", levels(larynx$yearc), col = 1:3, lwd = 3, bty = "n")
title("Survival in patients with larynx cancer", outer = TRUE)


# Function survreg() fits a parametric survival model
# ===================================================
args(survreg)

# The null model
# --------------
survreg(larsurv ~ 1) # The null model is: ln(T) = \mu + W\sigma (no covariates Z in addition to the intercept).
summary(survreg(larsurv ~ 1))
survreg(larsurv ~ 1, dist = "loglo")           # Log-logistic model
summary(survreg(larsurv ~ 1, dist = "loglo"))



# We continue with the Weibull model ...
weimod1 <- survreg(larsurv ~ 1)

# ... and include variable stage in the model
summary(survreg(larsurv ~ stage, larynx))
# What do you think about that model?

# Should instead regress with the function as a factor. 
larynx$stagec <- factor(larynx$stage, labels = paste("Stage", 1:4))

# A Weibull regression model for variable stage
# =============================================
weimod2 <- survreg(larsurv ~ stagec, larynx)
summary(weimod2)

# One possible effect size measure: 
# Acceleration factor: AF = e^{-\gamma}. Here: comparing survival in stage i with survival in stage 1.
# i.e. the median of the survival in stage 4 is accelerated by e^(1.57) \approx 4.8 compared 
# to the median of the survival in stage 1. (not only the median, but also any other quantile in the data).
# Of course, based on this model! The same can be done for comparing each of the other stages with stage 1 respectively.
# This interpretation is the same in any of the accelerated failure time models - does not depend on choice of distributions. 

# Another effect size measure DOES DEPEND on the distribution chosen for the model: 
# If T \sim Weibull --> this is a proportional hazards model. 
# Then: e^{-\gamma/\sigma} is the hazards ratio. 
# e.g. e^(1.5786/0.885) \approx 5.95, which is the ratio between the hazard functions. (in stage 4 divided by stage 1)
# In this proportional hazards model this proportion is assumed to be constant - only valid when T \sim Weibull. 
# We are comparing two patients in two different stages, with the same survival time. The hazard in that instance is approx 
# 6 times larger for stage 4 vs stage 1. A bit harder to interpret then the earlier measure, because of the 
# more difficult interpretation of the hazard functions. NOT A difference in probability for example!

# The survreg.object contains a lot of information:
# -------------------------------------------------
str(weimod2)

# The parameter estimates
weimod2$coefficient
round(summary(weimod2)$table, 3)
# The scale parameter
weimod2$scale

# Function anova() can be used to test the global hypothesis
# ----------------------------------------------------------
anova(weimod2)


## Pairwise comparisons (beyond the scope of our course, but nice to see anyway).
## --------------------
library(multcomp)
glht(weimod2, linfct = mcp(stagec = "Tukey"))
print(summary(glht(weimod2, linfct = mcp(stagec = "Tukey"))),
      signif.stars = FALSE)


# A Weibull regression model for variables stage and age
# ======================================================
(weimod3 <- update(weimod2, ~. + age)) # Can use this to update the fit, instead of doing it completely again!
# (weimod4 <- survreg(larsurv ~ stagec + age, larynx)) # The same as this!
summary(weimod3)
anova(weimod3)


## Little exercise:
## · Compute the HRs associated with Stages 2 to 4 with respect to Stage 1.
## · Compute the AFs associated with Stages 2 to 4 with respect to Stage 1.
## · Which is the hazard ratio associated to an age difference of 5 years?
## -----------------------------------------------------------------------------



## Comparison with the log-logistic model
## --------------------------------------
loglomod3 <- update(weimod3, dist = "loglogistic")
summary(loglomod3)
# Which is the interpretation of exp(gamma / sigma) under this model?
with(loglomod3, round(exp(coefficients[2:4] / scale), 2))
