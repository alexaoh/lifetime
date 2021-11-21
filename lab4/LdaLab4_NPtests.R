## =========================================================
## Lifetime Data Analysis, Course 2021/22
## 9.11.2021: Nonparametric tests to compare survival curves
## =========================================================
load("LdaLab4.RData")
library(survival)

## Two-sample hypothesis testing
## =============================
srem2t <- with(remission, Surv(time, cens) ~ treat)
svf2t <- survfit(srem2t)
svf2t

## The survival functions
## ----------------------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(svf2t, col = 3:4, xlab = "Time to relapse [Weeks]",
     lty = 1:2, ylab = expression(bold(hat(S)(t))), lwd = 3, yaxs = "i",
     bty = "l")
title("Survival functions among leukemia patients")
legend("bottomleft", levels(remission$treat), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)


## Comparison of both treatments with respect to the survival function
## -------------------------------------------------------------------


## The survdiff object
## -------------------
sdf <- survdiff(srem2t) # Default is rho = 0 --> logrank test. 
str(sdf)
sdf$obs; sdf$exp
sdf$chisq

## Further nonparametric tests
## ---------------------------
survdiff(srem2t) # logrank test. 
survdiff(srem2t, rho = 1) # Peto and Peto (and Prentice) modification of the Gehan-Wilcoxon test. Puts more weight in the beginning. 
survdiff(srem2t, rho = -1) # Puts more weight on later differences in the survival functions. 

## The Fleming-Harrington class
## ============================
library(KMsurv)
data(kidney)
head(kidney)
summary(kidney)
kidney$type <- factor(kidney$type, labels = c("Surgically", "Percutaneously"))

skid2t <- with(kidney, Surv(time, delta) ~ type)
survfit(skid2t)

## The figure on Slide 54/103
## --------------------------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(survfit(skid2t), col = 3:4, xlab = "Time to infection [Months]",
     lty = 1:2, ylab = expression(bold(hat(S)(t))), lwd = 3, yaxs = "i",
     bty = "l")
title("Survival functions according to catheter type")
legend("bottomleft", legend = levels(kidney$type), title = "Catheter placement",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)


## The table on Slide 55/103
## -------------------------
survdiff(skid2t)
survdiff(skid2t, rho = 1)

# install.packages("FHtest")
library(FHtest)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("Icens")

## The values of the table on Slide 55/103 (II)
## --------------------------------------------
FHtestrcc(skid2t)
FHtestrcc(skid2t, rho = 1)
FHtestrcc(skid2t, lambda = 1)
FHtestrcc(skid2t, rho = 1, lambda = 1)
FHtestrcc(skid2t, rho = 0.5, lambda = 2)


## What about the logrank test with left-truncation?
## =================================================
load("../lab3/LdaLab3.RData")

## Conditional survival function (T > 68)
## --------------------------------------
sub68 <- subset(channing, aged > 68)
(sch68 <- survfit(Surv(agee, aged, status) ~ sex, sub68))

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sch68, col = 2:3, lwd = 3, xlab = "Age [Years]",
     xlim = c(65, 102), ylab = expression(bolditalic(hat(S)(t))))
axis(1, at = seq(60, 100, 5))
title("Conditional survival functions of people older than 68 years")
legend("bottomleft", levels(channing$sex), col = 2:3, lwd = 3, bty = "n")

#survdiff(Surv(agee, aged, status) ~ sex, sub68) # Right censored data only!
summary(coxph(Surv(agee, aged, status) ~ sex, sub68)) # This works for truncated data (after fitting the cox model).
# This can hence be used to apply a non-parametric logrank test to left-truncated data.
## See: https://stat.ethz.ch/pipermail/r-help/2009-August/399999.html
## ==> Fit of a Cox model

## Stratified tests
## ================
? hodg
data(hodg)
head(hodg)
hodg$gtype <- factor(hodg$gtype, labels = c("Allogeneic", "Autologous"))
hodg$dtype <- factor(hodg$dtype, labels = c("Non Hodgkin lymphoma",
                                            "Hodgkins disease"))
summary(hodg)

shodg <- with(hodg, Surv(time, delta) ~ dtype + gtype)
survfit(shodg)

# Graphical representation (plot Slide 81/103)
# --------------------------------------------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(survfit(shodg), xlab = "Time to death or relapse [Days]",
     ylab = expression(bold(hat(S)(t))), col = rep(2:3, each = 2),
     lty = 1:2, lwd = 3, yaxs = "i", bty = "l")
title("Survival functions of time to death or relapse")
legend("bottomright", c("Group NHL allo", "Group NHL auto", "Group HOD allo",
                        "Group HOD auto"),
       bty = "n", col = rep(2:3, each = 2), lty = 1:2, lwd = 3)

# He says that it is not a good idea to use a stratified test for these data, but show how to fit it nonetheless. 
# Perhaps there are more comments in the R-file in Atenea. 
# He says there is an interaction effect, which means that they should be analyzed in each subgroup separetely.
# The interaction is apparent because: Groupd NHL allo is better then auto (in red), but it is the other way around for the green
# (in the green the auto is better than allo). Perhaps there is an interaction between the group and the allo/auto (if I understood correctly).
survdiff(Surv(time,delta) ~ gtype + strata(dtype), hodg) # strata gives stratified test. 
## With these data, a stratified test is not a good idea, because the effect of
## the graft type seems to be different among both diseases. We would better
## compare both graft types separately.
## Nontheless, let's see how to do a stratified test in R:

## Back to the "remission" data
## ----------------------------
svf <- survfit(Surv(time, cens) ~ status + treat, remission)
svf

windows(width = 10, height = 8)
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(svf, col = c(3, 3, 4, 4), xlab = "Time to relapse [Weeks]",
     lty = 1:2, ylab = expression(bold(hat(S)(t))), lwd = 3, yaxs = "i",
     bty = "l")
title("Survival functions of time to death or relapse")
legend("topright", names(svf$strata), bty = "n",
       col = c(3, 3, 4, 4), lty = 1:2, lwd = 3)

## Here, a stratified test is not necessary, because of the balanced study design:
with(remission, table(treat, status))
## Nonetheless ...
survdiff(Surv(time, cens) ~ treat + strata(status), remission)