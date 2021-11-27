# Exercise 4 (turnover)
library(FHtest) # Fleming Harrington tests.

setwd("/home/ajo/gitRepos/lifetime/Exercises2")
load("Assign2Exer4.RData")
head(turnover)
summary(turnover)

# a) Draw the survival functions of time until turnover for both men and women. Observations?
surv1 <- with(turnover, Surv(stag, event)) # Make survival object. 
summary(surv1)
sfit1 <- survfit(surv1 ~ gender, data = turnover) # Fit the survival curves (estimation with KM).
sfit1
summary(sfit1)

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sfit1, col = 3:4, xlab = "Time in company until turnover [Months]",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to gender")
legend("bottomleft", legend = levels(as.factor(turnover$gender)), title = "Gender",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)
axis(1, at = seq(0, 200, 50))
axis(2, at = seq(0, 1, 0.1))


# b) Survival curves separately for both genders of the supervisor. 
data.male.supervisor <- turnover[turnover$headgend == "Male", ]
surv.male.super <- with(data.male.supervisor, Surv(stag, event))
sfit.male.super <- survfit(surv.male.super ~ gender, data = data.male.supervisor)

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sfit.male.super, col = 3:4, xlab = "Time in company until turnover [Months]",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to gender with male supervisor")
legend("bottomleft", legend = levels(as.factor(turnover$gender)), title = "Gender",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)
axis(1, at = seq(0, 200, 50))
axis(2, at = seq(0, 1, 0.1))

data.female.supervisor <- turnover[turnover$headgend == "Female", ]
surv.female.super <- with(data.female.supervisor, Surv(stag, event))
sfit.female.super <- survfit(surv.female.super ~ gender, data = data.female.supervisor)

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sfit.female.super, col = 3:4, xlab = "Time in company until turnover [Months]",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1:2, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Survival functions according to gender with female supervisor")
legend("bottomleft", legend = levels(as.factor(turnover$gender)), title = "Gender",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)
axis(1, at = seq(0, 200, 50))
axis(2, at = seq(0, 1, 0.1))

# c) Test of hypothesis with Fleming Harrington. Testing the hypothesis:
# Time to turnover does not depend on the employee's gender.
skd <- with(turnover, Surv(stag, event) ~ gender)
FHtestrcc(skd) # rho = 0 and lambda = 0 is logrank.
FHtestrcc(skd, rho = 1) # Emphasizes earlier differences (lambda = 0) and positive rho.
FHtestrcc(skd, lambda = 1) # Emphasizes later differences (rho = 0) and positive lambda.
FHtestrcc(skd, rho = 1, lambda = 1) # Emphasized differences towards the middle (rho = lambda)
FHtestrcc(skd, rho = 0.5, lambda = 2) # Emphasizes differently also, more towards the end since lambda > rho I think. 
# Anyhow, none of these tests give significant results to any reasonable level. 

# d) Test the hypothesis with a stratified test. What is the conclusion now. 
survdiff(Surv(stag, event) ~ gender + strata(headgend), data = turnover) # strata gives stratified test. 
