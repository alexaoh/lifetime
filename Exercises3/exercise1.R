rm(list = ls())
setwd("/home/ajo/gitRepos/lifetime/Exercises3")
library(KMsurv)
library(rms)
library(survival)
library(FHtest)

data(tongue)
?tongue
head(tongue)
str(tongue)

# a) Draw the survival curves. 
s <- with(tongue, Surv(time, delta))
t <- npsurv(s ~ type, tongue)
t

par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
survplot(t, ylab = expression(bold(hat(S)(t))), col = 1:2, lwd = 3, conf = "none")
title("Survival functions according to tumour DNA profile")

# b) Logrank test. 
s2 <- with(tongue, Surv(time, delta) ~ type)
FHtestrcc(s2)

# c) Fit log-logistic model and use the model to test the same hypothesis. 
loglogistic <- survreg(s ~ type, data = tongue, dist = "loglogistic")
summary(loglogistic)
loglo.pred <- predict(loglogistic, type = "linear")
resids.loglo <- (log(tongue$time) - loglo.pred) / loglogistic$scale

par(font = 2, font.axis = 2, font.lab = 2, las = 1, mar = c(5, 5, 4, 2))
plot(survfit(Surv(resids.loglo, tongue$delta) ~ 1), col = c(1,2,2), xlab = "Residuals",
     ylab = expression(bolditalic(hat(S)(t))),
     lty = 1, lwd = 3, yaxs = "i", xaxs = "i", bty = "n")
title("Residuals of the log-logistic Regression Model")
curve(plogis(x, lower.tail = F), from = min(resids.loglo), to = max(resids.loglo), col = 3, lwd = 3,
      add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Logistic Distribution"),
       col = c(1, 2, 3), lty = c(1, 2, 1), lwd = 3, bty = "n")
