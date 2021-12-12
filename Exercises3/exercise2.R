# Generate survival times from a log-normal distribution and check whether the Weibull or the log-logistic distribution fits better to the data.

# a) 
RVT <- rlnorm(300, meanlog = log(2), sdlog = 1) # Gir ikke sd(RVT) = 1?

# b) 
RVC <- rexp(300, rate = 1/20)

# c) 
## Observed survival times and event indicator
obs <- pmin(RVT, RVC)
cens <- as.numeric(RVT <= RVC)
table(cens)

## Checking
head(cbind(RVT, RVC, obs, cens))

# d) Draw the cumulativ hazard plots for Weibull and log-logistic. Which model fits better?
library(GofCens)
#cumhazPlot(obs, cens, col = 4, distr = c("wei", "loglo"), ggplo = T)
cumhazPlot(obs, cens, col = 4, distr = c("wei", "loglo"))
# Fungerer ikke med libraryen nå!
# NOE GALT MED RVT!

# Manuelle plot: 
svf <- survfit(Surv(obs, cens) ~ 1, stype = 2, ctype = 1)
svf
summary(svf)

## Uncensored survival times
times <- summary(svf)$time

## Nelson-Aalen estimate of the cumulative hazard function
chaz <- -log(summary(svf)$surv)

## The probability plots
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
