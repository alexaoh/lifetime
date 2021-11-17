setwd("/home/ajo/gitRepos/lifetime/Exercises2")
data <- read.csv("Elderly.txt") # hard to manage.
# Made an Elderly2.txt, where the rows with text have been removed. 
d <- read.table("Elderly2.txt")
colnames(d) <- c("AgeEntry", "AgeExit", "Ind")
str(d)
head(d)
summary(d)

# a) Explain why the data are left-truncated. Explained in pdf-document. 

# b) Draw the number of people at risk of dying as a function of age.
Times <- seq(min(d$AgeEntry), max(d$AgeExit), 0.01)
Nrisk <- numeric(length(Times))
for (i in 1:length(Times)) {
  Nrisk[i] <- with(d, sum(AgeEntry <= Times[i] & AgeExit >= Times[i]))

}


par(font = 2, font.axis = 2, font.lab = 4, las = 1)
plot(65 + Times, Nrisk, type = "S", lwd = 3, xlab = "Age [Years]",
     ylab = "Individuals at risk", bty = "n", col = 2,
    yaxs = "i", xaxs = "i")
axis(1, at = seq(65, 45+65, 5))
axis(2, at = seq(0, 45, 5))
title("Elderly data: number of men at risk")
#legend("left", c("Men"), col = 2, lwd = 3, bty = "n")


# c) Conditional survival functions.
# Translation in x to make the next tasks a bit easier.
d2 <- d
d2$AgeEntry <- d2$AgeEntry + 65
d2$AgeExit <- d2$AgeExit + 65

sub70 <- subset(d2, AgeExit > 70)
sub85 <- subset(d2, AgeExit > 85)

sch70 <- survfit(Surv(AgeEntry, AgeExit, Ind) ~ 1, sub70)
sch85 <- survfit(Surv(AgeEntry, AgeExit, Ind) ~ 1, sub85)
summary(sch70)
summary(sch85)

# 2nd step: Graph
par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sch70, conf.int = F, lwd = 2, col = "red", lty = 1, xlab = "Age [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xlim = c(65, 108), ylim = c(0,1))
lines(sch85, conf.int = F, lwd = 2, col = "blue", lty = 2)
axis(1, at = seq(65, 105, 5))
title("Conditional survival functions of men older than 70 and 85 years")
legend("bottomleft", legend = c(" > 70", " > 85"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2, bty = "n")

# d) Corresponding probabilites of surviving 90 years when ignoring truncation.

# AgeEntry is not used (do not take truncation into account)
sch70.nontrunc <- survfit(Surv(AgeExit, Ind) ~ 1, sub70)
sch85.nontrunc <- survfit(Surv(AgeExit, Ind) ~ 1, sub85)

summary(sch70.nontrunc)
summary(sch85.nontrunc)

par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sch70, conf.int = F, lwd = 2, col = "red", lty = 1, xlab = "Age [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xlim = c(65, 108), ylim = c(0,1))
lines(sch85, conf.int = F, lwd = 2, col = "blue", lty = 1)
lines(sch70.nontrunc, conf.int = F, lwd = 2, col = "red", lty = 2)
lines(sch85.nontrunc, conf.int = F, lwd = 2, col = "blue", lty = 2)
axis(1, at = seq(65, 105, 5))
title("Conditional survival functions of men older than 70 and 85 years")
legend("bottomleft", legend = c(" > 70", " > 85", " > 70 non-trunc", " > 85 non-trunc"),
       col = c("red", "blue", "red", "blue"), lty = c(1, 1, 2, 2), lwd = 2, bty = "n")

