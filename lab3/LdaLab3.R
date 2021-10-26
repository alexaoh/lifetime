## ==================================================================
## Lifetime Data Analysis, Course 2021/22
## 26.10.2021: Nonparametric estimation of S(t) under left truncation
## ==================================================================
setwd("/home/ajo/gitRepos/lifetime/lab3")
load("LdaLab3.RData")
ls()

## The Channing house data from Palo Alto
## ======================================
str(channing)
head(channing)
summary(channing)       # All persons were older than 60 when entering the study.


## A graphical representation of the observation intervals
## -------------------------------------------------------
par(las = 1, font = 2, font.lab = 2, font.axis = 2)
with(channing, plot(c(agee[1], aged[1]), c(1, 1), type = "l", xlim = c(60, 100),
                    ylim = c(1, nrow(channing)), xlab = "Age [Years]",
                    ylab = "Subject"))
for (i in 2:nrow(channing)) {
  if (channing$status[i] == 1) {
    with(channing, segments(agee[i], i, aged[i], i))
  } else {
    with(channing, segments(agee[i], i, aged[i], i, col = 2))
  }
}
axis(1, at = seq(60, 100, 5))
axis(2, at = seq(0, 450, 25))
title("Channing house data: observation intervals")
legend("bottomleft", c("Death", "Right censoring"), col = 1:2, lwd = 2,
       bty = "n")

# Left truncation in this data come from the fact that all people are not observed since the age of 62 (e.g), 
# i.e. that they enter in any other ages. 

# According to gender and ordered by "agee"
# -----------------------------------------
channord <- dplyr::arrange(channing, sex, agee)
par(las = 1, font = 2, font.lab = 2, font.axis = 2)
with(channord, plot(c(agee[1], aged[1]), c(1, 1), type = "l", xlim = c(60, 100),
                    ylim = c(1, nrow(channord)), xlab = "Age [Years]",
                    ylab = "Subject"))
for (i in 2:nrow(channord)) {
  if (channord$sex[i] == "Male") {
    with(channord, segments(agee[i], i, aged[i], i))
  } else {
    with(channord, segments(agee[i], i, aged[i], i, col = 2))
  }
}
axis(1, at = seq(60, 100, 5))
axis(2, at = seq(0, 450, 25))
title("Channing house data: observation intervals")
legend("topleft", c("Men", "Women"), col = 1:2, lwd = 2,
       bty = "n")


## The Surv object with left-truncated data
## ========================================
library(survival)
with(channing, Surv(agee, aged, status))
scha <- with(channing, Surv(agee, aged, status))
svcha <- survfit(scha ~ sex, channing)
svcha
summary(svcha, times = seq(60, 100, 5))

# The explanation is the gap in the amount of men at risk around year 65, as discussed in class. 
# i.e. zero men at risk at this year, which gives a 0 product estimator (estimated survival function), 
# which later cannot be recovered, even though there will be people at risk later than this point 
# (i.e. the factors after this point will not be zero, but the whole product will stay at zero forever.)


## The Kaplan-Meier estimation of the survival function
## ====================================================
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(svcha, col = 3:2, lwd = 3, xlab = "Age [Years]", xlim = c(60, 105),
     ylab = expression(bolditalic(hat(S)(t))))
legend("topright", c("Women", "Men"), col = 2:3, lwd = 3, bty = "n")


## Number of individuals at risk over time (Slide 90/101)
## ======================================================
Times <- seq(min(channing$agee), max(channing$aged), 0.01)
NriskF <- NriskM <- numeric(length(Times))
for (i in 1:length(Times)) {
  NriskF[i] <- with(channing, sum(agee <= Times[i] & aged >= Times[i] &
                                  sex == "Female"))
  NriskM[i] <- with(channing, sum(agee <= Times[i] & aged >= Times[i] &
                                  sex == "Male"))
}


par(font = 2, font.axis = 2, font.lab = 4, las = 1)
plot(Times, NriskF, type = "S", lwd = 3, xlab = "Age [Years]",
                    ylab = "Individuals at risk", bty = "n", col = 2,
                    xlim = c(60, 105), ylim = c(0, 175), yaxs = "i", xaxs = "i")
lines(Times, NriskM, type = "S", lwd = 3, col = 3)
axis(1, at = seq(60, 100, 5))
axis(2, at = seq(0, 175, 25))
title("Channing House data: number of persons at risk")
legend("left", c("Women", "Men"), col = 2:3, lwd = 3, bty = "n")


## Conditional survival functions (Slide 93/101)
## ---------------------------------------------
sub68 <- subset(channing, aged > 68)
sub80 <- subset(channing, aged > 80)

sch68 <- survfit(Surv(agee, aged, status) ~ sex, sub68)
sch80 <- survfit(Surv(agee, aged, status) ~ sex, sub80)

# 2nd step: Graph
par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sch68, lwd = 2, col = 3:2, lty = 2, xlab = "Age [Years]",
     xlim = c(65, 102), ylab = expression(bolditalic(hat(S)(t))))
lines(sch80, lwd = 2, col = 3:2, lty = 1)
axis(1, at = seq(60, 100, 5))
title("Conditional survival functions of people older than 68 and 80 years")
legend("bottomleft", paste0(c("Women", "Men"), rep(c(" > 68", " > 80"), each = 2)),
       col = 2:3, lty = c(2, 2, 1, 1), lwd = 3, bty = "n")


## Exercise: Draw the conditional survival functions S(T | T > 68)
## considering and ignoring left truncation. What do you observe
## and which is the explanation for it?
## ---------------------------------------------------------------
par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sch68, lwd = 2, col = 3:2, lty = 2, xlab = "Years",
     xlim = c(65, 102), ylab = expression(bolditalic(hat(S)(t))))
lines(survfit(Surv(aged, status) ~ sex, sub68), lwd = 2, col = 3:2, lty = 1) # agee is not used (do not take truncation into account)
axis(1, at = seq(60, 100, 5))
title("Conditional survival functions of people older than 68 years")
legend("bottomleft", col = 2:3, lty = c(2, 2, 1, 1), lwd = 2, bty = "n",
       legend = paste0(c("Women", "Men"),
                       rep(c(": Considering ", ": Ignoring "), each = 2),
                       "truncation"))

# The survival function is overestimated when ignoring left truncation (because of the introduced bias), 
# as noted in the class. 


graphics.off()
