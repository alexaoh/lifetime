## ===============================================
## Lifetime Data Analysis, Course 2021/22
## 9.11.2021: Function survplot of the rms package
## ===============================================
setwd("/home/ajo/gitRepos/lifetime/lab4")
library(survival)
load("LdaLab4.RData")
ls.str()

srem <- with(remission, Surv(time, cens))
svf.treat <- survfit(srem ~ treat, remission)

# The plot we have seen already several times
# -------------------------------------------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
plot(svf.treat, col = 3:4, xlab = "Time to relapse [Weeks]", mark.time = TRUE,
     lty = 1:2, ylab = expression(bold(hat(S)(t))), lwd = 3, yaxs = "i",
     bty = "l", cex = 1.5)
title("Survival functions according to treatment")
legend("bottomleft", levels(remission$treat), title = "Treatment",
       bty = "n", col = 3:4, lty = 1:2, lwd = 3)

# The same curves drawn with function survplot of package rms
# -----------------------------------------------------------
library(rms)

svf.treat2 <- npsurv(srem ~ treat, remission)
svf.treat2
summary(svf.treat2)
str(svf.treat2)

# Use of function survplot() instead of plot.survfit()
# ====================================================
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
survplot(svf.treat2, ylab = expression(bold(hat(S)(t))), col = 1:2, lwd = 3)
title("Survival functions according to treatment")

# Example 2
# ---------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
survplot(svf.treat2, ylab = expression(bold(hat(S)(t))), col = 1:2, lwd = 3,
         col.fill = grey(c(0.75, 0.85)))
title("Survival functions according to treatment")

# Example 3
# ---------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
survplot(svf.treat2, ylab = expression(bold(hat(S)(t))), col = 1:2, lwd = 3,
         conf = "none")
title("Survival functions according to treatment")

# Example 4
# ---------
par(font = 2, font.axis = 2, font.lab = 4, las = 1, mar = c(5, 5, 4, 2))
survplot(svf.treat2, ylab = expression(bold(hat(S)(t))), col = 1:2, lwd = 3,
         col.fill = grey(c(0.75, 0.85)), n.risk = TRUE)
title("Survival functions according to treatment")

graphics.off()
