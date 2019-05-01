library(sharkbox)
library(RColorBrewer)
# par(mfrow=c(2,1))
#
# x = rbrokenbar(n=1000, S=4)
# y = rbrokenbar(n=1000, S=4, sequential=TRUE)
#
# barplot(t(x[order(x[,1], decreasing = TRUE),]), border=NA, space=0, col=1:4)
# barplot(t(y[order(y[,1], decreasing = TRUE),]), border=NA, space=0, col=1:4)

# 3 groups: dorado, billfish, sharks
# 7 species: dorado, 2 tunas, 2 billfish, 2 sharks


dat0 = read.csv("data/supersamples_demo.csv")
dat = split(dat0, f = dat0$trip)
dat = lapply(dat, FUN = function(x) x[, c(5,6)])

mod = smc(dat, S=13, L=3)
sim = predict(mod0)

S = nrow(mod$species$group_04)

#
# x = sim0(N=500)
# plot(x)
#
# mod = smcSim(G=G, S=S, L=L)
# sim = predict(mod)
#
# x = mod0
# names(x$groups$prop)
#
# x = mod0$species$group_04
# rownames(x)[rowSums(x)>0]

# np = 6
# par(mfrow=c(np,1), mar=c(3,1,1,1))
# for(i in 1:np) {
#   x = sim(N=N)
#   plot(x)
# }

# par(mfrow=c(2,1), mar=c(3,3,1,1))
# x = sim(N=N)
# plot(x)
# # table(x[,2])
# y = findLandingGroups(x[, 2], n=10)
# plot(y, type="s", col="red", pch=19, cex=0.3, xaxs="i", ylim=range(x[,1], y), lwd=2)
# lines(x[,1], type="s", col="blue", lwd=2)
# mtext(c("OBS", "SIM"), 3, line=c(-2,-3), adj=0.05, col=c("blue", "red"))

J = 5000
N = 500
X = list()
for(j in 1:J) X[[j]] = sim(N=N)
Y = lapply(X, FUN=function(x) x[, 2:3])
G = do.call(cbind, lapply(X, FUN=function(x) x[, 1]))

# x = X[[2]]
# par(mfrow=c(2,1), mar=c(3,3,1,1))
# plot(x)
# # table(x[,2])
# y = findLandingGroups(x[, 2])
# plot(y, type="s", col="red", pch=19, cex=0.3, xaxs="i", ylim=range(x[,1], y), lwd=2)
# lines(x[,1], type="s", col="blue", lwd=2)
# mtext(c("OBS", "SIM"), 3, line=c(-2,-3), adj=0.05, col=c("blue", "red"))

size = ceiling(10^seq(from=0.9, to=2.7, length=50))
rep = 50
err = array(dim=c(length(size), 2, rep))

for(i in seq_along(size)) {
  for(j in seq_len(rep)) {
    kali::DateStamp("Size =", size[i], ": rep =", j)
    ind = sample(length(X), size = size[i])
    fit0 = smc(data=Y[ind], S=S, L=L)
    fit1 = smc(data=X[ind], S=S, L=L)
    err[i, 1, j] = sum((mod$groups$prop - fit0$groups$prop)^2)
    err[i, 2, j] = sum((mod$groups$prop - fit1$groups$prop)^2)
  }
}

saveRDS(err, file="errorSMC.rds")

err0 = apply(err[, 1, ], 1, mean)
err0q = t(apply(err[, 1, ], 1, quantile, prob=c(0.05, 0.95)))
err1 = apply(err[, 2, ], 1, mean)
err1q = t(apply(err[, 2, ], 1, quantile, prob=c(0.05, 0.95)))

ylim = c(0, max(err0q, err1q))
plot.new()
plot.window(xlim=range(size), ylim=ylim)
matplot(size, cbind(err0, err1), type="l", lty=1, add=TRUE)
matplot(size, err0q, type="l", lty=3, add=TRUE, col="black")
title(xlab="Number of supersamples", ylab="error")
axis(1)
axis(2, las=1)
box()

barplot(rbind(mod$groups$prop, fit0$groups$prop, fit1$groups$prop), beside=TRUE)

# Old code ----------------------------------------------------------------

sapply(x$size, steadyStates)
sapply(x$species, steadyStates)
steadyStates(x$groups$TM(1000))




object = smc(data, ...) # create object of class xxx: estimate parameters for the model

predict(object, newdata) # create several replicates

DONE:
smcSim(G, S, L, ...) # create object of class xxx: simulate parameters for the model

sim = predict(object) # create a function: creates trajectories
sim(N=500) # create object of class xxy: a replicate of length 500


smc
smc_traj

library(kali)
kali::plot.map(domain="EPO")
