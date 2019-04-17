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

N = 300
G = 3
S = 7
L = 3

mod = smcSim(G=G, S=S, L=L)
sim = predict(mod)

# np = 6
# par(mfrow=c(np,1), mar=c(3,1,1,1))
# for(i in 1:np) {
#   x = sim(N=N)
#   plot(x)
# }

par(mfrow=c(2,1), mar=c(3,3,1,1))
x = sim(N=N)
plot(x)
# table(x[,2])
y = findLandingGroups(x[, 2], n=10)
plot(y, type="s", col="red", pch=19, cex=0.3, xaxs="i", ylim=range(x[,1], y), lwd=2)
lines(x[,1], type="s", col="blue", lwd=2)
mtext(c("OBS", "SIM"), 3, line=c(-2,-3), adj=0.05, col=c("blue", "red"))

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



