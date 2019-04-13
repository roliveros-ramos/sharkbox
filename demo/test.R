

par(mfrow=c(2,1))

x = rbrokenbar(n=1000, S=4)
y = rbrokenbar(n=1000, S=4, sequential=TRUE)

barplot(t(x[order(x[,1], decreasing = TRUE),]), border=NA, space=0, col=1:4)
barplot(t(y[order(y[,1], decreasing = TRUE),]), border=NA, space=0, col=1:4)

# 3 groups: dorado, billfish, sharks
# 7 species: dorado, 2 tunas, 2 billfish, 2 sharks

N = 500
G = 3
S = 7
L = 3

x = smcSim(G=G, S=S, L=L)

N = 100

y = .simMC(1000, x$groups$TM(1000))
plot(y)

sapply(x$size, steadyStates)
sapply(x$species, steadyStates)
steadyStates(x$groups$TM(1000))


x = x$groups
state = 1
sample(seq_len(nrow(x)), size=1, prob=x[state, ])


object = smc(data, ...) # create object of class xxx: estimate parameters for the model

sim = predict(object) # create a function: creates trajectories
sim(N=500) # create object of class xxy: a replicate of length 500
predict(object, newdata) # create several replicates

smcSim(G, S, L, ...) # create object of class xxx: simulate parameters for the model


smc
smc_traj



