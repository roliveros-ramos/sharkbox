

gp = interpret.smc(order|trip ~ b(block0) + mc(species, by=block0) + m(group, by=species) + mc(size, by=group))

term = lapply(gp$smooth.spec, FUN="[[", i="term")
by = lapply(gp$smooth.spec, FUN="[[", i="by")

exec = rapply(term, function(x) NA, how = "list")


x = as.factor(dat$species)
levels = levels(x)
x = as.numeric(x)
by = as.factor(dat$group)
byLevels = levels(by)
by = as.numeric(by)
yy = .getSimpleTM(x, levels, by, byLevels)




ff = interpret.smc(order|trip ~ mc(species, size))




mod = smc(dat, S=23, L=3)
sim = predict(mod)

x = sim(N=500)


spp = dat0[, c("sp_code", "species", "fao_code")]
spp = spp[!duplicated(spp), ]
spp = spp[order(spp$sp_code), ]



ff = interpret.smc(order|trip ~ species)

xx = glm(order ~ species, data=dat0)

.addLG = function(x) {
  x = cbind(landing_group=as.numeric(findLandingGroups(x[,1])), x)
  x = do.call(cbind, x)
  class(x) = c("smc_traj", class(x))
  return(x)
}

dat2 = lapply(dat, .addLG)

# spp$species = as.character(spp$species)
# spp$species[2] = "Mustelus sp."

N = 500
par(mfrow=c(20,1), mar=c(1,0,1,0), oma=c(6,2,1,2))
lapply(dat2[1:20], plot)
axis(1)
mtext("Number of unloaded fish", 1, outer=TRUE, line=2)
mtext("Número de peces descargados", 1, outer=TRUE, line=3.5)

colors = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))[-c(6,11)]
par(new=TRUE)
plot.window(xlim=c(0,1), ylim=c(0.55,0.95))
text(x = 0.3, y=seq(0.9,0.6, length=4), labels = spp$species[1:4], col=colors[1:4], cex=1.2, font=2)
text(x = 0.6, y=seq(0.9,0.6, length=4), labels = spp$species[5:8], col=colors[5:8], cex=1.2, font=2)
text(x = 0.9, y=seq(0.9,0.6, length=4), labels = spp$species[9:12], col=colors[9:12], cex=1.2, font=2)


plot(x)
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

size = c(1:9, seq(10, 100, by=2), seq(105, 200, 5))
rep = 200
err = array(dim=c(length(size), 2, rep))

mod_prop0 = setNames(numeric(S), nm=sprintf("group_%02d", seq_len(S)))
mod_prop = mod_prop0
mod_prop[names(mod$groups$prop)] = mod$groups$prop

for(i in seq_along(size)) {
  for(j in seq_len(rep)) {
    kali::DateStamp("Size =", size[i], ": rep =", j)
    ind = sample(length(X), size = size[i])

    fit0 = try(smc(data=Y[ind], S=S, L=L))
    if(inherits(fit0, "try-error")) next

    fit1 = try(smc(data=X[ind], S=S, L=L))
    if(inherits(fit1, "try-error")) next

    fit0_prop = fit1_prop = mod_prop0
    fit0_prop[names(fit0$groups$prop)] = fit0$groups$prop
    fit1_prop[names(fit1$groups$prop)] = fit1$groups$prop
    err[i, 1, j] = sum((mod_prop - fit0_prop)^2)
    err[i, 2, j] = sum((mod_prop - fit1_prop)^2)
  }
}

saveRDS(err, file="errorSMC2.rds")

err = readRDS(file="errorSMC2.rds")
err0 = apply(err[, 1, ], 1, median, na.rm=TRUE)
err0q = t(apply(err[, 1, ], 1, quantile, prob=c(0.05, 0.95), na.rm=TRUE))
err1 = apply(err[, 2, ], 1, median, na.rm=TRUE)
err1q = t(apply(err[, 2, ], 1, quantile, prob=c(0.05, 0.95), na.rm=TRUE))

ylim = c(0, max(err0q, err1q))
plot.new()
plot.window(xlim=c(1,50)*500, ylim=ylim, xaxs="i")
# matplot(size, cbind(err0, err1), type="l", lty=1, add=TRUE)
polygon(x=500*c(size, rev(size)), y=c(err0q[,1], rev(err0q[,2])),
        col="grey90", border=NA)
lines(size*500, err0, lwd=2)
# matplot(size*500, err0q, type="l", lty=3, add=TRUE, col="black")
title(xlab="Number of fish in super-samples - Número de peces en las supermuestras",
      ylab="Error")
mtext("Number of super-samples - Número de supermuestras", 3, line=2.5)
axis(1)
axis(2, las=1)
sax = sort(c(6*200, axTicks(3)))
axis(3, at=sax, labels = sax/200)
# abline(v=c(70, 80)*200, lty=2)
abline(v=1200, lty=2, col="red") # aprox
box()

barplot(rbind(mod$groups$prop, fit0$groups$prop, fit1$groups$prop), beside=TRUE)

# Old code ----------------------------------------------------------------

plot(dat[[1]])



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


x = mod$species$group_04
.getSpecies = function(x, species=NULL, merge=TRUE) {
  out = sapply(strsplit(colnames(x)[which(colSums(x) > 0)], split="_"), "[", i=2)
  out = as.numeric(out)
  if(!is.null(species)) {
    out = species[out]
    out = as.character(out)
    if(isTRUE(merge)) out = paste(out, collapse=", ")
  }
  return(out)
}

lapply(mod$species, .getSpecies, species=spp$species)

