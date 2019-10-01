library(colorful)

# Sampling design experiments

sp = "Carcharhinus falciformis"

# ad-hoc function to calculate size composition
getSizeProportion = function(x, sp, prob=TRUE) {
  .getSizeProportion = function(x, sp, prob) {
    y = factor(x$size, levels = 1:3)
    counts = table(y[x$species==sp])
    if(!isTRUE(prob)) return(as.numeric(counts))
    out = sharkbox:::.norm(counts)
    out = as.numeric(out)
    return(out)
  }
  out = t(sapply(x, FUN=.getSizeProportion, sp=sp, prob=prob))
  return(out)
}

# observed size composition
obs = getSizeProportion(unloading_data, sp=sp)

.multinomial = function(observed, prob) {
  out = numeric(nrow(observed))
  for(i in seq_len(nrow(observed))) {
    if(sum(observed[i,]<1)) out[i] = 1 else {
      out[i] = ExactMultinomialTest(observed = observed[i,], prob = prob[i,])$p.value
    }
  }
  return(out)
}

skip     = seq(0, 10, by=1)
start    = seq(0, 0.9, by=0.01)
fraction = seq(0.1, 1, by=0.01)
# sample size, cost

X = outer(start, fraction, FUN = "+")

error1 = array(NA, dim=c(length(start), length(fraction), length(skip)))
error2 = array(NA, dim=c(length(start), length(fraction), length(skip)))

for(k in seq_along(skip)) {
  kali::DateStamp(k)
  for(i in seq_along(start)) {
    for(j in seq_along(fraction)) {
      if(X[i,j] > 1) next
      xsim = getSizeProportion(sampling(unloading_data,
                                       skip=skip[k], start=start[i], fraction=fraction[j]), sp=sp, prob=FALSE)
      error1[i,j,k] = mean(.multinomial(xsim, obs))
      error2[i,j,k] = sum((obs - xsim)^2)
    }
  }
}

# diagnostics, size bin categories, NA for fraction


k=1
zlim = c(0,1)
col = biasPalette(zlim=zlim)

par(mfrow=c(2,2), mar=c(4,4,3,1))
for(k in c(0,1,3,5)+1)
image.plot(start, fraction, error1[,,k],
           main=sprintf("Skip = %s", k-1),
           zlim=zlim, col=col)


