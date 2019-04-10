
n = 100000
sizes = 1:10
maxs = matrix(nrow=n, ncol=length(sizes))
for(N in sizes)
  for(i in 1:n) maxs[i, N] = max(diff(sort(c(0, runif(N), 1))))

mins = matrix(nrow=n, ncol=length(sizes))
for(N in sizes)
  for(i in 1:n) mins[i, N] = min(diff(sort(c(0, runif(N), 1))))


par(mfrow=c(5,2), mar=c(4,4,1,1))
for(i in sizes) hist(maxs[, i], breaks=100)

par(mfrow=c(5,2), mar=c(4,4,1,1))
for(i in sizes) hist(mins[, i], breaks=100)

N = 4
plot(sort(diff(sort(c(0, runif(N), 1))), decreasing = TRUE), type="h")



par(mfrow=c(2,1))

x = rbrokenbar(n=1000, S=4)
y = rbrokenbar(n=1000, S=4, sequential=TRUE)

barplot(t(x[order(x[,1], decreasing = TRUE),]), border=NA, space=0, col=1:4)
barplot(t(y[order(y[,1], decreasing = TRUE),]), border=NA, space=0, col=1:4)








