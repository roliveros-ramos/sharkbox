

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

x = newSMM(G=G, S=S, L=L) # species greater than groups
pp = as.numeric(1/(N*rbrokenbar(n=1, S=G)))
diag(x$groups) = 1 - pp
x$groups[row(x$groups) == (col(x$groups) - 1)] = head(pp, -1)
x$groups = x$groups/rowSums(x$groups)

main_species = sample(x=S, size=G)

.sortSpecies = function(first, size) {
  c(first, sample(setdiff(seq_len(size), first)))
}

.coverage = function(A) {
  out = sapply(seq_along(A), FUN=function(x) length(unique(A[seq_len(x)])))
  return(which.max(out))
}

A = t(sapply(main_species, .sortSpecies, size=S))
A[tail(seq_along(A), -.coverage(A))] = NA
A = apply(A, 1, FUN=function(x) c(na.omit(x)))

for(i in seq_len(G)) x$species[[i]][A[[i]], A[[i]]] =
  rbrokenbar(n=length(A[[i]]), S=length(A[[i]]), sequential=TRUE)

for(i in seq_len(S)) x$size[[i]] =
  rbrokenbar(n=nrow(x$size[[i]]), S=nrow(x$size[[i]]))



markovchain:::.ctmcEigen(x$size$species_1, transpose = FALSE)

function(object) {
  out<-.ctmcEigen(matr=transMatr, transpose=FALSE)
  if(is.null(out)) {
    warning("Warning! No steady state")
    return(NULL)
  }
   out <- - out / diag(object@generator)
  if(transposeYN==TRUE){
    out <- out / rowSums(out)
  }
  else{
    out <- out / colSums(out)
  }
  return(out)
}
