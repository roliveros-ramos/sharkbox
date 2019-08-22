.getTM = function(x, S=NULL, simplify=FALSE, diagonal2zero=FALSE) {
  if(is.null(S)) S = max(x)
  x = factor(x, levels=seq_len(S))
  out = table(head(x, -1), tail(x, -1))
  if(isTRUE(simplify)) {
    cs = colSums(out)
    rs = rowSums(out)
    ind = which((cs==0) & (rs==0))
    if(length(ind)>0) {
      out = out[-ind, ]
      out = out[, -ind]
    }
  }
  if(isTRUE(diagonal2zero)) diag(out) = 0
  return(out)
}

.countTraj_old = function(x, S, L=3, thr=0.03) {
  # x has two columns: species and size-class
  if(ncol(x)>2) {
    y = x[, 1]
    x = x[, 2:3]
  } else {
    y = findLandingGroups(x[, 1], S=S, thr=thr)
  }
  x0 = tapply(x[,1], INDEX = y, FUN = .getTM, S=S)
  x1 = tapply(x[,2], INDEX = x[,1], FUN = .getTM, S=L)
  mainSp = apply(table(y, factor(x[,1], levels=seq_len(S))), 1, which.max)
  nn = newSMC(G=S, S=S, L=L)
  nn$species[mainSp] = x0
  nn$groups$prop[mainSp] = table(y)
  nn$groups$jump[mainSp, mainSp] = .getTM(x=y, simplify=TRUE)
  diag(nn$groups$jump) = 0
  nn$size[as.numeric(names(x1))] = x1
  return(nn)
}
