
.completeTails = function(x) {
  nn = length(x)
  nona = which(!is.na(x))
  h = head(nona, 1)
  t = tail(nona, 1)
  x[seq_len(h-1)] = x[h]
  x[t+seq_len(nn-t)] = x[t]
  return(x)
}

.getBlocks = function(x, n) {
  x1 = 0 + (.getRates(x=x, n=n) > 0.5)
  x1 = .completeTails(x1)
  return(x1)
}

.getRates = function(x, n) {
  x0 = as.numeric(filter(x, filter=rep(1,n)/n))
  return(x0)
}

.removeBumps = function(x, thr) {

  x = rle(x)

  px = x$lengths/sum(x$lengths)
  x$values[px<thr] = NA
  validGroups = unique(na.omit(x$values))
  validGroups = setNames(seq_along(validGroups), nm=validGroups)
  validGroups = validGroups[as.character(x$values)]
  names(validGroups) = NULL
  x$values = validGroups
  x = inverse.rle(x)
  x = round(approx(x=x, xout = seq_along(x), rule=2)$y, 0)

  x = rle(x)
  if(length(x$values)<3) {
    x$values = seq_along(x$values)
    x = inverse.rle(x)
    return(x)
  }

  while(any(table(x$values)>1)) {
    bumps = numeric(length(x$values))
    for(i in seq_len(length(x$values)-2)) {
      # if(x$values[i]==x$values[i+2]) x$values[i+1] = x$values[i]
      if(x$values[i]==x$values[i+2]) bumps[i+1] = 1
    }
    thisFirst = which.min(x$lengths/bumps)
    x$values[thisFirst] = x$values[thisFirst+1]
    x = inverse.rle(x)
    x = rle(x)
  }
  # x = inverse.rle(x)
  # x = rle(x)
  x$values = seq_along(x$values)
  x = inverse.rle(x)

  return(x)
}

