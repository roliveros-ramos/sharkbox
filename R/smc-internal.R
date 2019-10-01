#' @export
predict.multinomial.TM = function(par, mf) {

  N = nrow(mf)
  initial = NULL
  states  = dimnames(par)[[1]]

  if(identical(attr(par, "by"), "NA"))
    return(states[.sim_multinom(N=N, A=par)])

  fac_by = .getFactorLevels(mf=mf, vars=attr(par, "by"), asCharacter = TRUE)
  Ns  = table(fac_by)
  lev = names(Ns)

  out = numeric(N)

  for(i in seq_along(Ns)) {
    out[fac_by==lev[i]] = states[.sim_multinom(N=Ns[i], A=par[ ,lev[i]])]
  }

  return(out)
}


#' @export
predict.mc.TM = function(par, mf) {

  N = nrow(mf)
  initial = attr(par, "initial")
  states  = names(initial)

  if(identical(attr(par, "by"), "NA"))
    return(states[.sim_mc(N=N, A=par, initial=initial)])

  fac_by = .getFactorLevels(mf=mf, vars=attr(par, "by"), asCharacter = TRUE)
  Ns  = table(fac_by)
  lev = names(Ns)

  out = numeric(N)

  for(i in seq_along(Ns)) {
    out[fac_by==lev[i]] = states[.sim_mc(N=Ns[i], A=par[, , lev[i]], initial=initial)]
  }

  return(out)
}

#' @export
predict.block.TM = function(object, mf) {

  N = nrow(mf)
  initial = object$prop
  states  = names(initial)

  A       = object$TM(N) # TM for size n

  out = states[.sim_mc(N=N, A=A, initial=initial)]

  return(out)
}

.sim_mc = function(N, A, initial) {
  out     = numeric(N)
  states  = seq_len(nrow(A))
  out[1] = sample(states, size=1, prob = initial)
  for(i in seq_len(N-1)) {
    prob = if(sum(A[out[i],])==0) initial else A[out[i],]
    out[i+1] = sample(states, size=1, prob = prob)
  }
  return(out)
}

.sim_multinom = function(N, A) {
  return(sample(x=length(A), size=N, prob=A, replace = TRUE))
}

# what to do if rows sum zero? always must sum one! Uniform? Initial?
