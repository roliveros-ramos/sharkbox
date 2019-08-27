.makeLabel = function(term, special) {
  dim = length(term)
  full.call <- paste0(special, "(", term[1])
  if (dim > 1)
    for (i in 2:dim) full.call <- paste0(full.call, ",", term[i])
  label <- paste0(full.call, ")")
  return(label)
}

mc = function(..., by = NA) {
  vars   = as.list(substitute(list(...)))[-1]
  dim    = length(vars)
  by.var = deparse(substitute(by), backtick = TRUE)
  term   = character(dim)
  for(i in 1:dim) {
    term[i] = attr(terms(reformulate(deparse(vars[[i]], backtick = TRUE))),
                    "term.labels")
  }
  if(length(unique(term)) != dim)
    stop("Repeated variables as arguments are not permitted.")
  label = .makeLabel(term, "mc")
  out = list(term = term, by = by.var, label = label, dim = dim)
  class(out) = "mc.spec"
  return(out)
}

m = function(..., by = NA) {
  vars   = as.list(substitute(list(...)))[-1]
  dim    = length(vars)
  by.var = deparse(substitute(by), backtick = TRUE)
  term   = character(dim)
  for(i in 1:dim) {
    term[i] = attr(terms(reformulate(deparse(vars[[i]], backtick = TRUE))),
                   "term.labels")
  }
  if(length(unique(term)) != dim)
    stop("Repeated variables as arguments are not permitted.")
  label = .makeLabel(term, "m")
  out = list(term = term, by = by.var, label = label, dim = dim)
  class(out) = "multinomial.spec"
  return(out)
}

b = function(...) {
  vars   = as.list(substitute(list(...)))[-1]
  dim    = length(vars)
  term   = character(dim)
  for(i in 1:dim) {
    term[i] = attr(terms(reformulate(deparse(vars[[i]], backtick = TRUE))),
                   "term.labels")
  }
  if(length(unique(term)) != dim)
    stop("Repeated variables as arguments are not permitted.")
  label = .makeLabel(term, "b")
  out = list(term = term, by = "NA", label = label, dim = dim)
  class(out) = "block.spec"
  return(out)
}
