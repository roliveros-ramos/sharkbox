#' @export
interpret.smc = function (formula, textra = NULL, extra.special = NULL) {
  p.env = environment(formula)
  tf = terms.formula(formula, specials = c("b", "mc",
                                       "m", extra.special))
  terms = attr(tf, "term.labels")
  attr(tf, "intercept") = 1 # no intercept
  nt = length(terms)
  if (attr(tf, "response") > 0) {
    response = as.character(attr(tf, "variables")[2])
    if(grepl(x=response, pattern="\\|")) {
      response = gsub(x=response, pattern=" ", replacement = "")
      xr = unlist(strsplit(x=response, split="\\|"))
      response = xr[1]
      splits = xr[-1]
      if(length(splits)>1) stop("Only one split factor is currently allowed.")
    } else splits = NULL
  }
  else {
    response = NULL
  }
  doSplit = length(splits)>0
  sp = attr(tf, "specials")$b
  tp = attr(tf, "specials")$mc
  tip = attr(tf, "specials")$m
  zp = if (is.null(extra.special))
    NULL
  else attr(tf, "specials")[[extra.special]]
  off = attr(tf, "offset")
  vtab = attr(tf, "factors")
  if (length(sp) > 0)
    for (i in 1:length(sp)) {
      ind = (1:nt)[as.logical(vtab[sp[i], ])]
      sp[i] = ind
    }
  if (length(tp) > 0)
    for (i in 1:length(tp)) {
      ind = (1:nt)[as.logical(vtab[tp[i], ])]
      tp[i] = ind
    }
  if (length(tip) > 0)
    for (i in 1:length(tip)) {
      ind = (1:nt)[as.logical(vtab[tip[i], ])]
      tip[i] = ind
    }
  if (length(zp) > 0)
    for (i in 1:length(zp)) {
      ind = (1:nt)[as.logical(vtab[zp[i], ])]
      zp[i] = ind
    }
  k = kt = kti = kt2 = ks = kz = kp = 1
  len.sp = length(sp)
  len.tp = length(tp)
  len.tip = length(tip)
  len.zp = length(zp)
  ns = len.sp + len.tp + len.tip + len.zp
  pav = av = rep("", 0)
  smooth.spec = list()
  sbns = loadNamespace("sharkbox")
  if(nt) {
    for (i in seq_len(nt)) {
      if (k <= ns && ((ks <= len.sp && sp[ks] == i) ||
                      (kt <= len.tp && tp[kt] == i) || (kz <= len.zp &&
                                                        zp[kz] == i) || (kti <= len.tip && tip[kti] ==
                                                                         i))) {
        st = try(eval(parse(text = paste("sharkbox::",
                                         terms[i], sep = "")), envir = p.env),
                 silent = TRUE)
        if (inherits(st, "try-error"))
          st = eval(parse(text = terms[i]), enclos = p.env,
                    envir = sbns)
        if (!is.null(textra)) {
          pos = regexpr("(", st$lab, fixed = TRUE)[1]
          st$label = paste(substr(st$label, start = 1,
                                  stop = pos - 1), textra, substr(st$label,
                                                                  start = pos, stop = nchar(st$label)), sep = "")
        }
        smooth.spec[[k]] = st
        if (ks <= len.sp && sp[ks] == i)
          ks = ks + 1
        else if (kt <= len.tp && tp[kt] == i)
          kt = kt + 1
        else if (kti <= len.tip && tip[kti] == i)
          kti = kti + 1
        else kz = kz + 1
        k = k + 1
      }
      else {
        av[kp] = terms[i]
        kp = kp + 1
      }
    }
  }

  if (!is.null(off)) {
    av[kp] = as.character(attr(tf, "variables")[1 +
                                                   off])
    kp = kp + 1
  }
  pf = paste(response, "~", paste(av, collapse = " + "))
  if (attr(tf, "intercept") == 0) {
    pf = paste(pf, "-1", sep = "")
    if (kp > 1)
      pfok = 1
    else pfok = 0
  }
  else {
    pfok = 1
    if (kp == 1) {
      pf = paste(pf, "1")
    }
  }
  fake.formula = pf
  if (length(smooth.spec) > 0)
    for (i in 1:length(smooth.spec)) {
      nt = length(smooth.spec[[i]]$term)
      ff1 = paste(smooth.spec[[i]]$term[1:nt], collapse = "+")
      fake.formula = paste(fake.formula, "+", ff1)
      if (smooth.spec[[i]]$by != "NA") {
        fake.formula = paste(fake.formula, "+",
                              smooth.spec[[i]]$by)
        av = c(av, smooth.spec[[i]]$term, smooth.spec[[i]]$by)
      }
      else av = c(av, smooth.spec[[i]]$term)
    }
  if(length(splits)>0) {
    for(i in seq_along(splits)) {
      fake.formula = paste(fake.formula, "+", splits[i])
    }
  }
  fake.formula = as.formula(fake.formula, p.env)
  if(length(av)) {
    pred.formula = as.formula(paste("~", paste(av,
                                                collapse = "+")))
    pav = all.vars(pred.formula)
    pred.formula = reformulate(pav)
  }
  else pred.formula = ~ 1
  if(length(splits)) {
    split.formula = as.formula(paste("~", paste(splits,
                                               collapse = "+")))
    sav = all.vars(split.formula)
    split.formula = reformulate(sav)
  } else split.formula = ~ 1

  sort.formula = reformulate(c(response, splits))

  xsmooths = sapply(smooth.spec, class)
  if(sum(xsmooths=="block.spec")>1)
    stop("Only one block specification is allowed.")

  names(smooth.spec) = sapply(smooth.spec, FUN="[[", i="label")

  ret = list(pf = as.formula(pf, p.env), split.names=splits, split.formula=split.formula,
             pfok = pfok, smooth.spec = smooth.spec, fake.formula = fake.formula,
             response = response, fake.names = c(av, splits), pred.names = pav,
             pred.formula = pred.formula, sort.formula=sort.formula,
             doSplit=doSplit)
  class(ret) = "split.smc.formula"
  return(ret)
}
