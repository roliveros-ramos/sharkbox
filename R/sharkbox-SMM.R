

newSMM = function(G, S, L=3) {

  .generateNames = function(n, label) {
    if(is.character(n)) return(list(n, n))
    sp_code = sprintf("%s_%%0%dd", label, ceiling(log10(n)))
    species_names = sprintf(sp_code, seq_len(n))
    return(list(species_names, species_names))
  }

  if(length(L)==1) L = rep(L, S)
  if(length(L)!=S)
    stop("Number of size groups (L) must be equal to species number (S).")

  species_names = .generateNames(S, "species")
  group_names   = .generateNames(G, "group")

  if(is.character(S)) S = length(S)
  if(is.character(G)) G = length(G)

  out = list()
  out$groups = matrix(nrow=G, ncol=G, dimnames = group_names)
  out$species = list()
  out$size    = list()

  for(i in seq_len(G)) {
    out$species[[i]] = matrix(nrow=S, ncol=S, dimnames = species_names)
  }
  names(out$species) = group_names[[1]]

  for(i in seq_len(S)) {
    size_names = .generateNames(L[i], "size")
    if(L[i]==2) size_names = .generateNames(c("small", "large"))
    if(L[i]==3) size_names = .generateNames(c("small", "medium", "large"))
    out$size[[i]]    = matrix(nrow=L[i], ncol=L[i], dimnames = size_names)
  }
  names(out$size) = species_names[[1]]

  return(out)

}
