#' @title Create parameters
#' @description Create parameters
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return parameters
#' @export
create_pars <- function(
 lambdas,
 mus,
 ks = c(Inf, Inf)
) {
 n_clades <- length(lambdas)
 testit::assert(n_clades > 0)
 testit::assert(length(lambdas) == n_clades)
 testit::assert(length(mus) == n_clades)
 testit::assert(length(ks) == n_clades)
 
 testit::assert(all(lambdas >= 0))
 testit::assert(all(mus >= 0))
 testit::assert(all(ks > 0))
 
 pars <- vector("list", n_clades)
 for (clade in 1:n_clades) {
  pars[[clade]] <- c(
   lambda = lambdas[clade],
   mu = mus[clade],
   K = ks[clade]
  )
  names(pars[[clade]]) <- c(
   paste0(get_param_names(), "_", clade)
  )
 }
 pars
}