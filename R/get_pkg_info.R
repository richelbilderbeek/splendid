#' @title Conditionings
#' @author Giovanni Laudanno
#' @description Gives the conditionings
#' @inheritParams default_params_doc
#' @return the conditionings
#' @export
get_conds <- function() {
 conds <- c(as.numeric(sub(
  cond_functions(),
  pattern = "cond_function_",
  replacement = ""
 )))
 conds
}

#' @title Starting species
#' @author Giovanni Laudanno
#' @description Gives the amount of starting species
#' @inheritParams default_params_doc
#' @return the possible n_0s
#' @export
get_n_0s <- function() {
 n_0s <- c(2)
 n_0s
}

#' @title Loglik functions
#' @author Giovanni Laudanno
#' @description Gives the available loglik functions
#' @inheritParams default_params_doc
#' @return the loglik function
#' @export
get_logliks <- function() {
 fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
 loglik_functions <- fun_list[sapply(
  fun_list, function(x)
   any(grepl("loglik", x))
 )]
 loglik_functions
}

#' Get the names of the parameters used in the model
#' @author Giovanni Laudanno
#' @export
get_param_names <- function() {
 c("lambda_m", "mu_m", "lambda_s", "mu_s")
}
