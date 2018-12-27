#---- package specific functions

#' @title Check input
#' @author Giovanni Laudanno
#' @description It checks the inputs
#' @inheritParams default_params_doc
#' @return nothing
#' @export
check_input <- function(
  brts_m,
  cond,
  n_0,
  t_0
) {
  if (length(brts_m) <= 0) {
    stop("main clade branching times cannot be an empty vector")
  }
  if (length(brts_s) <= 0) {
    stop("sub clade branching times cannot be an empty vector")
  }
  if (any(brts_m < 0)) {
    stop("all the branching times for the main clade have to be non negative")
  }
  if (any(brts_s < 0)) {
    stop("all the branching times for the sub clade have to be non negative")
  }
  if (!(cond %in% conds())) {
   stop("this conditioning is not implemented")
  }
  if (!(n_0 %in% n_0s())) {
   stop("this n_0 is not implemented")
  }
  if (n_max <= 0) {
   stop("it's not going to work with maximum species set to 0 or less")
  }
}

#' @title Get package name
#' @author Giovanni Laudanno
#' @description Get package name
#' @inheritParams default_params_doc
#' @return Package name
#' @export
get_pkg_name <- function() {
 pkg_name <- "splendid"
 pkg_name
}

#---- general functions
#' @title cat2
#' @author Giovanni Laudanno
#' @description If verbose == TRUE cats the message, otherwise stays silent
#' @inheritParams default_params_doc
#' @return prints on screen
#' @export
cat2 <- function(
  message,
  verbose
) {
 if (verbose == TRUE) {
   cat(message)
 } else {
   return()
 }
}

#' @title Transform parameters
#' @description Transform parameters according to y = x / (1 + x)
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @return transformed parameters
#' @export
pars_transform_forward <- function(pars) {
  pars <- as.numeric(unlist(pars))
  pars_transformed <- pars / (1 + pars)
  pars_transformed[which(pars == Inf)] <- 1
  pars_transformed
}

#' @title Transform parameters back
#' @description Transform parameters back according to x = y / (1 + y)
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @return the original parameters
#' @export
pars_transform_back <- function(pars_transformed) {
  pars_transformed <- as.numeric(unlist(pars_transformed))
  pars <- pars_transformed / (1 - pars_transformed)
  pars
}

#' @title Cut word "loglik" from a name
#' @author Giovanni Laudanno
#' @description Cut word "loglik" from a name
#' @inheritParams default_params_doc
#' @return clean name
#' @export
get_model_name <- function(
  function_name
) {
  if (grepl("loglik", function_name)) {
    model_name <- gsub(
      "_loglik",
      "",
      gsub(
        "loglik_",
        "",
        function_name
      )
    )
  } else {
    model_name <- NA
  }
  model_name
}

#' @title Get function names
#' @author Giovanni Laudanno
#' @description Get function names
#' @inheritParams default_params_doc
#' @return function names
#' @export
get_function_names <- function(
  models
) {
  pkg_name <- get_pkg_name() # nolint internal function
  fun_list <- ls(paste0("package:", pkg_name))
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    pkg_name,
    "!"
  )

  if (is.vector(models)) {
    function_names <- model_names <- which_function <- rep(NA, length(models))
    for (m in seq_along(models)) {
      fun <- eval(models[m])[[1]]
      if (is.character(models[m])) {
        if (length(
          (find_function <- which(fun_list == models[m]))
        ) == 0) {
          stop(error_message)
        }
        which_function[m] <- find_function
      } else {
        for (i in seq_along(fun_list)) {
          if (all.equal(get(fun_list[i]), fun) == TRUE) {
            which_function[m] <- i
          }
        }
      }
      if (is.null(which_function[m]) | is.na(which_function[m])) {
        stop(error_message)
      }
      function_names[m] <- toString(fun_list[which_function[m]])
      model_names[m] <- get_model_name(function_names[m]) # nolint internal function
    }
  } else {
    fun <- eval(models)
    if (is.character(models)) {
      if (length(
        (find_function <- which(fun_list == models))
      ) == 0) {
        stop(error_message)
      }
      which_function <- find_function
    } else {
      which_function <- 0
      for (i in seq_along(fun_list)) {
        if (all.equal(get(fun_list[i]), fun) == TRUE) {
          which_function <- i
        }
      }
    }
    if (is.null(which_function) | is.na(which_function)) {
      stop(error_message)
    }
    function_names <- toString(fun_list[which_function])
    model_names <- get_model_name(function_names) # nolint internal function
  }

  if (any(is.na(model_names))) {
    stop(error_message)
  }
  invisible(function_names)
}

#' @title Check if provided models make sense
#' @author Giovanni Laudanno
#' @description Check if provided models make sense
#' @inheritParams default_params_doc
#' @return models names
#' @export
get_model_names <- function(
  function_names,
  verbose = FALSE
) {
  model_names <- function_names
  for (m in seq_along(function_names)) {
    model_names[m] <- get_model_name(function_names[m]) # nolint internal function
    if (is.null(model_names[m]) | is.na(model_names[m])) {
      stop(paste0(
        "This is not a likelihood function provided by ",
        get_pkg_name(), # nolint internal function
        "!"
      ))
    }
  }
  if (verbose == TRUE) {
    cat("You are using the functions:", model_names)
  }
  model_names
}

#' @title Convert l_0 to brts
#' @author Giovanni Laudanno
#' @description Convert l_0 to brts
#' @inheritParams default_params_doc
#' @return branching times
#' @export
l_2_brts <- function(
 l_0,
 n_0,
 dropextinct = TRUE
) {
 
 n_row_l <- nrow(l_0)
 n_col_l <- ncol(l_0)
 brts <- NULL
 l_0 <- l_0[order(abs(l_0[, 3])), 1:n_col_l]
 dim(l_0) <- c(n_row_l, n_col_l)
 age <- l_0[1, 1]
 l_0[, 1] <- age - l_0[, 1]
 l_0[1, 1] <- -1
 not_min_1 <- l_0[, 4] != -1 & l_0[, 5] == 0
 extant <- abs(l_0[, 3])[!not_min_1]
 if (any(l_0[, 5] > 0)) {
  shifted <- which(l_0[, 5] > 0)
  shift_times <- rep(0, length(extant))
  shift_times[extant == shifted] <- l_0[shifted, 4]
 } else {
  shifted <- 0
  shift_times <- rep(0, length(extant))
 }
 
 l_0[not_min_1, 4] <- age - l_0[not_min_1, 4]
 if (dropextinct == TRUE) {
  sall <- abs(l_0[extant, 3])
  t_end <- age * (sall != shifted) + shift_times * (sall == shifted)
 } else {
  sall <- which(l_0[, 4] >= -1)
  t_end <- (l_0[, 4] == -1) * age + (l_0[, 4] > -1) * l_0[, 4]
 }
 l_0_redux <- l_0[, -4]
 dim(l_0_redux) <- c(n_row_l, n_col_l - 1)
 df <- data.frame(
  matrix(
   l_0_redux[sall, ],
   nrow = length(sall),
   ncol = n_col_l - 1
  )
 )
 lin_list <- cbind(
  df,
  paste0("t", abs(l_0_redux[sall, 3])),
  t_end
 )
 names(lin_list) <- c(
  "birth",
  "parent",
  "id",
  "shift",
  "label",
  "t_end"
 )
 lin_list$label <- as.character(lin_list$label)
 if (nrow(lin_list) == 1) {
  brts <- age
 } else {
  done <- 0
  while (done == 0) {
   j <- which.max(lin_list$birth)
   # daughter <- lin_list$id[j] # nolint
   parent <- lin_list$parent[j]
   parent_j <- which(parent == lin_list$id)
   parent_in_list <- length(parent_j)
   if (parent_in_list == 1) {
    spec1 <- paste0(
     lin_list$label[parent_j],
     ":",
     lin_list$t_end[parent_j] - lin_list$birth[j]
    )
    spec2 <- paste0(
     lin_list$label[j],
     ":",
     lin_list$t_end[j] - lin_list$birth[j]
    )
    lin_list$label[parent_j] <- paste0("(", spec1, ",", spec2, ")")
    lin_list$t_end[parent_j] <- lin_list$birth[j]
    brts <- c(brts, lin_list$birth[j])
    lin_list <- lin_list[-j, ]
   } else {
    lin_list[j, 1:3] <- l_0_redux[which(l_0_redux[, 3] == parent), 1:3]
   }
   if (nrow(lin_list) == 1) {
    done <- 1
   }
  }
  brts <- rev(sort(age - brts))
  if (n_0 == 1) {
   brts <- c(age, brts)
  }
 }
 brts
}

