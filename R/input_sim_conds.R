#' @title Check if inputs for conditioning make sense
#' @description Check if inputs for conditioning make sense
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return nothing
#' @export
sim_cond_check <- function(
 data,
 l_2,
 cond
) {
 if (nrow(l_2) > 2) {
  stop("Currently this only works for 1 shift!")
 }
 if (cond == 2) {
  stop("Cond 2 is not supported!")
 }
}

#' @title Sim conditioning 0
#' @description Sim conditioning 0
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sim_cond_0 <- function(
 data,
 l_2
) {
 1
}

#' @title Sim conditioning 1
#' @description Sim conditioning 1
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sim_cond_1 <- function(
 data,
 l_2
) {
 l_1 <- data$l_1; clade <- 1
 for (clade in l_2$clade_id) {
  l_0 <- l_1[[clade]]
  if (!is.matrix(l_0) && !is.null(l_0)) {
   dim(l_0) <- c(1, 5)
  }
  shifts <- sim_get_shifts_info(l_2 = l_2, clade = clade) # nolint internal function
  shifts_times <- shifts$when
  shifted_id <- l_0[(l_0[, 5] != 0), 3]
  
  t_p <- 0; ts <- shifts_times[1]
  if (clade == 1) {
   # subclade does NOT start from here (M1)
   coords_left  <- sign(l_0[, 3]) == -sign(shifted_id)
   
   #subclade does start from here!!! (M2)
   coords_right <- sign(l_0[, 3]) ==  sign(shifted_id)
   
   l_0_left  <- l_0[coords_left, ]
   dim(l_0_left)  <- c(sum(coords_left), 5)
   l_0_right <- l_0[coords_right, ]
   dim(l_0_right) <- c(sum(coords_right), 5)
   if (nrow(l_0_right) > 0) {
    testit::assert(
     any(l_0_right[, 5] > 0)
    )
   }
   colnames(l_0_left) <- colnames(l_0_right) <- colnames(l_0)
   
   # lineages in M1 born before tshift
   coords_left_cs  <- (l_0_left[, 1]  > shifts_times[1])
   l_0_left_cs  <- l_0_left[coords_left_cs, ]
   dim(l_0_left_cs)  <- c(sum(coords_left_cs), 5)
   
   # lineages in M2 born before tshift
   coords_right_cs <- (l_0_right[, 1] > shifts_times[1])
   l_0_right_cs <- l_0_right[coords_right_cs, ]
   dim(l_0_right_cs) <- c(sum(coords_right_cs), 5)
   
   surv_left_cp  <- sim_check_survival(
    l_0 = l_0_left,
    final_time = t_p
   ) #M1 survives from c to p
   surv_right_cp <- sim_check_survival(
    l_0 = l_0_right,
    final_time = t_p
   ) #M2 survives from c to p
   surv_left_cs  <- sim_check_survival(
    l_0 = l_0_left_cs,
    final_time = ts
   ) #M1 survives from c to s
   surv_right_cs <- sim_check_survival(
    l_0 = l_0_right_cs,
    final_time = ts
   ) #M2 survives from c to s
   
   testit::assert(surv_left_cs  >= surv_left_cp)
   testit::assert(surv_right_cs >= surv_right_cp)
  } else {
   l_0 <- l_1[[clade]]
   surv_s <- sim_check_survival(
    l_0 = l_0,
    final_time = t_p
   ) #S survives from s to p
  }
 }
 
 cond <- (surv_left_cp && surv_right_cs && surv_s)
 cond
}

#' @title Sim conditioning 2
#' @description Sim conditioning 2
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sim_cond_2 <- function(
 data,
 l_2
) {
 l_1 <- data$l_1; clade <- 1
 for (clade in l_2$clade_id) {
  l_0 <- l_1[[clade]]
  if (!is.matrix(l_0) && !is.null(l_0)) {
   dim(l_0) <- c(1, 5)
  }
  shifts <- sim_get_shifts_info(l_2 = l_2, clade = clade) # nolint internal function
  shifts_times <- shifts$when
  shifted_id <- l_0[(l_0[, 5] != 0), 3]
  
  t_p <- 0; ts <- shifts_times[1]
  if (clade == 1) {
   # subclade does NOT start from here (M1)
   coords_left  <- sign(l_0[, 3]) == -sign(shifted_id)
   
   #subclade does start from here!!! (M2)
   coords_right <- sign(l_0[, 3]) ==  sign(shifted_id)
   
   l_0_left  <- l_0[coords_left, ]
   dim(l_0_left)  <- c(sum(coords_left), 5)
   l_0_right <- l_0[coords_right, ]
   dim(l_0_right) <- c(sum(coords_right), 5)
   if (nrow(l_0_right) > 0) {
    testit::assert(
     any(l_0_right[, 5] > 0)
    )
   }
   colnames(l_0_left) <- colnames(l_0_right) <- colnames(l_0)
   
   # lineages in M1 born before tshift
   coords_left_cs  <- (l_0_left[, 1]  > shifts_times[1])
   l_0_left_cs  <- l_0_left[coords_left_cs, ]
   dim(l_0_left_cs)  <- c(sum(coords_left_cs), 5)
   
   # lineages in M2 born before tshift
   coords_right_cs <- (l_0_right[, 1] > shifts_times[1])
   l_0_right_cs <- l_0_right[coords_right_cs, ]
   dim(l_0_right_cs) <- c(sum(coords_right_cs), 5)
   
   surv_left_cp  <- sim_check_survival(
    l_0 = l_0_left,
    final_time = t_p
   ) #M1 survives from c to p
   surv_right_cp <- sim_check_survival(
    l_0 = l_0_right,
    final_time = t_p
   ) #M2 survives from c to p
   surv_left_cs  <- sim_check_survival(
    l_0 = l_0_left_cs,
    final_time = ts
   ) #M1 survives from c to s
   surv_right_cs <- sim_check_survival(
    l_0 = l_0_right_cs,
    final_time = ts
   ) #M2 survives from c to s
   
   testit::assert(surv_left_cs  >= surv_left_cp)
   testit::assert(surv_right_cs >= surv_right_cp)
  } else {
   l_0 <- l_1[[clade]]
   surv_s <- sim_check_survival(
    l_0 = l_0,
    final_time = t_p
   ) #S survives from s to p
  }
 }
 
 cond <- (surv_left_cp && surv_right_cp && surv_s)
 cond
}

#' @title Sim conditioning 3
#' @description Sim conditioning 3
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sim_cond_3 <- function(
 data,
 l_2
) {
 l_1 <- data$l_1; clade <- 1
 for (clade in l_2$clade_id) {
  l_0 <- l_1[[clade]]
  if (!is.matrix(l_0) && !is.null(l_0)) {
   dim(l_0) <- c(1, 5)
  }
  shifts <- sim_get_shifts_info(l_2 = l_2, clade = clade) # nolint internal function
  shifts_times <- shifts$when
  shifted_id <- l_0[(l_0[, 5] != 0), 3]
  
  t_p <- 0; ts <- shifts_times[1]
  if (clade == 1) {
   # subclade does NOT start from here (M1)
   coords_left  <- sign(l_0[, 3]) == -sign(shifted_id)
   
   #subclade does start from here!!! (M2)
   coords_right <- sign(l_0[, 3]) ==  sign(shifted_id)
   
   l_0_left  <- l_0[coords_left, ]
   dim(l_0_left)  <- c(sum(coords_left), 5)
   l_0_right <- l_0[coords_right, ]
   dim(l_0_right) <- c(sum(coords_right), 5)
   if (nrow(l_0_right) > 0) {
    testit::assert(
     any(l_0_right[, 5] > 0)
    )
   }
   colnames(l_0_left) <- colnames(l_0_right) <- colnames(l_0)
   
   # lineages in M1 born before tshift
   coords_left_cs  <- (l_0_left[, 1]  > shifts_times[1])
   l_0_left_cs  <- l_0_left[coords_left_cs, ]
   dim(l_0_left_cs)  <- c(sum(coords_left_cs), 5)
   
   # lineages in M2 born before tshift
   coords_right_cs <- (l_0_right[, 1] > shifts_times[1])
   l_0_right_cs <- l_0_right[coords_right_cs, ]
   dim(l_0_right_cs) <- c(sum(coords_right_cs), 5)
   
   surv_left_cp  <- sim_check_survival(
    l_0 = l_0_left,
    final_time = t_p
   ) #M1 survives from c to p
   surv_right_cp <- sim_check_survival(
    l_0 = l_0_right,
    final_time = t_p
   ) #M2 survives from c to p
   surv_left_cs  <- sim_check_survival(
    l_0 = l_0_left_cs,
    final_time = ts
   ) #M1 survives from c to s
   surv_right_cs <- sim_check_survival(
    l_0 = l_0_right_cs,
    final_time = ts
   ) #M2 survives from c to s
   
   testit::assert(surv_left_cs  >= surv_left_cp)
   testit::assert(surv_right_cs >= surv_right_cp)
  } else {
   l_0 <- l_1[[clade]]
   surv_s <- sim_check_survival(
    l_0 = l_0,
    final_time = t_p
   ) #S survives from s to p
  }
 }
 
 cond <- (surv_left_cp && surv_right_cs && surv_s) ||
  (surv_left_cp && surv_right_cp)
 cond
}