#' @title Events
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data for the clade
#' @export
sim_events <- function(
 data,
 pars,
 l_2,
 clade = 1
) {
 event_names <- c(
  "speciation",
  "extinction",
  "shift",
  "end"
 )
 rate <- c(
  pars[[clade]][1],
  pars[[clade]][2],
  0,
  0
 )
 time <- c(
  max(l_2$birth_time),
  max(l_2$birth_time),
  l_2$birth_time[l_2$motherclade == clade],
  0
 )
 priority <- c(2, 2, 1, 1)
 per_capita <- c(
  TRUE,
  TRUE,
  FALSE,
  FALSE
 )
 occurred <- rep(FALSE, length(priorities))
 total_rate <- length(data$pools[[clade]]) ^ per_capita * rates
 events <- rbind(
  rates,
  times,
  priorities,
  per_capita,
  total_rate,
  occurred
 )
 colnames(events) <- event_names
 events <- data.frame(t(events))
 events
}

#' @title Sample deltas for Doob-Gillespie
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return delta_t
#' @export
sim_sample_delta_t <- function(
 data,
 clade,
 pars,
 l_2
) {
 events <- sim_events(
  data = data,
  pars = pars,
  l_2 = l_2,
  clade = clade
 )
 l_1 <- data$l_1
 pools <- data$pools
 pool <- pools[[clade]]
 pars_clade <- pars[[clade]]
 
 # n <- length(pool)
 total_rate <- sum(events$total_rate)
 testit::assert(total_rate >= 0)
 if (total_rate > 0) {
  delta_t <- (total_rate > 0) *
   stats::rexp(1, rate = total_rate + (total_rate == 0)) +
   (total_rate == 0) * l_1[[1]][1, 1]
 } else {
  delta_t <- 1e5
 }
 delta_t <- unname(delta_t)
 delta_t
}

#' @title Determines the event
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the event
#' @export
sim_decide_event2 <- function(
 data,
 clade,
 delta_t,
 l_2
) {
 t <- data$t[[clade]]
 l_1 <- data$l_1
 l_0 <- l_1[[clade]]
 already_shifted <- any(l_0[, 5] > 0)
 tshifts <- sim_get_shifts_info(l_2 = l_2, clade = clade)
 if (nrow(tshifts) > 1) {
  stop("Check the function if you want to implement more than 1 shift!")
 }
 p <- 0; event <- NULL
 while (is.null(event)) {
  p <- p + 1
  sub_events <- events[events$priorities == p, ]
  not_occurred_yet <- sub_events[sub_events$occurred == 0, ]
  earliest <- not_occurred_yet[
   not_occurred_yet$times == max(not_occurred_yet$times) &
    t - delta_t < not_occurred_yet$times, ]
  if (nrow(earliest) > 0) {
   event <- sample(
    rownames(earliest),
    size = 1,
    prob = earliest$total_rate
   )
   if (p == 1) {
    events$occurred <- events$occurred + (rownames(events) == event)
   }
  }
 }
 testit::assert(!is.null(event))
 event
}

#' @title Update data and time given the event
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data
#' @export
sim_use_event2 <- function(
 data,
 clade,
 l_2,
 event,
 delta_t
) {
 
 shifts <- sim_get_shifts_info(l_2 = l_2, clade = clade) # nolint internal function
 t <- data$t[[clade]] - delta_t
 l_0 <- data$l_1[[clade]]
 pool <- data$pools[[clade]]; pool
 n <- length(pool)
 n_max <- data$n_max[[clade]]
 shifted <- 0
 
 if (event == "shift") {
  where <- shifts$where
  t <- shifts$when #time becomes the shift point
  if (n > 1) {
   shifted <- sample(pool, replace = FALSE, size = 1)
  } else {
   shifted <- pool
  }
  l_0[abs(shifted), 4] <- t #remove the shifted from the l_0 table
  pool <- pool[pool != shifted]
  l_0[abs(l_0[, 3]) == abs(shifted), 5] <- where #register if shift occurs
 }
 
 if (event == "speciation") {
  data <- sim_adapt_l_matrix_size(
   data = data,
   clade = clade
  )
  l_0 <- data$l_1[[clade]]
  pool <- data$pools[[clade]]; pool
  n <- length(pool)
  
  if (n > 1) {
   parents <- sample(pool, replace = FALSE, size = 1)
  } else {
   parents <- pool
  }
  n_max <- n_max + 1
  
  new_line <- c(
   t,
   parents,
   abs(n_max) * sign(parents),
   -1,
   0
  )
  dim(new_line) <- c(1, 5)
  l_0[n_max, ] <- new_line
  pool <- c(pool, abs(n_max) * sign(parents))
 }
 
 if (event == "extinction") {
  if (n > 1) {
   dead <- sample(pool, replace = FALSE, size = 1)
  } else {
   dead <- pool
  }
  l_0[abs(dead), 4] <- t
  pool <- pool[pool != dead]
 }
 
 if (event == "end" | length(pool) == 0) {
  t <- 0
  l_02 <- sim_cut_l_matrix(l_0) # nolint internal function
  l_0 <- l_02
 }
 
 # store output
 t <- unname(t)
 data2 <- data
 if (is.null(l_0)) {
  data2$l_1[clade] <- list(l_0)
 } else {
  data2$l_1[[clade]] <- l_0
 }
 if (is.null(pool)) {
  data2$pools[clade] <- list(pool)
 } else {
  data2$pools[[clade]] <- pool
 }
 if (is.null(pool)) {
  data2$n_max[clade] <- list(n_max)
 } else {
  data2$n_max[[clade]] <- n_max
 }
 data2$t[[clade]] <- t
 
 return(
  data = data2
 )
}

#' @title Apply the event to the simulation
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the output of the event
#' @export
sim_event2 <- function(
 data,
 clade,
 l_2,
 delta_t
) {
 # decide the event
 event <- sim_decide_event2(
  data = data,
  clade = clade,
  l_2 = l_2,
  delta_t = delta_t
 ); event
 
 # modify data accordingly
 output <- sim_use_event2(
  data = data,
  clade = clade,
  l_2 = l_2,
  event = event,
  delta_t = delta_t
 ); output
 output
}