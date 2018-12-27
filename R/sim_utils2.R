#' @title Events
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data for the clade
#' @export
sim_clade_events <- function(
 data,
 pars,
 l_2,
 clade,
 events
) {
 if ("shift" %in% colnames(events)) {
  n_shifts <- sum(l_2$motherclade == clade)
  pippo <- colnames(events)
  shift_coord <- which(pippo == "shift")
  pippo_1 <- pippo[seq_along(pippo) < shift_coord]
  pippo_2 <- pippo[seq_along(pippo) > shift_coord]
  event_names <- c(
   pippo_1,
   paste0("shift_", 1:n_shifts),
   pippo_2
  )
  events_1 <- events[pippo_1]
  events_2 <- events[pippo_2]
  events_11 <- matrix(
   unlist(events_1),
   nrow = nrow(events),
   ncol = ncol(events_1)
  )
  events_22 <- matrix(
   unlist(events_2),
   nrow = nrow(events),
   ncol = ncol(events_2)
  )
  priority_coord <- which(rownames(events_1) == "priority")
  priority <- c(
   events_11[priority_coord, ],
   rep(1, n_shifts),
   events_22[priority_coord, ]
  )
  per_capita_coord <- which(rownames(events_1) == "per_capita")
  per_capita <- c(
   events_11[per_capita_coord, ],
   rep(FALSE, n_shifts),
   events_22[per_capita_coord, ]
  )
 } else {
  event_names <- colnames(events)
  priority <- events[rownames(events) == "priority", ]
  per_capita <- events[rownames(events) == "per_capita", ]
 }
 rate <- rate_names <- rep(NA, length(event_names))
 for (i in seq_along(event_names)) {
  if (priority[i] == 2) {
   pippo <- events[event_names[i]]
   rate_names[i] <- levels(droplevels(
    pippo[which(rownames(pippo) == "rate_name"), ]
   ))
   rate[i] <- pars[[clade]][which(names(pars[[clade]]) == rate_names[i])]
  }
 }
 names(rate) <- event_names
 
 shift_times <- l_2$birth_time[l_2$motherclade == clade]
 j <- 1
 time <- rep(NA, length(event_names))
 for (i in seq_along(event_names)) {
  if (priority[i] == 1) {
   if ("shift" == substr(event_names[i], 1, 5)) {
    time[i] <- shift_times[j]
    j <- j + 1
   }
   if ("end" %in% event_names[i]) {
    time[i] <- 0
   }
  }
  if (priority[i] == 2) {
   time[i] <- max(l_2$birth_time)
  }
 }
 total_rate <- length(data$pools[[clade]]) ^ (per_capita == TRUE) * rates
 occurred <- rep(FALSE, length(per_capita))
 clade_events <- rbind(
  rate,
  time,
  priority,
  per_capita,
  total_rate,
  occurred
 )
 colnames(clade_events) <- event_names
 clade_events <- data.frame(t(clade_events))
 clade_events
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
 l_1 <- data$l_1
 pools <- data$pools
 pool <- pools[[clade]]
 pars_clade <- pars[[clade]]
 events <- data$events[[clade]]
 
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
 events <- data$events[[clade]]
 if (nrow(tshifts) > 1) {
  stop("Check the function if you want to implement more than 1 shift!")
 }
 p <- 1; event <- NULL
 while (is.null(event)) {
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
  p <- p + 1
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

#' @title Initialize data for a new clade
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data for the clade
#' @export
sim_initialize_data_new_clade2 <- function(
 data,
 clade,
 pars,
 l_2 = sim_get_standard_l_2(),
 events
) {
 if (clade == 0) {
  return(
   list(
    l_1 = vector("list", length(l_2$clade_id)),
    pools = vector("list", length(l_2$clade_id)),
    n_max = vector("list", length(l_2$clade_id)),
    t = vector("list", length(l_2$clade_id))
   )
  )
 }
 
 clade_events <- sim_clade_events(
  data = data,
  pars = pars,
  l_2 = l_2,
  clade = clade,
  events = events
 )
 
 l_matrix_size <- 1e2
 l_1 <- data$l_1
 n_0 <- sim_get_n_0(l_2 = l_2, clade = clade) # nolint internal function
 t_0 <- sim_get_t_0(l_2 = l_2, clade = clade) # nolint internal function
 motherclade <- sim_get_motherclade(l_2 = l_2, clade = clade) # nolint internal function
 n_max <- n_0
 
 cladeborn <- 1
 if (clade > 1) {
  l_mother <- l_1[[motherclade]]
  motherspecies <- l_mother[l_mother[, 5] == clade, 3]
  cladeborn <- length(motherspecies) != 0
 }
 
 if (cladeborn) {
  t <- t_0
  l_0 <- matrix(0, nrow = l_matrix_size, 5)
  l_0[, 5] <- 0
  l_0[, 4] <- -1
  l_0[, 3] <- 0
  l_0[1, 1:4] <- c(t, 0,  1, -1)
  if (n_0 == 2) {
   l_0[2, 1:4] <- c(t, 1, -2, -1)
  }
  if (clade > 1) {
   l_0[1, 3] <- sign(motherspecies)
  }
  colnames(l_0) <- c(
   "birth_time",
   "parent",
   "id",
   "death_time",
   "shifted_to"
  )
  pool <- l_0[1:n_0, 3]
  data$l_1[[clade]] <- l_0
  data$pools[[clade]] <- pool
  data$n_max[[clade]] <- n_max
  data$events[[clade]] <- clade_events
 } else {
  data$l_1[clade] <- list(NULL)
  data$pools[clade] <- list(NULL)
  data$n_max[clade] <- list(NULL)
  data$events[clade] <- list(NULL)
 }
 
 data$t[[clade]] <- t_0
 return(
  data
 )
}