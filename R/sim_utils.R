# MAIN COMPONENTS ----

#' @title Creates the standard l_2
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return l_2
#' @export
sim_get_standard_l_2 <- function(
 crown_age = 10,
 n_0 = 2,
 shift_time = 4
) {
 testit::assert(crown_age > shift_time)
 l_2 <- as.data.frame(matrix(0, nrow = 2, ncol = 4))
 l_2[, 1] <- c(0, shift_time)
 l_2[, 2] <- c(0, 1)
 l_2[, 3] <- c(1, 2)
 l_2[, 4] <- c(n_0, 1)
 l_2[1, 1] <- crown_age
 colnames(l_2) <- c("birth_time", "motherclade", "clade_id", "n_0")
 
 t_0s <- l_2[, 1]
 motherclades <- l_2[, 2]
 n_0s <- l_2[, 4]
 testit::assert(
  all(motherclades < seq(from = 1, to = length(motherclades)))
 )
 testit::assert(all(n_0s >= 1))
 if (any(n_0s[-1] != 1)) {
  stop("Every subclade should start with 1 species!")
 }
 testit::assert(all(t_0s > 0))
 
 l_2
}

#' @title Initialize data for a new clade
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data for the clade
#' @export
sim_initialize_data_new_clade <- function(
 data,
 clade,
 pars,
 l_2 = sim_get_standard_l_2(),
 model_events
) {
 if (clade == 0) {
  return(
   list(
    l_1 = vector("list", length(l_2$clade_id)),
    pools = vector("list", length(l_2$clade_id)),
    n_max = vector("list", length(l_2$clade_id)),
    events = vector("list", length(l_2$clade_id)),
    t = vector("list", length(l_2$clade_id))
   )
  )
 }
 
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
  data$events[[clade]] <- sim_clade_events(
   data = data,
   pars = pars,
   l_2 = l_2,
   clade = clade,
   model_events = model_events
  )
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
 model_events
) {
 if ("shift" %in% colnames(model_events)) {
  n_shifts <- sum(l_2$motherclade == clade)
  pippo <- colnames(model_events)
  shift_coord <- which(pippo == "shift")
  pippo_1 <- pippo[seq_along(pippo) < shift_coord]
  pippo_2 <- pippo[seq_along(pippo) > shift_coord]
  event_names <- c(
   pippo_1,
   paste0("shift_", 1:n_shifts),
   pippo_2
  )
  events_1 <- model_events[pippo_1]
  events_2 <- model_events[pippo_2]
  events_11 <- matrix(
   unlist(events_1),
   nrow = nrow(model_events),
   ncol = ncol(events_1)
  )
  events_22 <- matrix(
   unlist(events_2),
   nrow = nrow(model_events),
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
  event_names <- colnames(model_events)
  priority <- model_events[rownames(model_events) == "priority", ]
  per_capita <- model_events[rownames(model_events) == "per_capita", ]
 }
 rate <- rate_names <- rep(0, length(event_names))
 for (i in seq_along(event_names)) {
  if (priority[i] == 2) {
   pippo <- model_events[event_names[i]]
   rate_names[i] <- levels(droplevels(
    pippo[which(rownames(pippo) == "rate_name"), ]
   ))
   rate[i] <- pars[[clade]][
    which(
     grepl(
      rate_names[i], names(pars[[clade]])
     )
    )
    ]
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
 total_rate <- length(data$pools[[clade]]) ^ (per_capita == TRUE) * rate # PN: still not working, but fixed name of rates to rate
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
 clade_events <- data.frame(
  t(clade_events),
  stringsAsFactors = FALSE
 )
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

 total_rate <- sum(as.numeric(events$total_rate), na.rm = TRUE)
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
sim_decide_event <- function(
 data,
 clade,
 delta_t,
 l_2
) {
 t <- data$t[[clade]]
 l_1 <- data$l_1
 l_0 <- l_1[[clade]]
 events <- data$events[[clade]]
 if (nrow(l_2) > 2) {
  stop("Check the function if you want to implement more than 1 shift!")
 }
 p_level <- 1; event <- NULL
 while (is.null(event)) {
  p_level_events <- events[events$priority == p_level, ]
  not_occurred_yet <- p_level_events[p_level_events$occurred == FALSE, ]
  times <- as.numeric(as.matrix(not_occurred_yet$time))
  earliest <- not_occurred_yet[
   times == max(times) &
    t - delta_t < times,
   ]
  if (nrow(earliest) > 0) {
   probability <- as.numeric(as.matrix(earliest$total_rate))
   if (any(is.na(probability))) {
    probability <- rep(1, length(probability))
    if (p_level > 1) {
     print("there is a NA probability!")
    }
   }
   event <- sample(
    rownames(earliest),
    size = 1,
    prob = probability
   )
  }
  p_level <- p_level + 1
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
sim_use_event <- function(
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

 # store event with priority 1
 temp <- data$events[[clade]]
 temp_matrix <- as.matrix(temp)
 if (temp_matrix[rownames(temp_matrix) == event, "priority"] == "1") {
  temp_matrix[rownames(temp_matrix) == event, "occurred"] <- "TRUE"
 }
 data$events[[clade]] <- as.data.frame(
  temp_matrix,
  stringsAsFactors = FALSE
 )
 
 # update total rates
 total_rates <- rep(length(data$pools[[clade]]), length(data$events[[clade]]$per_capita)) ^ 
  as.numeric(as.logical(data$events[[clade]]$per_capita)) * 
  as.numeric(data$events[[clade]]$rate)
 data$events[[clade]]$total_rate <- total_rates
 
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
sim_event <- function(
 data,
 clade,
 l_2,
 delta_t
) {
 # decide the event
 event <- sim_decide_event(
  data = data,
  clade = clade,
  l_2 = l_2,
  delta_t = delta_t
 ); print(event)
 
 # modify data accordingly
 output <- sim_use_event(
  data = data,
  clade = clade,
  l_2 = l_2,
  event = event,
  delta_t = delta_t
 ); output
 output
}

#' @title Check if the conditioning is fulfilled
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a boolean
#' @export
sim_conditioning <- function(
 data,
 l_2,
 cond
) {
 
 #check
 sim_cond_check(
  data = data,
  l_2 = l_2,
  cond = cond
 )
 
 #conditioning
 keep_the_sim <- get(cond_functions()[cond])(
  data = data,
  l_2 = l_2
 )
 keep_the_sim
}

# UTILITIES ----

#' @title Initialize n for a new clade
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return n0 for the clade
#' @export
sim_get_n_0 <- function(
 l_2 = sim_get_standard_l_2(),
 clade
) {
 n_0 <- l_2[l_2[, 3] == clade, 4]
 n_0
}

#' @title Initialize t for a new clade
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return t_0 for the clade
#' @export
sim_get_t_0 <- function(
 l_2 = sim_get_standard_l_2(),
 clade
) {
 t_0 <- l_2[l_2[, 3] == clade, 1]
 t_0
}

#' @title Get the motherclade of a given clade
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return the motherclade
#' @export
sim_get_motherclade <- function(
 l_2 = sim_get_standard_l_2(),
 clade
) {
 motherclade <- l_2[l_2[, 3] == clade, 2]
 motherclade
}

#' @title Get shifts' information
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return when and where the shifts occurs
#' @export
sim_get_shifts_info <- function(
 l_2 = sim_get_standard_l_2(),
 clade
) {
 when  <- l_2[l_2[, 2] == clade, 1]
 where <- l_2[l_2[, 2] == clade, 3]
 
 info <- data.frame(when = when, where = where)
 info
}

#' @title Increase the size of the l_0 table if needed
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a bigger l_0 table
#' @export
sim_adapt_l_matrix_size <- function(
 data,
 clade
) {
 n_max <- data$n_max[[clade]]
 l_0 <- data$l_1[[clade]]
 if (n_max >= nrow(l_0) - 2) {
  append_l_0 <- matrix(0, nrow = nrow(l_0), ncol = ncol(l_0))
  append_l_0[, 4] <- -1
  l_02 <- rbind(l_0, append_l_0)
  data$l_1[[clade]] <- l_02
 }
 return(data)
}

#' @title Cuts the l_0 table only to the useful rows
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return a clean l_0 table
#' @export
sim_cut_l_matrix <- function(
 l_0
) {
 if (is.null(l_0)) {
  return(NULL)
 }
 if (is.matrix(l_0)) {
  n_max <- max(abs(l_0[, 3]))
  l_02 <- l_0[1:n_max, ]
  cols_l_0 <- ncol(l_0)
 } else {
  n_max <- max(abs(l_0[3]))
  l_02 <- l_0
  cols_l_0 <- length(l_0)
 }
 dim(l_02) <- c(n_max, cols_l_0)
 return(l_02)
}
