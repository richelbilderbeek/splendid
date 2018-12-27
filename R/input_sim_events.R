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
 event_rates <- c(
  pars[[clade]][1],
  pars[[clade]][2],
  0,
  0
 )
 event_times <- c(
  -1,
  -1,
  l_2$birth_time[l_2$motherclade == clade],
  0
 )
 event_priorities <- c(3, 3, 1, 2)
 events <- rbind(
  event_rates,
  event_times,
  event_priorities
 )
 colnames(events) <- event_names
 events
}