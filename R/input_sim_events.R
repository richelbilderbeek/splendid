#' @title Events
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return events that can happen in the phylogeny
#' @export
sim_events <- function() {
 event_names <- c(
  "speciation",
  "extinction",
  "shift",
  "end"
 )
 event_priorities <- c(2, 2, 1, 1)
 event_rate_names <- c("lambda", "mu", "", "")
 events <- rbind(
  event_priorities,
  event_rate_names
 )
 colnames(events) <- event_names
 events
}
