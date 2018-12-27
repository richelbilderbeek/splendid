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
 priority <- c(2, 2, 1, 1)
 rate_name <- c("lambda", "mu", "", "")
 per_capita <- c(
  TRUE,
  TRUE,
  FALSE,
  FALSE
 )
 events <- data.frame(rbind(
  priority,
  rate_name,
  per_capita
 ))
 colnames(events) <- event_names
 events
}
