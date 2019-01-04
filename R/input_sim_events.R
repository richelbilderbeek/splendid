#' @title Events
#' @description sim module
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return data for the clade
#' @export
sim_model_events <- function() {
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
 model_events <- data.frame(rbind(
  priority,
  rate_name,
  per_capita
 ), stringsAsFactors = FALSE)
 colnames(model_events) <- event_names
 model_events
}