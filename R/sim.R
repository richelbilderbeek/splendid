#' @title Simulate an sls process
#' @description Simulate an sls process
#' @inheritParams default_params_doc
#' @return l_1 table and brts
#' @author Giovanni Laudanno
#' @export
sim <- function(
 pars,
 cond = 1,
 l_2 = sim_get_standard_l_2()
) {

  # check the parameters

  good_sim <- 0
  while (!good_sim) {

    # initialize data
    data <- sim_initialize_data_new_clade(clade = 0, l_2 = l_2); clade <- 1; # nolint internal function
    for (clade in l_2$clade_id) {

      # initialize data for the clade
      data <- sim_initialize_data_new_clade(
        data = data,
        clade = clade,
        pars = pars,
        l_2 = l_2
      )
      while (data$t[[clade]] > 0) {

        # sample delta_n and delta_t
        deltas <- sim_sample_deltas(
          data = data,
          clade = clade,
          pars = pars
        ); deltas

        # simulate the event
        output <- sim_event(
          data = data,
          clade = clade,
          l_2 = l_2,
          deltas = deltas
        ); output

        data <- output
      }
    }

    # is the simulation in agreement with the conditioning?
    good_sim <- sim_conditioning(
      data = data,
      l_2 = l_2,
      cond = cond
    ); good_sim
  }

  # retrieve branching times info
  brts <- sim_get_brts(
    data = data,
    l_2 = l_2
  )

  return(
    list(
      l_tables = data$l_1,
      brts = brts
    )
  )
}
