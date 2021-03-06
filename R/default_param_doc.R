#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param age the age of the phylogeny
#' @param brts branchin times
#' @param brts_m branching times for the Main-clade
#' @param brts_s branching times for the Sub-clade
#' @param brts_clade branching times for a single clade
#' @param clade the id of the clade
#' @param cond type of conditioning:
#' \itemize{
#'   \item cond = 0 no conditiong;
#'   \item cond = 1 conditions on the survival of crown descendents;
#'   \item cond = 2 not available;
#'   \item cond = 3 conditions on the survival of subclade and
#'   on the other crown descendents in the main clade;
#'   \item cond = 4 conditions on the survival of the subclade and
#'   both crown descendents in the main clade;
#' }
#' @param crown_age the age of the phylogeny
#' @param data contains all the information about the simulated process
#' @param data_folder The data folder insider the project folder.
#' @param deltas in the Doob-Gillespie algorithm,
#' the collection of delta_n and delta_t, which are, respectively,
#' the change in number of species and
#' the waiting time for the next event to occur
#' @param delta_n in the Doob-Gillespie algorithm,
#' the change in number of species
#' @param delta_t in the Doob-Gillespie algorithm,
#' the waiting time for the next event to occur
#' @param dropextinct TRUE if you want to remove the dead species from the
#' tree. FALSE otherwise.
#' @param d_0 starting value for BiSSE's D function
#' @param d_0s starting values for BiSSE's D functions
#' @param e_0 starting value for BiSSE's E function
#' @param event the event occurring in the simulated process at a given time
#' @param final_time the final time that you want to consider for the survival
#' of the species considered in the l table
#' @param fun a function
#' @param fun1 a function
#' @param fun2 another function
#' @param function_name function name
#' @param function_names function names
#' @param functions_names function names
#' @param k frequencies in the Discrete Fourier Transform (DFT)
#' @param ks carrying capacities, for all the clades
#' @param l_0 the l table, for a single clade
#' @param l_1 the collection of all the l tables, for all the clades
#' @param l_2 the matrix containing the information about how the subclades are
#' nested into the main clade. See sls_sim.get_standard_l_2() for more info.
#' @param l_matrix the l table
#' @param l_matrix_size the initial length of the l matrix. It will be
#' increased if necessary
#' @param lambda speciation rate
#' @param lambdas speciation rates, for all the clades
#' @param lambdaterms set it to TRUE if you desire the powers of lambda
#' in the likelihood
#' @param log_scale set it to TRUE if you desire the output in log form
#' @param loglik_function the loglik function you want to use
#' @param lx size of the matrix
#' @param matrix_size size of the matrix
#' @param message the message to print
#' @param models the models you want to use to define the likelihood
#' @param mu extinction rate
#' @param mus extinction rate, for all the clades
#' @param n number of lineages
#' @param n_0 starting number of lineages
#' @param n_0s starting number of lineages for all the clades
#' @param n_max maximum number of lineages to consider
#' @param optim_ids ids of the parameters you want to optimize.
#' @param output the mle output
#' @param pars parameters of the likelihood functions:
#' \itemize{
#'   \item pars[1] is lambda_m, i.e. speciation rate of the main clade;
#'   \item pars[2] is mu_m, i.e. extinction rate of the main clade;
#'   \item pars[3] is lambda_s, i.e. speciation rate of the sub clade;
#'   \item pars[4] is mu_s, i.e. extinction rate of the sub clade;
#' }
#' @param pars_m parameters for the main clade (lambda, mu)
#' @param pars_s parameters for the sub clade (lambda, mu)
#' @param pars_clade parameters for a single clade
#' @param pars_transformed parameters of the likelihood functions, transformed
#' according to y = x / (1 + x)
#' @param project_folder the folder when you want to save data and results
#' @param results mle results
#' @param results_folder The results folder insider the project folder.
#' @param seed the seed
#' @param sim the results of a sim run
#' @param t time
#' @param t_0 starting time
#' @param t_0s starting time for each clade
#' @param tf ending time
#' @param t_f ending time
#' @param t_c crown time
#' @param td decoupling time
#' @param t_d decoupling time
#' @param tds decoupling times
#' @param t_ds decoupling times
#' @param t_p present time
#' @param ts times
#' @param times times
#' @param tbar time left from shift time to the present
#' @param true_pars true parameter values when running the ml process.
#' @param shift_time the time of the shift
#' @param sim_pars parameters of the simulation
#' @param verbose set it to TRUE if you want to see the outputs on screen
#' @param start_pars parameters to start from for the search of the likelihood
#' maximum
#' @param vec a vector or a matrix to be transformed
#' @param missnumspec number of missing (unseen) species in the phylogeny
#' @param ddep see DDD package
#' @param trparsopt see DDD package
#' @param trparsfix see DDD package
#' @param idparsopt see DDD package
#' @param idparsfix see DDD package
#' @param idparsnoshift see DDD package
#' @param initparsopt see DDD package
#' @param optimmethod see DDD package
#' @param tolerance see DDD package
#' @param pars2 see DDD package

default_params_doc <- function(
  age,
  brts,
  brts_m,
  brts_s,
  brts_clade,
  cond,
  clade,
  crown_age,
  d_0,
  d_0s,
  data,
  data_folder,
  delta_n,
  delta_t,
  deltas,
  dropextinct,
  e_0,
  event,
  final_time,
  fun,
  fun1,
  fun2,
  function_name,
  function_names,
  functions_names,
  k,
  lambda,
  lambdas,
  mu,
  mus,
  ks,
  l_0,
  l_1,
  l_2,
  l_matrix_size,
  l_matrix,
  log_scale,
  lambdaterms,
  loglik_function,
  lx,
  matrix_size,
  message,
  missnumspec,
  models,
  n,
  n_0,
  n_0s,
  n_max,
  optim_ids,
  output,
  pars,
  pars_m,
  pars_s,
  pars_clade,
  pars_transformed,
  project_folder,
  results,
  results_folder,
  shift_time,
  sim,
  sim_pars,
  start_pars,
  t,
  t_c,
  t_d,
  ts,
  t_0,
  t_0s,
  tf,
  t_f,
  td,
  t_p,
  tds,
  t_ds,
  times,
  tbar,
  true_pars,
  vec,
  verbose,
  ddep,
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  idparsnoshift,
  initparsopt,
  optimmethod,
  tolerance,
  pars2,
  seed
) {
  # Nothing
}
