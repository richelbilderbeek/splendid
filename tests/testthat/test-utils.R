context("sls_utils")

test_that("p_transition_matrix", {
  lambda <- 0.5
  mu <- 0.1
  matrix_size <- 20
  test <- p_transition_matrix(
    lambda = lambda,
    mu = mu,
    matrix_size = matrix_size
  )
  testthat::expect_true(
    all(dim(test) == rep(matrix_size, 2))
  )
})

test_that("cat2", {
  testthat::expect_output(
    cat2(
      message = "ciao",
      verbose = TRUE
    )
  )
  testthat::expect_silent(
    cat2(
      message = "ciao",
      verbose = FALSE
    )
  )
})

test_that("check_input", {

  brts_m <- c(6, 3, 2)
  brts_s <- c(2.5, 1)
  start_pars <- c(0.5, 0.3, 0.5, 0.3)
  cond <- 3
  n_0 <- 2
  n_max <- 1e2

  # use
  testthat::expect_silent(
    check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    )
  )

  #abuse
  testthat::expect_error(
    test <- check_input(
      brts_m = c(),
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    ),
    "main clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- check_input(
      brts_m = brts_m,
      brts_s = c(),
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    ),
    "sub clade branching times cannot be an empty vector"
  )
  testthat::expect_error(
    test <- check_input(
      brts_m = c(-2, -1),
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    ),
    "all the branching times for the main clade have to be non negative"
  )
  testthat::expect_error(
    test <- check_input(
      brts_m = brts_m,
      brts_s = c(-2, -1),
      cond = cond,
      n_0 = n_0,
      n_max = n_max
    ),
    "all the branching times for the sub clade have to be non negative"
  )
  testthat::expect_error(
    test <- check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = n_0,
      n_max = -1
    ),
    "it's not going to work with maximum species set to 0 or less"
  )
  testthat::expect_error(
    test <- check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = 17,
      n_0 = n_0,
      n_max = n_max
    ),
    "this conditioning is not implemented"
  )
  testthat::expect_error(
    test <- check_input(
      brts_m = brts_m,
      brts_s = brts_s,
      cond = cond,
      n_0 = 10,
      n_max = n_max
    ),
    "this n_0 is not implemented"
  )
})

test_that("conds", {
  testthat::expect_true(
    is.numeric(conds())
  )
})

test_that("n_0s", {
  testthat::expect_true(
    is.numeric(n_0s())
  )
})

test_that("get_pkg_name", {
  testthat::expect_true(
    get_pkg_name() == "sls"
  )
})

test_that("get_function_names & get_model_names", {
  #use
  testthat::expect_true(
    length(
      get_function_names(
        models = get_logliks()
      )
    ) > 0
  )
  testthat::expect_silent(
    get_model_names(
      function_names = get_logliks(),
      verbose = FALSE
    )
  )
  testthat::expect_output(
    get_model_names(
      function_names = get_logliks(),
      verbose = TRUE
    ),
    "You are using the functions: sls_p sls_q"
  )
  #abuse
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    get_pkg_name(),
    "!"
  )
  testthat::expect_error(
    get_function_names(
      models = "nonsense"
    ),
    error_message
  )
  testthat::expect_error(
    get_function_names(
      models = c("nonsense1", "nonsense2")
    ),
    error_message
  )
  testthat::expect_error(
    get_function_names(
      models = grepl
    ),
    error_message
  )
  testthat::expect_error(
    get_function_names(
      models = c(exp, grepl)
    ),
    error_message
  )
  testthat::expect_error(
    get_model_names(
      function_names = "nonsense"
    ),
    error_message
  )
})

test_that("cat2", {
  testthat::expect_output(
    cat2(
      message = "test",
      verbose = TRUE
    )
  )
})

test_that("print_info", {
  brts <- c(3, 2, 1)
  n_0 <- 2
  cond <- 1
  testthat::expect_output(
    print_info(
      brts = brts,
      n_0 = n_0,
      cond = cond,
      verbose = TRUE
    )
  )
  brts <- list(c(3, 2, 1), c(2.5, 1.5, 0.5))
  n_0s <- c(2, 1)
  cond <- 1
  testthat::expect_output(
    print_info(
      brts = brts,
      n_0 = n_0,
      cond = cond,
      verbose = TRUE
    )
  )
  testthat::expect_silent(
    print_info(
      brts = brts,
      n_0 = n_0,
      cond = cond,
      verbose = FALSE
    )
  )
})
