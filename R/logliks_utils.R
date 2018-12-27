#' @title Transition matrix builder
#' @author Giovanni Laudanno
#' @description Builds the transition matrix to integrate the differential equations of the P-equation
#' @inheritParams default_params_doc
#' @return The transition matrix
#' @export
p_transition_matrix <- function(
 lambda,
 mu,
 matrix_size
) {
 nvec <- 0:(matrix_size - 1)
 m <- matrix(0, nrow = matrix_size, ncol = matrix_size)
 m[row(m) == col(m) + 1] <- m[row(m) == col(m) + 1] +
  lambda * (nvec[1:(matrix_size - 1)])
 m[row(m) == col(m) - 1] <- m[row(m) == col(m) - 1] +
  mu * (nvec[2:matrix_size])
 m[row(m) == col(m)] <- m[row(m) == col(m)] -
  (lambda + mu) * (nvec[1:matrix_size])
 m[matrix_size, matrix_size] <- - m[matrix_size - 1, matrix_size]; m
 testit::assert(colSums(m) < 1e-10)
 return(m)
}