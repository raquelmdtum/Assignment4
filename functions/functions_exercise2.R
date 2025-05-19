kf_logLik_dt <- function(par, df, return_residuals = FALSE) {
  # par: vector of parameters
  # df: data frame with observations and inputs as columns (Y, Ta, S, I)
  # par: Could be on the form c(A11, A12, A21, A22, B11, B12, B21, B22, Q11, Q12, Q22)
  A <- matrix(par[1], 1, 1) # transition matrix
  B <- matrix(par[2:4], 1, 3) # input matrix
  Sigma1lt <- matrix(par[5], 1, 1) # lower-triangle of system covariance matrix
  Sigma1   <- Sigma1lt %*% t(Sigma1lt) # THAT IS!!! The system covariance matrix is given by Qlt %*% t(Qlt) (and is this symmetric positive definite)
  C <- matrix(1, 1, 1) # observation matrix
  Sigma2 <- matrix(par[6]^2, 1, 1) # observation noise covariance matrix
  X0 <- matrix(par[7], 1, 1) # initial state

  # Variables
  obs_cols <- c("Y") # observation column names
  input_cols <- c("Ta","S","I") # input column names

  # pull out data
  Y  <- as.matrix(df[, obs_cols])     # m×T
  U  <- as.matrix(df[, input_cols])   # p×T
  Tn <- nrow(df)

  # init
  n      <- nrow(A)
  x_est  <- matrix(Y[1,], n, 1)            # start state from first obs
  x_est <- X0 
  P_est  <- diag(1e1, n)                   # X0 prior covariance
  logLik <- 0

  residuals <- numeric(Tn)

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% t(U[t, , drop = FALSE])
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    y_pred  <- C %*% x_pred
    S_t     <- C %*% P_pred %*% t(C) + Sigma2
    innov   <- Y[t] - y_pred

    residuals[t] <- innov  # store residual

    logLik <- logLik - 0.5 * (log(2 * pi) + log(det(S_t)) + t(innov) %*% solve(S_t, innov))

    K_t    <- P_pred %*% t(C) %*% solve(S_t)
    x_est  <- x_pred + K_t %*% innov
    P_est  <- (diag(n) - K_t %*% C) %*% P_pred
  }

  if (return_residuals) {
    return(residuals)
  } else {
    return(as.numeric(logLik))
  }

  as.numeric(logLik)
}

# Optimizer wrapper
estimate_dt <- function(start_par, df, lower=NULL, upper=NULL) {
  negLL <- function(par){ -kf_logLik_dt(par, df) }
  optim(
    par    = start_par, fn = negLL,
    method = "L-BFGS-B",
    lower  = lower, upper = upper,
    control= list(maxit=1000, trace=1)
  )
}
