source("functions/functions_exercise2.R")
df <- read.csv("transformer_data.csv")

# ---------- Exercise 2.3 ----------

kf_logLik_dt <- function(par, df, return_residuals = FALSE) {
  A <- matrix(par[1:4], 2, 2)
  B <- matrix(par[5:10], 2, 3)
  Qlt <- matrix(0, 2, 2)
  Qlt[1, 1] <- par[11]
  Qlt[2, 1] <- par[12]
  Qlt[2, 2] <- par[13]
  Sigma1 <- Qlt %*% t(Qlt)
  C <- matrix(c(1, 0), nrow = 1)
  Sigma2 <- matrix(par[14]^2, 1, 1)
  X0 <- matrix(par[15:16], 2, 1)

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

start_par <- c(
  0.9, 0, 0, 0.9,     # A (4)
  0.01, 0.01, 0.01,   # B row 1 (3)
  0.01, 0.01, 0.01,   # B row 2 (3)
  1.0, 0.01, 1.0,     # Qlt (3)
  1.0,                # sigma2
  20.0, 20.0          # initial X0
)
lower <- rep(-10, length(start_par))
upper <- rep(10, length(start_par))
lower[11] <- lower[13] <- 1e-3   # Q diagonal elements must be >0
lower[14] <- 1e-3                # sigma2 > 0

result_2d <- estimate_dt(start_par, df, lower, upper)
params_2d <- result_2d$par

residuals_2d <- kf_logLik_dt(params_2d, df, return_residuals = TRUE)

pdf("images/exer23.pdf", width = 10, height = 6)
par(mfrow = c(2, 2))

plot(residuals_2d, type = 'l', main = 'Residuals', ylab = "Residuals", xlab = "Time")
acf(residuals_2d, main = 'ACF of Residuals')
pacf(residuals_2d, main = 'PACF of Residuals')
qqnorm(residuals_2d); qqline(residuals_2d)

logLik_2d <- -kf_logLik_dt(params_2d, df)
n <- nrow(df)
k <- length(params_2d)
AIC_2d <- -2 * logLik_2d + 2 * k
BIC_2d <- -2 * logLik_2d + log(n) * k
cat("AIC (2D):", AIC_2d, "\nBIC (2D):", BIC_2d, "\n")

dev.off()

