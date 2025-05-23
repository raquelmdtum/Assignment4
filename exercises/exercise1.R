#-------1.1-----
# Set a random seed to ensure reproducibility
set.seed(42)

# Parameter setup
a <- 0.9           # State transition coefficient
b <- 1             # Bias term
sigma <- 1         # Standard deviation of noise
X0 <- 5            # Initial value of the process
n <- 100           # Number of time steps
num_trajectories <- 5  # Number of independent simulations

# Create a matrix to store all trajectories (rows: time steps, columns: trajectories)
trajectories <- matrix(NA, nrow = n + 1, ncol = num_trajectories)

# Simulate each trajectory
for (j in 1:num_trajectories) {
  X <- numeric(n + 1)  # Create a vector to hold the state values
  X[1] <- X0           # Set initial value
  for (t in 2:(n + 1)) {
    e <- rnorm(1, mean = 0, sd = sigma)     # Generate random noise
    X[t] <- a * X[t - 1] + b + e            # Apply the state update equation
  }
  trajectories[, j] <- X  # Store this trajectory in the matrix
}

# Plotting
time <- 0:n                        # Time vector
colors <- rainbow(num_trajectories)  # Generate different colors for each trajectory

# Plot the first trajectory to set up the figure
windows(width = 12, height = 6)   
plot(time, trajectories[,1], type = "l", col = colors[1], ylim = range(trajectories),
     xlab = "Time", ylab = expression(X[t]), main = "5 Trajectories of X_t")

# Add the remaining trajectories
for (j in 2:num_trajectories) {
  lines(time, trajectories[,j], col = colors[j])
}

# Add a legend to label each trajectory
legend("topleft", legend = paste("Trajectory", 1:num_trajectories),
       col = colors, lty = 1, cex = 0.8)


#-------1.2---------
# Set random seed for reproducibility
set.seed(42)

# Parameters
a <- 0.9
b <- 1
sigma1 <- 1  # Standard deviation for process noise
sigma2 <- 1  # Standard deviation for observation noise
X0 <- 5
n <- 100

# Initialize vectors
X <- numeric(n + 1)  # Latent state
Y <- numeric(n + 1)  # Observations
X[1] <- X0           # Initial state

# Simulate process
for (t in 2:(n + 1)) {
  e1 <- rnorm(1, mean = 0, sd = sigma1)
  X[t] <- a * X[t - 1] + b + e1
}

# Generate observations
e2 <- rnorm(n + 1, mean = 0, sd = sigma2)
Y <- X + e2

# Plot
time <- 0:n
windows(width = 12, height = 6)   
plot(time, X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y)),
     xlab = "Time", ylab = "Value", main = "True State vs Noisy Observations")
lines(time, Y, col = "red", lty = 2, lwd = 2)

legend("topleft", legend = c("True state (X_t)", "Observation (Y_t)"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)



#--------1.3-----

# ------------set the parameters ------------
set.seed(42)
a <- 0.9
b <- 1
sigma1 <- 1
sigma2 <- 1
X0 <- 5
n <- 100
time <- 0:n



# ------------ Kalman Filter definition ------------
myKalmanFilter <- function(y, theta, R, x_prior = 0, P_prior = 10) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  N <- length(y)
  
  x_pred  <- numeric(N)
  P_pred  <- numeric(N)
  x_filt  <- numeric(N)
  P_filt  <- numeric(N)
  innovation     <- numeric(N)
  innovation_var <- numeric(N)
  
  # initial prediction
  x_pred[1] <- x_prior
  P_pred[1] <- P_prior
  innovation[1] <- y[1] - x_pred[1]
  innovation_var[1] <- P_pred[1] + R
  K_1 <- P_pred[1] / innovation_var[1]
  x_filt[1] <- x_pred[1] + K_1 * innovation[1]
  P_filt[1] <- (1 - K_1) * P_pred[1]
  
  # main loop
  for (t in 2:N) {
    x_pred[t] <- a * x_filt[t - 1] + b
    P_pred[t] <- a^2 * P_filt[t - 1] + sigma1^2
    
    innovation[t] <- y[t] - x_pred[t]
    innovation_var[t] <- P_pred[t] + R
    K_t <- P_pred[t] / innovation_var[t]
    x_filt[t] <- x_pred[t] + K_t * innovation[t]
    P_filt[t] <- (1 - K_t) * P_pred[t]
  }
  
  return(list(
    x_pred = x_pred,
    P_pred = P_pred,
    x_filt = x_filt,
    P_filt = P_filt,
    innovation = innovation,
    innovation_var = innovation_var
  ))
}

# ------------ apply Kalman Filter  ------------
theta <- c(a, b, sigma1)
R <- sigma2^2
kf <- myKalmanFilter(Y, theta, R, x_prior = 0, P_prior = 10)

# confidence interval
lower <- kf$x_pred - 1.96 * sqrt(kf$P_pred)
upper <- kf$x_pred + 1.96 * sqrt(kf$P_pred)

# set the size
y_min <- min(c(X, Y, kf$x_pred, kf$x_filt))
y_max <- max(c(X, Y, kf$x_pred, kf$x_filt))
par(mar = c(5, 4, 4, 8), xpd = TRUE)  # 增加右边距，并允许越界绘图

# open a new window
windows(width = 12, height = 10)

# set the size
par(mar = c(5, 4, 4, 8), xpd = TRUE)

# the limir of CI
ci_min <- min(lower)
ci_max <- max(upper)

# add the margin
y_min <- min(c(X, Y, kf$x_pred, kf$x_filt))
y_max <- max(c(X, Y, kf$x_pred, kf$x_filt))

# correct the range of axis y
buffer <- 0.5  # buffer
ylim_low <- min(y_min, ci_min) - buffer
ylim_high <- max(y_max, ci_max) + buffer

# Open a new plotting window and set margins
windows(width = 12, height = 10)
par(mar = c(5, 4, 4, 2) + 0.1)

# Plotting
plot(time, X, type = "l", col = "black", lwd = 2,
     ylim = c(ylim_low, ylim_high), xlim = c(min(time), max(time)),
     xlab = "Time", ylab = "Value", main = "Kalman Filter Prediction")

# Add confidence interval shading
polygon(c(time, rev(time)), c(lower, rev(upper)), 
        col = rgb(0, 0, 1, 0.1), border = NA)

# Add additional lines
lines(time, Y, col = "red", lty = 2, lwd = 1.5)
lines(time, kf$x_pred, col = "blue", lwd = 2)
lines(time, kf$x_filt, col = "darkgreen", lty = 4, lwd = 2)

# Legend in the bottom right corner
legend("bottomright",
       legend = c("True State (X)", "Observation (Y)", "Predicted (X̂_t+1|t)",
                  "Filtered (X̂_t|t)", "95% CI"),
       col = c("black", "red", "blue", "darkgreen", rgb(0, 0, 1, 0.1)),
       lty = c(1, 2, 1, 4, NA), lwd = c(2, 1.5, 2, 2, NA),
       fill = c(NA, NA, NA, NA, rgb(0, 0, 1, 0.1)),
       border = NA, bty = "n")


#-------1.4------for (1,0.9,1),(1,0.9,5)

myLogLikFun <- function(theta, y, R, x_prior = 0, P_prior = 10) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  
  # apply Kalman filter
  kf_result <- myKalmanFilter(y, theta, R, x_prior, P_prior)
  
  err <- kf_result$innovation       # predict the error
  S <- kf_result$innovation_var     # variance of error
  
  # delete outliers
  if (any(S <= 0)) return(Inf)
  
  # calculate the log-likelihood value
  logL <- -0.5 * sum(log(2 * pi * S) + (err^2) / S)
  
  return(-logL)  
}

# ------------ Simulation and Estimation ------------
set.seed(42)
n <- 100
R <- 1   # sigma2^2
true_theta <- c(1, 0.9, 5)  # a, b, sigma1
N_sim <- 100

results <- matrix(NA, nrow = N_sim, ncol = 3)
colnames(results) <- c("a", "b", "sigma1")

for (i in 1:N_sim) {
  # simulate latent X and observations Y
  X <- numeric(n + 1)
  X[1] <- 5
  for (t in 2:(n + 1)) {
    X[t] <- true_theta[1] * X[t - 1] + true_theta[2] + rnorm(1, 0, true_theta[3])
  }
  Y <- X + rnorm(n + 1, 0, sqrt(R))
  
  # initial guess for parameters
  theta_init <- c(0.5, 0.5, 0.5)
  fit <- optim(theta_init, myLogLikFun, y = Y, R = R, method = "L-BFGS-B",
               lower = c(-Inf, -Inf, 1e-4))
  
  results[i, ] <- fit$par
}

# ------------ Visualization ------------
windows(width = 10, height = 6)  # open a new window
boxplot(results, names = c("a", "b", "sigma1"), main = "Estimated Parameters from 100 Simulations")
points(x = 1:3, y = true_theta, col = "red", pch = 18, cex = 1.5)
legend("topright", legend = "True values", col = "red", pch = 18)


##-------1.4 for (5,0.9,1)------

myLogLikFun <- function(theta, y, R, x_prior = 0, P_prior = 10) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  
  # Run the Kalman filter
  kf_result <- myKalmanFilter(y, theta, R, x_prior, P_prior)
  
  err <- kf_result$innovation       # Prediction errors (innovations)
  S <- kf_result$innovation_var     # Innovation variances
  
  # Return Inf if any variance is non-positive (to avoid invalid likelihood)
  if (any(S <= 0)) return(Inf)
  
  # Compute the total log-likelihood (Gaussian log-density for each time step)
  logL <- -0.5 * sum(log(2 * pi * S) + (err^2) / S)
  
  return(-logL)  # Return negative log-likelihood for minimization
}


# ------------ Simulation and Estimation ------------
set.seed(42)
n <- 100
R <- 1   # sigma2^2
true_theta <- c(5, 0.9, 1)  # a, b, sigma1
N_sim <- 100

results <- matrix(NA, nrow = N_sim, ncol = 3)
colnames(results) <- c("a", "b", "sigma1")

for (i in 1:N_sim) {
  # simulate latent X and observations Y
  X <- numeric(n + 1)
  X[1] <- 5
  for (t in 2:(n + 1)) {
    X[t] <- true_theta[1] * X[t - 1] + true_theta[2] + rnorm(1, 0, true_theta[3])
  }
  Y <- X + rnorm(n + 1, 0, sqrt(R))
  
  # initial guess for parameters
  theta_init <- c(0.5, 0.5, 0.5)
  fit <- optim(theta_init, myLogLikFun, y = Y, R = R,
               method = "L-BFGS-B", lower = c(-2, -Inf, 1e-3), upper = c(2, Inf, 10))
  
  
  results[i, ] <- fit$par
}

# ------------ Visualization ------------
windows(width = 10, height = 6)  # 打开新窗口
boxplot(results, names = c("a", "b", "sigma1"), main = "Estimated Parameters from 100 Simulations")
points(x = 1:3, y = true_theta, col = "red", pch = 18, cex = 1.5)
legend("topright", legend = "True values", col = "red", pch = 18)

#----1.5--
# ----- Parameter settings -----
set.seed(42)
n <- 100                   # Length of the time series
a <- 1                    # State transition coefficient
b <- 0.9                  # Constant term in state update
sigma1 <- 1              # Standard deviation of process noise
sigma2 <- 1              # Standard deviation of observation noise
nu_values <- c(100, 5, 2, 1)  # Degrees of freedom for the t-distribution

# ---- Simulate one realization of latent state X and observation Y (using t-distributed process noise) ----
lambda <- rt(n + 1, df = 2)   # Example with df = 2 (can change to other values)
X <- numeric(n + 1)
X[1] <- 0
for (t in 2:(n + 1)) {
  X[t] <- a * X[t - 1] + b + sigma1 * lambda[t]   # State update equation with t noise
}
Y <- X + rnorm(n + 1, 0, sigma2)  # Add Gaussian observation noise to generate Y

# ---- Plot density: t-distribution vs standard normal ----
x_vals <- seq(-10, 10, length.out = 1000)  # X-axis for density plot
windows(width = 10, height = 6)            # Open new plotting window

# Plot standard normal density
plot(x_vals, dnorm(x_vals), type = "l", lwd = 2, col = "black",
     ylab = "Density", main = "t-distributions vs Standard Normal")

# Set line colors and types
cols <- c("blue", "red", "green", "purple")
ltys <- c(2, 3, 4, 5)

# Overlay t-distribution densities with various degrees of freedom
for (i in seq_along(nu_values)) {
  lines(x_vals, dt(x_vals, df = nu_values[i]), col = cols[i], lty = ltys[i], lwd = 2)
}

# Add legend
legend("topright", legend = c("Normal",
                              paste0("t(df=", nu_values, ")")),
       col = c("black", cols), lty = c(1, ltys), lwd = 2)

#-----1.5----second part
set.seed(42)

n <- 100                             # Length of the time series
R <- 1                               # Observation noise variance
true_theta <- c(1, 0.9, 1)           # True parameters: a, b, sigma1
nu_values <- c(100, 5, 2, 1)         # Degrees of freedom for the t-distribution
N_sim <- 100                         # Number of Monte Carlo simulations

results_list <- list()              # To store parameter estimates for each nu

for (nu in nu_values) {
  sim_results <- matrix(NA, nrow = N_sim, ncol = 3)  # Each row stores (a, b, sigma1)
  
  for (i in 1:N_sim) {
    # Generate system noise using t-distribution with current nu
    lambda <- rt(n + 1, df = nu)
    X <- numeric(n + 1)
    X[1] <- 0
    for (t in 2:(n + 1)) {
      X[t] <- true_theta[1] * X[t - 1] + true_theta[2] + true_theta[3] * lambda[t]
    }
    
    # Generate noisy observations
    Y <- X + rnorm(n + 1, mean = 0, sd = sqrt(R))
    
    # Robust initial guess for optimization
    theta_init <- c(runif(1, 0.5, 1.5), runif(1, 0.5, 1.5), runif(1, 0.5, 1.5))
    
    # Perform maximum likelihood estimation (under Gaussian assumption)
    fit <- optim(theta_init, myLogLikFun, y = Y, R = R,
                 method = "L-BFGS-B", lower = c(-Inf, -Inf, 1e-4))
    
    # Check convergence and print diagnostic info if needed
    if (fit$convergence != 0 || any(is.na(fit$par))) {
      cat("⚠️ Optim failed: iter =", i, "nu =", nu, "\n")
      cat("  theta_init =", paste(round(theta_init, 2), collapse = ", "), "\n")
      cat("  fit$par    =", paste(round(fit$par, 4), collapse = ", "), "\n")
      cat("  logLik     =", fit$value, "\n\n")
    }
    
    sim_results[i, ] <- fit$par  # Store estimated parameters
  }
  
  results_list[[paste0("nu_", nu)]] <- sim_results  # Save results for this nu
}

#------------- Visualization of parameter estimates for each nu ---------------
true_theta <- c(1, 0.9, 5)                      # True values for plotting reference
param_names <- c("a", "b", "sigma1")            # Axis labels
nu_labels <- paste0("nu = ", nu_values)         # Plot titles

windows(width = 12, height = 6)                 # Open new plotting window
par(mfrow = c(1, length(nu_values)), mar = c(5, 4, 4, 2))  # Layout settings

for (k in seq_along(nu_values)) {
  nu <- nu_values[k]
  sim_results <- results_list[[paste0("nu_", nu)]]
  
  # Draw boxplot for the estimated parameters
  boxplot(sim_results, names = param_names,
          main = nu_labels[k], ylim = c(0, 2), col = "white")
  
  # Add red points for true parameter values
  points(x = 1:3, y = true_theta, col = "red", pch = 18, cex = 1.5)
}
