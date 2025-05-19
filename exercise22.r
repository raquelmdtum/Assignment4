source("functions/functions_exercise2.R")
df <- read.csv("transformer_data.csv")

# Run parameter estimation
start_par <- c(0.9, 0.01, 0.01, 0.01, 1.0, 1.0, 20.0)
lower <- c(-1.5, -1, -1, -1, 1e-3, 1e-3, 10)
upper <- c(1.5, 1, 1, 1, 10, 10, 30)

result <- estimate_dt(start_par, df, lower, upper)
params <- result$par
print(params)

# Residuals
residuals <- kf_logLik_dt(params, df, return_residuals = TRUE)

# Plots
par(mfrow=c(2,2))
plot(residuals, type='l', main='Residuals')
acf(residuals, main='ACF of Residuals')
pacf(residuals, main='PACF of Residuals')
qqnorm(residuals); qqline(residuals)

# AIC and BIC
logLik <- -kf_logLik_dt(params, df)
n <- nrow(df)
k <- length(params)
AIC <- -2 * logLik + 2 * k
BIC <- -2 * logLik + log(n) * k
cat("AIC:", AIC, "\nBIC:", BIC, "\n")

# Residual diagnostics
pdf("images/exer22.pdf", width = 10, height = 6)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Residuals plot
plot(residuals, type = 'l', main = 'Residuals', ylab = "Residuals", xlab = "Time")

# ACF plot
acf(residuals, main = 'Autocorrelation Function (ACF)')

# PACF plot
pacf(residuals, main = 'Partial Autocorrelation Function (PACF)')

# QQ plot
qqnorm(residuals, main = 'Normal Q-Q Plot')
qqline(residuals)

dev.off()