import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# First bulletpoint
# Load the data
df = pd.read_csv("C:/Users/40731/Desktop/Assignment4/transformer_data.csv")

# Extract variables
Y = df['Y'].values.reshape(-1, 1)
U = df[['Ta', 'S', 'I']].values
T = len(Y)

# Replace with your estimated values from R (params_2d)
A = np.array([[0.95, 0.05],
              [0.02, 0.90]])  # <-- Estimated A matrix (2x2)
B = np.array([[0.01, 0.02, 0.01],
              [0.015, 0.01, 0.02]])  # <-- Estimated B matrix (2x3)
Qlt = np.array([[1.0, 0.0],
                [0.01, 1.0]])  # <-- Lower triangular for Sigma1
Sigma1 = Qlt @ Qlt.T
C = np.array([[1.0, 0.0]])  # Only observe the first state
Sigma2 = np.array([[1.0]])  # Variance of measurement noise
X0 = np.array([[20.0],
               [20.0]])  # Initial state estimate
P0 = np.eye(2) * 10

# Kalman Filter
X_filtered = np.zeros((T, 2))
x_est = X0
P_est = P0

for t in range(T):
    # Prediction
    x_pred = A @ x_est + B @ U[t].reshape(-1, 1)
    P_pred = A @ P_est @ A.T + Sigma1

    # Update
    y_pred = C @ x_pred
    S = C @ P_pred @ C.T + Sigma2
    K = P_pred @ C.T @ np.linalg.inv(S)
    innov = Y[t] - y_pred
    x_est = x_pred + K @ innov
    P_est = (np.eye(2) - K @ C) @ P_pred

    # Store filtered states
    X_filtered[t] = x_est.ravel()

# Plot the two states
plt.figure(figsize=(10, 5))
plt.plot(X_filtered[:, 0], label='State 1', color='blue')
plt.plot(X_filtered[:, 1], label='State 2', color='pink')
plt.title('Reconstructed latent states over time')
plt.xlabel('Time (hours)')
plt.ylabel('State value')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
 # Add input variables for comparison
Ta = df['Ta'].values
S = df['S'].values
I = df['I'].values


# Second bulletpoint
# Create subplots
plt.figure(figsize=(12, 10))

# Plot State 1
plt.subplot(5, 1, 1)
plt.plot(X_filtered[:, 0], color='blue')
plt.title('State 1')
plt.ylabel('State 1')
plt.grid(True)

# Plot State 2
plt.subplot(5, 1, 2)
plt.plot(X_filtered[:, 1], color='orange')
plt.title('State 2')
plt.ylabel('State 2')
plt.grid(True)

# Plot Outdoor Air Temperature (Ta)
plt.subplot(5, 1, 3)
plt.plot(Ta, color='pink')
plt.title('Outdoor air temperature ($T_{a,t}$)')
plt.ylabel('°C')
plt.grid(True)

# Plot Solar Radiation (S)
plt.subplot(5, 1, 4)
plt.plot(S, color='red')
plt.title('Solar radiation ($\\Phi_{s,t}$)')
plt.ylabel('W/m²')
plt.grid(True)

# Plot Transformer Load (I)
plt.subplot(5, 1, 5)
plt.plot(I, color='purple')
plt.title('Transformer load ($\\Phi_{I,t}$)')
plt.xlabel('Time (hours)')
plt.ylabel('kA')
plt.grid(True)

# Adjust layout
plt.tight_layout()
plt.show()

