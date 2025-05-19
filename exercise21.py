# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset using the correct file path
file_path = r"C:\Users\40731\Desktop\Assignment4\transformer_data.csv"
data = pd.read_csv(file_path)

# Display the first few rows to confirm data is loaded correctly
print("First few rows of the dataset:")
print(data.head())

# Print column names for reference
print("\nColumn names in the dataset:")
print(data.columns)

# Extract the columns using the correct names
time = data['time']        # Time in hours
Y = data['Y']              # Transformer station temperature (°C)
Ta = data['Ta']            # Outdoor air temperature (°C)
S = data['S']              # Horizontal global solar radiation (W/m²)
I = data['I']              # Load on the transformer station (kA)

# Set up the plot
plt.figure(figsize=(12, 10))  # Set figure size for better visibility

# Plot 1: Transformer station temperature (Y)
plt.subplot(4, 1, 1)
plt.plot(time, Y, color='pink')
plt.title('Transformer station temperature (Y)')
plt.ylabel('Temperature (°C)')
plt.xlabel('Time (hours)', labelpad=0)
plt.grid(True)

# Plot 2: Outdoor air temperature (Ta)
plt.subplot(4, 1, 2)
plt.plot(time, Ta, color='blue')
plt.title('Outdoor air temperature (Ta)')
plt.ylabel('Temperature (°C)')
plt.xlabel('Time (hours)', labelpad=0)
plt.grid(True)

# Plot 3: Horizontal global solar radiation (S)
plt.subplot(4, 1, 3)
plt.plot(time, S, color='orange')
plt.title('Horizontal global solar radiation (S)')
plt.ylabel('Radiation (W/m²)')
plt.xlabel('Time (hours)', labelpad=0)
plt.grid(True)

# Plot 4: Transformer load (I)
plt.subplot(4, 1, 4)
plt.plot(time, I, color='red')
plt.title('Transformer load (I)')
plt.ylabel('Load (kA)')
plt.xlabel('Time (hours)', labelpad=0)
plt.grid(True)

# Adjust layout to avoid overlapping of plots
plt.tight_layout()
plt.subplots_adjust(hspace=0.6)

# Show all the plots
plt.show()
