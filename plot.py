import numpy as np
import matplotlib.pyplot as plt

# Define parameters
k_0 = 1.0 # Example value for k_0 
alpha = 0.5 # Example value for alpha 
c = 3*10**8 # Example value for c

# Define the function C(1) 
def C(l, k_0, alpha, c):
    prefactor = -np.pi / c**2 
    cosine_term = (k_0 / alpha) * np.cos(k_0 * l)
    sine_term= np.sin(k_0 * l)
    exponential = np.exp(-alpha * l)
    return prefactor * (cosine_term - sine_term) * exponential

# Generate l values
l = np.linspace(0, 15, 1000) # 1 ranges from 0 to 10

# Compute C(l)
C_values = C(l, k_0, alpha, c)

# Plot the function
plt.figure(figsize=(9, 6))
plt.plot(l, C_values, label=r"$C(l)$", color='r') 
# Extend y-axis limits
y_min, y_max = C_values.min(), C_values.max()
y_range = y_max - y_min
plt.ylim([y_min - 0.1 * y_range, y_max + 0.5 * y_range])  # Extend by 10% on both sides

# plt.title("$C(L)$", fontsize=14)
plt.xlabel("$L$", fontsize=12)
plt.ylabel("$C(L)$", fontsize=12)
# plt.grid(True)
plt.legend (fontsize=12)
plt.show()