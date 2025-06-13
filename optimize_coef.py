import json
import numpy as np
from scipy.optimize import minimize

# === Load coefficient matrices from JSON ===
with open("all_coefs.json", "r") as f:
    data = json.load(f)

num_times = len(data)
num_q = len(data[0])
num_start = len(data[0][0])
num_target = len(data[0][0][0])

# Build complex coefficient tensor: shape [T][d+1][N][N]
degree = num_q  # number of powers (e.g., c0 to cd)
C = np.zeros((num_times, num_q, num_start, num_target), dtype=np.complex128)

for t in range(num_times):
    for q in range(num_q):
        for i in range(num_start):
            for j in range(num_target):
                real = data[t][q][i][j]['real']
                imag = data[t][q][i][j]['imag']
                C[t, q, i, j] = real + 1j * imag

print(f"Loaded coefficient tensor C with shape: {C.shape}")
# === Prepare target vector and start vector ===

# Load target data
with open("amplitudes_e_chg_pulse3amp0.002000_q10.json", "r") as f:
    target_data = json.load(f)

num_time = len(target_data)
num_target = len(target_data[0]) if num_time > 0 else 0


target_vector = np.zeros(num_target, dtype=np.complex128)
for j in range(num_target):
    real = target_data[-1][j]['real']  # take final time index
    imag = target_data[-1][j]['imag']
    target_vector[j] = real + 1j * imag

# Print target vector information
print(f"Target vector shape: {target_vector.shape}")
print(f"Target vector: {target_vector}")

# Normalize the target vector
target_vector /= np.linalg.norm(target_vector)

# Print normalized target vector
print(f"Normalized target vector: {target_vector}")

# Create the start vector
start_vector = np.zeros(num_start, dtype=np.complex128)
start_vector[0] = 1.0 + 0.0j  # (1, 0, 0, ..., 0)

# Print start vector information
print(f"Start vector shape: {start_vector.shape}")
print(f"Start vector: {start_vector}")

# === Polynomial evaluation ===
def P(x, C_t):
    """Evaluate matrix polynomial at scalar x using coefficients C_t[q]"""
    result = np.zeros_like(C_t[0], dtype=np.complex128)
    for q in range(len(C_t)):
        result += C_t[q] * x**q
    return result

# === Objective function: 1 - |<target | M(x) | start>|^2 ===
def objective(x_vals):
    M = np.eye(num_target, dtype=np.complex128)
    for x in x_vals:
        M = M @ P(x, C[-1])
    psi = M @ start_vector
    psi /= np.linalg.norm(psi)  # normalize intermediate state
    fidelity = np.abs(np.vdot(target_vector, psi))**2
    return 1.0 - fidelity

# === Run optimization ===
x0 = np.zeros(4)  # or any number of steps

result = minimize(objective, x0, method='BFGS', options={'disp': True})

# === Print results ===
print("\nOptimal x:", result.x)
print("Minimum cost:", result.fun)
print("Final fidelity:", 1 - result.fun)
print("num iterations:", result.nit)
print("Success:", result.success)

# print final state vector
final_M = np.eye(num_target, dtype=np.complex128)
for x in result.x:
    final_M = final_M @ P(x, C[-1])
final_state = final_M @ start_vector
final_state /= np.linalg.norm(final_state)  # normalize final state
print("Final state vector:", final_state)
print("Final state vector norm:", np.linalg.norm(final_state))