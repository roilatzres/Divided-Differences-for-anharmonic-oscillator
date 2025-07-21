import json
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

# === 1. load coefficient‐matrices from JSON ===
def load_C(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    T   = len(data)
    d1  = len(data[0])
    Nc  = len(data[0][0])
    assert Nc == len(data[0][0][0]), "cavity matrices must be square"
    C = np.zeros((T, d1, Nc, Nc), dtype=np.complex128)
    for t in range(T):
        for q in range(d1):
            for i in range(Nc):
                for j in range(Nc):
                    r = data[t][q][i][j]["real"]
                    im= data[t][q][i][j]["imag"]
                    C[t, q, i, j] = r + 1j*im
    return C

# load your two gate‐sets (replace with your actual file paths)
Ce = load_C("all_coefs_e.json")    # coefficients for Ue
Cg = load_C("all_coefs_g.json")    # coefficients for Ug

# === 2. dimensions ===
T, d1, dim_cavity, _ = Cg.shape
n_steps     = 10
dim_qubit   = 2
dim_total   = dim_qubit * dim_cavity

# === 3. qubit projectors & paulis ===
P_g = np.array([[1, 0], [0, 0]], dtype=complex)
P_e = np.array([[0, 0], [0, 1]], dtype=complex)
sigma_x  = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y  = np.array([[0, -1j], [1j, 0]], dtype=complex)

# === 4. cavity‐only polynomial gate ===
def P_poly(x, C_t):
    M = np.zeros((dim_cavity, dim_cavity), dtype=complex)
    for q in range(len(C_t)):
        M += C_t[q] * x**q
    return M

def controlled_U(x, Cg_t, Ce_t):
    Ug = P_poly(x, Cg_t)
    Ue = P_poly(x, Ce_t)
    return np.kron(P_g, Ug) + np.kron(P_e, Ue)

# === 5. arbitrary‐axis qubit rotation ===
def R_axis(phi, theta):
    axis_op = sigma_x*np.cos(phi) + sigma_y*np.sin(phi)
    return expm(-1j * theta/2 * axis_op)

def rotation_on_total(phi, theta):
    Rq = R_axis(phi, theta)
    return np.kron(Rq, np.eye(dim_cavity, dtype=complex))

# === 6. one‐step propagator ===
def U_step(x, phi, theta, Cg_t, Ce_t):
    return controlled_U(x, Cg_t, Ce_t) @ rotation_on_total(phi, theta)

# === 7. build start & target vectors in joint space ===
start = np.zeros(dim_total, dtype=complex)
start[0] = 1.0  # |g>⊗|0>

# target = np.zeros(dim_total, dtype=complex)
# m = 2  # desired cavity Fock level
# target[dim_cavity + m] = 1.0  # |e>⊗|m>
# target /= np.linalg.norm(target)
# 1) build the cavity Fock projector |n><n|

n = 1  # the Fock level you want
fock_n = np.zeros((dim_cavity,), dtype=complex)
fock_n[n] = 1.0
P_cav = np.outer(fock_n, fock_n.conj())    # shape (dim_cavity,dim_cavity)

# 2) lift to total space
P_total = np.kron(np.eye(dim_qubit, dtype=complex), P_cav)
# now P_total has shape (dim_total, dim_total)


# === 8. objective over (x_i, φ_i, θ_i) for each step ===
num_vars = n_steps * 3  # each step has (x, φ, θ)

def objective(vars):
    M = np.eye(dim_total, dtype=complex)
    for i in range(n_steps):
        x_i     = vars[3*i]
        phi_i   = vars[3*i + 1]
        theta_i = vars[3*i + 2]
        M = M @ U_step(x_i, phi_i, theta_i, Cg[-1], Ce[-1])
    psi = M @ start
    psi /= np.linalg.norm(psi)

    F = np.real(np.vdot(psi, P_total @ psi))
    return 1 - F  # fidelity squared, we want to minimize 1 - fidelity TODO: maybe need F**2?
    # return 1 - np.abs(np.vdot(target, psi))**2

# === 9. initial guess & bounds ===
x0 = np.zeros(num_vars, dtype=float)
for i in range(n_steps):
    x0[3*i + 0] = 0.0        # x_i
    x0[3*i + 1] = 0.1        # φ_i
    x0[3*i + 2] = np.pi/4    # θ_i

bounds = []
for _ in range(n_steps):
    bounds += [(None, None),    # x_i
               (0.0, 2*np.pi),  # φ_i
               (0.0, 2*np.pi)]  # θ_i

# === 10. run optimization ===
res = minimize(
    objective,
    x0,
    method='L-BFGS-B',
    bounds=bounds,
    options={'disp': True}
)

# === 11. unpack & report ===
opt = res.x
opt_x     = opt[0::3]
opt_phi   = opt[1::3]
opt_theta = opt[2::3]
opt_fid   = 1 - res.fun

print("Optimized x values:     ", opt_x)
print("Optimized φ values (rad):", opt_phi)
print("Optimized θ values (rad):", opt_theta)
print("Final fidelity:          ", opt_fid)
print("Success?                 ", res.success)
