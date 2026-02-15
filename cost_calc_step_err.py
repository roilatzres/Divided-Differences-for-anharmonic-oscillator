import numpy as np
import matplotlib.pyplot as plt
from scipy.special import loggamma

def get_optimized_cost_random_walk_error(n, d, total_time=200, eps_base=1e-6):
    """
    Finds the minimum cost assuming Random Walk Error accumulation.
    Constraint: eps_step = eps_base / sqrt(T)
    """
    
    # We scan candidate q values
    # q=1 can cause division by zero in (q-0.5) if we aren't careful, 
    # but mathematically we usually start higher for convergence.
    q_candidates = np.arange(2, 20)
    
    costs = []
    
    # Pre-calculate the numerator term K = d * sqrt(n) * TotalTime
    K = d * np.sqrt(n) * total_time
    
    for q in q_candidates:
        # 1. Calculate Minimum Steps (T)
        # Formula: T = [ K^q / (q! * eps_base) ] ^ (1 / (q - 0.5))
        
        # Using Log-Space for stability
        log_factorial_q = loggamma(q + 1)
        
        # LHS_log = q*ln(K) - ln(q!) - ln(eps_base)
        log_LHS = (q * np.log(K)) - log_factorial_q - np.log(eps_base)
        
        # Power is 1 / (q - 0.5)
        log_T = log_LHS / (q - 0.5)
        
        # T must be at least 1
        T_steps = np.max([1.0, np.exp(log_T)])
        
        # 2. Calculate Cost
        # C = d^q * 2^(q+1) * q^2 * T * n^2
        log_cost = (q * np.log(d)) + \
                   ((q + 1) * np.log(2)) + \
                   (2 * np.log(q)) + \
                   np.log(T_steps) + \
                   (2 * np.log(n))
                   
        costs.append(np.exp(log_cost))
        
    return np.min(costs)

# --- Simulation Parameters ---
N_RANGE = np.arange(1, 3000, 1)  # Photon numbers 1 to 30
EPS_BASE = 1e-4                # Base error 10^-4
TOTAL_TIME = 200

HAMILTONIANS = {
    2: "Dispersive (d=2)",
    4: "Jaynes-Cummings (d=4)",
    6: "Rabi (d=6)",
    10: "more terms (d=10)",
    14: "even more terms (d=12)",
    20: "extreme case (d=20)",
    25: "very extreme (d=25)",
    30: "ultra extreme (d=30)"
}

# --- Plotting ---
plt.figure(figsize=(10, 7))

# Baselines for context
ONE_HOUR_OPS = 1.3e12 * 3600    # RTX 4090 FP64
ONE_DAY_OPS = ONE_HOUR_OPS * 24

for d, label in HAMILTONIANS.items():
    optimized_costs = []
    print(f"Calculating costs for Hamiltonian with branching factor d={d}...")
    for n in N_RANGE:
        min_cost = get_optimized_cost_random_walk_error(n, d, TOTAL_TIME, EPS_BASE)
        optimized_costs.append(min_cost)
    
    plt.plot(N_RANGE, optimized_costs, 'o--', linewidth=2, label=label)

# # --- Visual Baselines ---
# plt.axhline(y=ONE_HOUR_OPS, color='gray', linestyle=':', alpha=0.5)
# plt.text(1, ONE_HOUR_OPS*1.5, "1 Hour (RTX 4090 FP64)", color='gray', fontsize=9)

# plt.axhline(y=ONE_DAY_OPS, color='gray', linestyle=':', alpha=0.5)
# plt.text(1, ONE_DAY_OPS*1.5, "1 Day (RTX 4090 FP64)", color='gray', fontsize=9)

# plt.axhline(y=1e20, color='red', linestyle='--', alpha=0.3)
# plt.text(1, 1e20*1.5, "Ph.D. Limit (100 ExaFLOPs)", color='red', fontsize=9)

# --- Formatting ---
plt.yscale('log')
plt.xlabel(r'Max Photon Number ($n_{max}$)', fontsize=14)
plt.ylabel('Optimized Complexity Cost (FLOPs)', fontsize=14)
plt.title(f'Optimized Complexity with Diffusive Error\n$\epsilon_{{step}} = 10^{{-4}} / \sqrt{{T}}$', fontsize=16)
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.legend(fontsize=12)

# Info Box
info_text = (
    r"$\bf{Constraint:}$" + "\n"
    r"Error accumulates diffusively" + "\n"
    r"$\epsilon_{total} \propto \sqrt{T} \cdot \epsilon_{step}$" + "\n"
    r"Allows for slightly looser" + "\n"
    r"tolerances at large $T$."
)
plt.text(0.95, 0.05, info_text, transform=plt.gca().transAxes, 
         fontsize=10, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", alpha=0.8))

plt.tight_layout()
plt.show()