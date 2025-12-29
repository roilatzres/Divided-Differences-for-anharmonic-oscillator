import json
import subprocess
import os
from itertools import product

# Path to the compiled C++ executable
EXECUTABLE = "./main_jc_multi_piece_base"  # adjust if needed

# Base parameter template
base_params = {
    "final_t": 192,
    "multiplier": 4,
    "sigma": 48,
    "timestep": 1,
    "chi": -279e-6,
    "alpha": 200e-3,
    "delta": 2,
    "amplitude": 0,
    "q_max": 4,
    "max_target_cavity": 4,
    "max_target_qubit": 2,
    "start_state_cavity": 0,
    "start_state_qubit": 1,
    "num_pulses": 192,
    "out_dir": "jc_base/",
    "max_total_cavity": 300,
    "max_total_qubit": 2
}

# Define parameter ranges to sweep over
amplitudes = [1e-7]
q_max_values = [4,7]
num_pulses_list = [192]
max_target_cavity = [3,6,9]
max_target_qubit = [2]
sq_list = [0,1]
# # Define parameter ranges to sweep over
# amplitudes = [0.0001, 0.0005]
# q_max_values = [5, 6, 7, 8]
# num_pulses_list = [192, 384, 640 ,960]
# max_target_cavity = [4,5,6,7,8]
# max_target_qubit = [2]

# Create output directory if needed
if not os.path.exists(base_params["out_dir"]):
    os.makedirs(base_params["out_dir"])

# Sweep over all combinations
for amp, q_max, num_pulses, n_cav, n_qub,sq in product(amplitudes, q_max_values, num_pulses_list, max_target_cavity, max_target_qubit, sq_list):
    print(f"\nRunning with amplitude={amp}, q_max={q_max}, pulses={num_pulses}, cavity={n_cav}, qubit={n_qub}, sq={sq}")
    # Update parameters
    params = base_params.copy()
    params["amplitude"] = amp
    params["q_max"] = q_max
    params["num_pulses"] = num_pulses
    params["max_target_cavity"] = n_cav
    params["max_target_qubit"] = n_qub
    params["start_state_qubit"] = sq
    # Optionally set dynamic output directory
    params["out_dir"] = f"jc_base/a{amp}_q{q_max}_p{num_pulses}_sq{sq}/"
    os.makedirs(params["out_dir"], exist_ok=True)

    # Write to params.json
    with open("params.json", "w") as f:
        json.dump(params, f, indent=4)

    # Run the C++ executable
    result = subprocess.run([EXECUTABLE], capture_output=True, text=True)

    # Log output
    log_file = os.path.join(params["out_dir"], f"run_a{amp}_q{q_max}_nc{n_cav}_nq{n_qub}_p{num_pulses}_sq{sq}.log")
    with open(log_file, "w") as logf:
        logf.write(result.stdout)
        logf.write(result.stderr)

    print("Simulation finished. Output saved.")
