import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.signal import gaussian, square
from warnings import warn
from qutip import *


def show_amp_json_dispersive(amp_json):
    """
    Show transition amplitudes from JSON files.
    """
    # print(amp_json)
    if not amp_json:
        warn("No amplitude JSON files provided.")
        return

    # Check if the files exist
    for file in amp_json:
        if not Path(file).is_file():
            warn(f"File {file} does not exist.")
            return

    # Load and process each file
    for file in amp_json:
        pulse_num = file.split('pulse')[1].split('amp')[0]

        with open(file, 'r') as f:
            data = json.load(f)

        all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]
        N_plot = len(all_transition_amplitudes[0])
        # print(N_plot)

        psi_t = []
        # create Qobj for each time frame in the transition amplitudes
        for t,cn_t in enumerate(all_transition_amplitudes):
            psi = basis(N_plot, 0)*0
            for n,cn in enumerate(cn_t):
                psi += cn*basis(N_plot, n)
            psi_t.append(psi)
            
        # print(psi_t)

        a = destroy(N_plot)
        alpha_t = expect(a, psi_t)
        

        print("original_alpha_t")
        print(alpha_t[-1])
        
       
        # print()
    combined_alpha_t_real = []
    combined_alpha_t_imag = []
    labels = []

    for file in amp_json:
        pulse_num = file.split('pulse')[1].split('amp')[0]
        with open(file, 'r') as f:
            data = json.load(f)
        all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]
        N_plot = len(all_transition_amplitudes[0])
        psi_t = []
        for t, cn_t in enumerate(all_transition_amplitudes):
            psi = basis(N_plot, 0)*0
            for n, cn in enumerate(cn_t):
                psi += cn*basis(N_plot, n)
            psi_t.append(psi)
        a = destroy(N_plot)
        alpha_t = expect(a, psi_t)
        combined_alpha_t_real.append(alpha_t.real)
        combined_alpha_t_imag.append(alpha_t.imag)
        labels.append(f"pulse {pulse_num}")

    print("combined_alpha_t_real")
    print(combined_alpha_t_real)
    print("combined_alpha_t_imag")
    print(combined_alpha_t_imag)

    # Combined plot
    plt.figure(figsize=(10, 6))
    for i in range(len(combined_alpha_t_real)):
        plt.plot(combined_alpha_t_real[i], combined_alpha_t_imag[i],
             marker='o', linestyle='none', label=labels[i])
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Combined alpha_t trajectories for all pulses")
    # plt.legend()
    plt.grid(True)
    plt.show()






def generate_amp_json(prefix, state, num_pulses, amp, q, n):
    """
    Generate filenames for amplitudes JSON files.
    Example: amplitudes_e_chg_pulse0amp0.025000_q18_n8_p192.json
    """
    return [
        f"{prefix}_{state}_chg_pulse{i}amp{amp:.6f}_q{q}_n{n}_p{num_pulses}.json"
        # f"{prefix}_{state}_chg_pulse{i}amp{amp:.6f}_q{q}_n{n}.json"
        for i in range(num_pulses)
    ]

def generate_amp_jc_json(prefix, amp, q, nc, nq, num_pulses):
    """
    Generate filenames for amplitudes JSON files.
    Example: amplitudes_e_chg_pulse0amp0.025000_q18_n8_p192.json
    """
    return f"{prefix}_{amp:.6f}_q{q}_nc{nc}_nq{nq}_p{num_pulses}.json"
    

# show only half of the amplitudes

def show_amp(amp_json, amp, q_max, sq, chi):
    """
    Show transition amplitudes from single JSON files.
    """
    if not amp_json:
        warn("No amplitude JSON files provided.")
        return
    # Check if the files exist
    if not Path(amp_json).is_file():
        warn(f"File {amp_json} does not exist.")
        return
    

    with open(amp_json, 'r') as f:
            data = json.load(f)

    all_transition_amplitudes = []

    num_pulses = len(data)
    num_inner_time = len(data[0])
    num_cavity = len(data[0][0])
    print(f"num_pulses: {num_pulses}, num_inner_time: {num_inner_time}, num_cavity: {num_cavity}")
    
    dir_path = f"plots_folder_base/a{amp}_qmax{q_max}_num_p{num_pulses}_nc{num_cavity}_sq{sq}"
    os.makedirs(dir_path,  exist_ok=True)
    
    # all_transition_amplitudes = [[[[[complex(amp['real'], amp['imag']) for amp in cavity] 
    #                             for cavity in inner_time]
    #                             for inner_time in qubit] 
    #                             for qubit in pulses] 
    #                             for pulses in data]
    for pulse in data:
        for t_index in range(len(pulse[0])):  # inner_time loop (assume same length for all qubits)
            time_step = []
            for cavity in pulse[t_index]:
                time_step.append(complex(cavity['real'], cavity['imag']))
            all_transition_amplitudes.append(time_step)


    # nc = len(all_transition_amplitudes[0][0])
    # # print(N_plot)
    # nq = len(all_transition_amplitudes[0])
    # # print(q_plot)

    print("finished loading amplitudes")
    # create Qobj for each time frame in the transition amplitudes
    psi_t = []
    for cn_t in all_transition_amplitudes:
        psi = 0
        for q in range(nq):
            for c in range(nc):
                cn = cn_t[q][c]
                psi += cn * tensor(basis(nq, q), basis(nc, c))
        psi_t.append(psi)

    proj_g = basis(nq, 0) * basis(nq, 0).dag()   # projector onto ground
    proj_e = basis(nq, 1) * basis(nq, 1).dag()   # projector onto excited
    
    P_g = tensor(proj_g, qeye(nc))
    P_e = tensor(proj_e, qeye(nc))


    cavity_states_cond_g = []
    cavity_states_cond_e = []

    for psi in psi_t:
        # Project onto qubit ground
        projected_g = P_g * psi
        norm_g = projected_g.norm()
        if norm_g > 1e-10:  # Avoid division by zero
            projected_g = projected_g / norm_g
            cavity_state_g = projected_g.ptrace(1)  # trace out qubit
            cavity_states_cond_g.append(cavity_state_g)

        # Project onto qubit excited
        projected_e = P_e * psi
        norm_e = projected_e.norm()
        if norm_e > 1e-10:
            projected_e = projected_e / norm_e
            cavity_state_e = projected_e.ptrace(1)  # trace out qubit
            cavity_states_cond_e.append(cavity_state_e)
    
    a = destroy(nc)
    alpha_t_cond_g = [expect(a, rho) for rho in cavity_states_cond_g]
    alpha_t_cond_e = [expect(a, rho) for rho in cavity_states_cond_e]

    pop_g = expect(P_g, psi_t)
    pop_e = expect(P_e, psi_t)

    print("loading plots", dir_path)

    # Plot pop
    plt.figure(figsize=(10, 6))
    plt.plot(pop_g, label='Ground State Population', color='blue')
    plt.plot(pop_e, label='Excited State Population', color='red')
    plt.xlabel("Time Step")
    plt.ylabel("Population")
    plt.title("Qubit State Populations")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "populations.png"))
    plt.show()




    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(alpha_t_cond_g, label='Ground State', color='blue')
    plt.plot(alpha_t_cond_e, label='Excited State', color='red')
    plt.xlabel("Time Step")
    plt.ylabel("Cavity Amplitude (alpha_t)")
    plt.title("Cavity Amplitude Trajectories")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "cavity_amplitude_trajectories.png"))
    plt.show()

    
    # plt.figure(figsize=(10, 6))
    # plt.plot(alpha_t_cond_g[:-2000], label='Ground State', color='blue')
    # plt.plot(alpha_t_cond_e[:-2000], label='Excited State', color='red')
    # plt.xlabel("Time Step")
    # plt.ylabel("Cavity Amplitude (alpha_t)")
    # plt.title("Cavity Amplitude Trajectories")
    # plt.legend()
    # plt.grid(True)
    # plt.savefig(os.path.join(dir_path, "avity_amplitude_trajectories.png"))
    # plt.show()


    # plot the real and imaginary parts of alpha_t_cond_g over phase space
    plt.figure(figsize=(10, 6))
    plt.plot(np.real(alpha_t_cond_g), np.imag(alpha_t_cond_g), label='Ground State', color='blue')
    # Mark start, middle, and end points
    start_idx = 0
    mid_idx = len(alpha_t_cond_g) // 2
    end_idx = len(alpha_t_cond_g) - 1
    plt.scatter(np.real(alpha_t_cond_g[start_idx]), np.imag(alpha_t_cond_g[start_idx]), color='green', marker='o', s=100, label='Start')
    plt.scatter(np.real(alpha_t_cond_g[mid_idx]), np.imag(alpha_t_cond_g[mid_idx]), color='orange', marker='s', s=100, label='Middle')
    plt.scatter(np.real(alpha_t_cond_g[end_idx]), np.imag(alpha_t_cond_g[end_idx]), color='red', marker='*', s=150, label='End')
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Cavity Amplitude Phase Space Trajectories")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "phase_space_g.png"))
    plt.show()

        
    # plot the real and imaginary parts of alpha_t_cond_e over phase space
    plt.figure(figsize=(10, 6))
    plt.plot(np.real(alpha_t_cond_e), np.imag(alpha_t_cond_e), label='Excited State', color='red')
    # Mark start, middle, and end points
    start_idx = 0
    mid1_idx = len(alpha_t_cond_e) // 4
    mid2_idx = len(alpha_t_cond_e) // 2
    mid3_idx = len(alpha_t_cond_e) * 3 // 4
    end_idx = len(alpha_t_cond_e) - 1
    plt.scatter(np.real(alpha_t_cond_e[start_idx]), np.imag(alpha_t_cond_e[start_idx]), color='green', marker='o', s=100, label='Start')
    plt.scatter(np.real(alpha_t_cond_e[mid1_idx]), np.imag(alpha_t_cond_e[mid1_idx]), color='orange', marker='1', s=100, label='qurter')
    plt.scatter(np.real(alpha_t_cond_e[mid2_idx]), np.imag(alpha_t_cond_e[mid2_idx]), color='orange', marker='2', s=100, label='Middle')
    plt.scatter(np.real(alpha_t_cond_e[mid3_idx]), np.imag(alpha_t_cond_e[mid3_idx]), color='orange', marker='3', s=100, label='3qurters')
    plt.scatter(np.real(alpha_t_cond_e[end_idx]), np.imag(alpha_t_cond_e[end_idx]), color='red', marker='*', s=150, label='End')
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Cavity Amplitude Phase Space Trajectories")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "phase_space_e.png"))
    plt.show()


    # #  plot the real and imaginary parts of alpha_t_cond_g over phase space
    # #color every x seconds in different color
    # plt.figure(figsize=(10, 6))
    # time_interval = num_inner_time  # Change color every x time steps
    # colors = ['red', 'blue', 'green', 'purple']
    # for i in range(len(alpha_t_cond_e)//4):
    #     start_idx = i * time_interval
    #     end_idx = min((i + 1) * time_interval, len(alpha_t_cond_e))
    #     plt.plot(np.real(alpha_t_cond_e[start_idx:end_idx]), np.imag(alpha_t_cond_e[start_idx:end_idx]), color=colors[i%4])
    # plt.xlabel("Re(alpha_t)")
    # plt.ylabel("Im(alpha_t)")
    # plt.title("Cavity Amplitude Phase Space Trajectories with Time Segments")
    # plt.colorbar(plt.cm.ScalarMappable(cmap='viridis'), label='Time Segments')
    # plt.grid(True)
    # plt.savefig(os.path.join(dir_path, "phase_space_e_sample.png"))
    # plt.show()

    if(sq == 1):
        alpha_for_disp = alpha_t_cond_e
    if(sq == 0):
        alpha_for_disp = alpha_t_cond_g
    
    # plot first quarter of the real and imaginary parts of alpha_for_disp over phase space
    plt.figure(figsize=(10, 6))
    start_idx = 0
    end_idx = len(alpha_for_disp) // 4
    plt.plot(np.real(alpha_for_disp[start_idx:end_idx]), np.imag(alpha_for_disp[start_idx:end_idx]), color='purple')
    #mark statr and end
    plt.scatter(np.real(alpha_for_disp[start_idx]), np.imag(alpha_for_disp[start_idx]), color='green', marker='o', s=100, label='Start')
    plt.scatter(np.real(alpha_for_disp[end_idx-1]), np.imag(alpha_for_disp[end_idx-1]), color='red', marker='*', s=150, label='End')
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Cavity Amplitude Phase Space Trajectories (First Quarter)")
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "phase_space_first_quarter.png"))
    plt.show()


    # plot second quarter of the real and imaginary parts of alpha_for_disp over phase space
    plt.figure(figsize=(10, 6))
    start_idx = len(alpha_for_disp) // 4
    end_idx = len(alpha_for_disp) // 2
    plt.plot(np.real(alpha_for_disp[start_idx:end_idx]), np.imag(alpha_for_disp[start_idx:end_idx]), color='blue')
    # mark statr and end
    plt.scatter(np.real(alpha_for_disp[start_idx]), np.imag(alpha_for_disp[start_idx]), color='green', marker='o', s=100, label='Start')
    plt.scatter(np.real(alpha_for_disp[end_idx-1]), np.imag(alpha_for_disp[end_idx-1]), color='red', marker='*', s=150, label='End')
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Cavity Amplitude Phase Space Trajectories (Second Quarter)")
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "phase_space_second_quarter.png"))
    plt.show()

    # plot third quarter of the real and imaginary parts of alpha_for_disp over phase space
    plt.figure(figsize=(10, 6))
    start_idx = len(alpha_for_disp) // 2
    end_idx = len(alpha_for_disp) * 3 // 4
    plt.plot(np.real(alpha_for_disp[start_idx:end_idx]), np.imag(alpha_for_disp[start_idx:end_idx]), color='green')
    # mark statr and end
    plt.scatter(np.real(alpha_for_disp[start_idx]), np.imag(alpha_for_disp[start_idx]), color='green', marker='o', s=100, label='Start')
    plt.scatter(np.real(alpha_for_disp[end_idx-1]), np.imag(alpha_for_disp[end_idx-1]), color='red', marker='*', s=150, label='End')
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Cavity Amplitude Phase Space Trajectories (Third Quarter)")
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "phase_space_third_quarter.png"))
    plt.show()

    # plot fourth quarter of the real and imaginary parts of alpha_for_disp over phase space
    plt.figure(figsize=(10, 6))
    start_idx = len(alpha_for_disp) * 3 // 4
    end_idx = len(alpha_for_disp)
    plt.plot(np.real(alpha_for_disp[start_idx:end_idx]), np.imag(alpha_for_disp[start_idx:end_idx]), color='red')
    # mark statr and end
    plt.scatter(np.real(alpha_for_disp[start_idx]), np.imag(alpha_for_disp[start_idx]), color='green', marker='o', s=100, label='Start')
    plt.scatter(np.real(alpha_for_disp[end_idx-1]), np.imag(alpha_for_disp[end_idx-1]), color='red', marker='*', s=150, label='End')
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Cavity Amplitude Phase Space Trajectories (Fourth Quarter)")
    plt.grid(True)
    plt.savefig(os.path.join(dir_path, "phase_space_fourth_quarter.png"))
    plt.show()
    print("finished plots")

    
    



amplitudes = [0.0000001]
q_max_values = [4,7]
num_pulses_list = [960]
nc_list = [3, 4]
chi_list = [-2.79e-6]
sq_list = [0,1]

# amp_json = generate_amp_jc_json("jc_runs-neg-alpha/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json("jc_no_amp/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json("jc_runs/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json("jc_runs/a1e-05_q6_p6400/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json(f"jc_runs/a2e-07_q{q}_p{num_pulses}/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json(f"jc_runs/a0.0001_q{q}_p{num_pulses}/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json("jc_runs_py_scrupt/a0.0001_q5_p192/amplitudes",  amp, q, nc, nq, num_pulses)
# amp_json = generate_amp_jc_json("jc_tests/amplitudes",  amp, q, nc, nq, num_pulses)
# print(amp_json)
# show_amp(amp_json, amp, q)

# for amp in amplitudes:
#     for q in q_max_values:
#         for num_pulses in num_pulses_list:
#             for nc in nc_list:
#                 for chi in chi_list:
#                     for sq in sq_list:
#                         # amp_json = generate_amp_jc_json(f"fresh/0_0/a{amp}_q{q}_p{num_pulses}/amplitudes",  amp, q, nc, nq, num_pulses)
#                         # amp_json = generate_amp_jc_json(f"jc_fresh2/a1e-07_q{q}_p{num_pulses}_chi{chi}_sq{sq}_delta{delta}/amplitudes",  amp, q, nc, nq, num_pulses)
#                         amp_json = generate_amp_jc_json(f"jc_base/a1e-07_q{q}_p{num_pulses}_sq{sq}/amplitudes",  amp, q, nc, nq, num_pulses)
#                         print(amp_json)
#                         show_amp(amp_json, amp, q, sq, chi)

num_pulses = 192
amp = 0.002
# amp = 0.0001
q = 4
n = 4
nc = 5
nq = 2
amp_json = generate_amp_json("sin_test3/amplitudes", "e", num_pulses, amp, q, n)
# print(amp_json)
show_amp_json_dispersive(amp_json)
amp_json = generate_amp_json("sin_test3/amplitudes", "g", num_pulses, amp, q, n)
show_amp_json_dispersive(amp_json)