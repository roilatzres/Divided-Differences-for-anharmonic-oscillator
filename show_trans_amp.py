import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.signal import gaussian, square
from warnings import warn
from qutip import *


# amp_json =['amplitudes_e_chg_pulse0amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse1amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse2amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse3amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse4amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse5amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse6amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse7amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse8amp0.025000_q18_n8.json',
#            'amplitudes_e_chg_pulse9amp0.025000_q18_n8.json',
#            'amplitudes_g_chg_pulse0amp0.015000_q18_n10.json',
#            'amplitudes_g_chg_pulse1amp0.015000_q18_n10.json',
#            'amplitudes_g_chg_pulse2amp0.015000_q18_n10.json',
#            'amplitudes_g_chg_pulse3amp0.015000_q18_n10.json']

def generate_amp_json(prefix, state, num_pulses, amp, q, n):
    """
    Generate filenames for amplitudes JSON files.
    Example: amplitudes_e_chg_pulse0amp0.025000_q18_n8_p192.json
    """
    return [
        f"{prefix}_{state}_chg_pulse{i}amp{amp:.6f}_q{q}_n{n}_p{num_pulses}.json"
        for i in range(num_pulses)
    ]


# show only half of the amplitudes

def show_amp_json(amp_json):
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
    # for i, file in enumerate(amp_json):
    # for file in all_amp_g_json:
        # load transition amplitudes
        # what is the pulse number in the filename?
        pulse_num = file.split('pulse')[1].split('amp')[0]

        with open(file, 'r') as f:
            data = json.load(f)

        all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]
        # print(len(all_transition_amplitudes))
        
        # # cut all_transition_amplitudes to half
        # all_transition_amplitudes = all_transition_amplitudes[len(all_transition_amplitudes)//2:]
        
        # print all_transition_amplitudes
        # print(all_transition_amplitudes)

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
        # plt.figure()
        # plt.plot(alpha_t.real, alpha_t.imag)
        # if '_e_' in file:
        #     # add pulse number to the title
        #     plt.title(f'e, pulse {pulse_num}')
        #     # plt.title('e')
        # else:
        #     plt.title('g')

        # # output all alpha_t to a file
        # with open('alpha_t.json', 'w') as f:
        #     json.dump(alpha_t, f)

        # #print all alpha_t to txt file
        # with open('alpha_t.txt', 'w') as f:
        #     for item in alpha_t:
        #         f.write("%s\n" % item)
        
        #print last alpha_t
        # print(file)
        # print("original_alpha_t")
        # print(alpha_t[-1])
        # # print abs square of alpha_t
        # print("abs_square_alpha_t")
        # print(np.abs(alpha_t[-1])**2 )


        # print("new alpha_t")
        # print(alpha_t[-1])

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

    # Combined plot
    plt.figure(figsize=(10, 6))
    for i in range(len(combined_alpha_t_real)):
        plt.plot(combined_alpha_t_real[i], combined_alpha_t_imag[i], label=labels[i])
    plt.xlabel("Re(alpha_t)")
    plt.ylabel("Im(alpha_t)")
    plt.title("Combined alpha_t trajectories for all pulses")
    # plt.legend()
    plt.grid(True)
    plt.show()




# num_pulses = 64
# amp = 0.004
# q = 8
# n = 10

num_pulses = 64
amp = 0.045
q = 7
n = 25
amp_json = generate_amp_json("sin_pulse_runs/amplitudes", "e", num_pulses, amp, q, n)
# amp_json = generate_amp_json("new_runs/amplitudes", "e", num_pulses, amp, q, n)
show_amp_json(amp_json)
amp_json = generate_amp_json("sin_pulse_runs/amplitudes", "g", num_pulses, amp, q, n)
# amp_json = generate_amp_json("new_runs/amplitudes", "g", num_pulses, amp, q, n)
show_amp_json(amp_json)
