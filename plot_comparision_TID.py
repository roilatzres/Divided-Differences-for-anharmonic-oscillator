import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.signal import gaussian, square
from warnings import warn
from qutip import *
from itertools import product



# load alph_t data from json
def load_alpha_t(file):
    with open(file, 'r') as f:
        data = json.load(f)

    g_abs_square = data['e_abs_square']
    # ensure we always return a 1D numpy array (handles scalar or list in json)
    return np.atleast_1d(np.array(g_abs_square))

#load all abs^2 data for all amplitudes in single plot with ampltiude x axis
def plot_all_amp_antySym(amp_list, ax=None, color='C0', marker='o'):
    # draw base (experimental) points on given axes (no show)
    if ax is None:
        _, ax = plt.subplots()
    for a in amp_list:
        file_e = f'data_TID/antisym_results/amp_{a:.4f}.json'
        print(file_e)
        try:
            g_abs_square = load_alpha_t(file_e)
        except FileNotFoundError:
            warn(f"file not found: {file_e}")
            continue
        g_abs_square = np.atleast_1d(np.array(g_abs_square))
        x = np.repeat(a, g_abs_square.size)
        ax.plot(x, g_abs_square, linestyle='--', marker=marker, color=color,
                markeredgecolor='k', markersize=6, label=None)
    return ax


def load_alpha_t_sim(file):
    combined_alpha_t_real = []
    combined_alpha_t_imag = []
    labels = []

    with open(file, 'r') as f:
        data = json.load(f)
    for pulse_idx, pulse in enumerate(data):
        if not pulse:
            continue
        # loop over all qubits present in this pulse
        for qubit_idx, qubit_data in enumerate(pulse):
            if not qubit_data:
                continue

            # qubit_data: list of times; each time is list of targets {real,imag}
            all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in time] for time in qubit_data]
            if not all_transition_amplitudes or not all_transition_amplitudes[0]:
                continue

            N_plot = len(all_transition_amplitudes[0])
            psi_t = []
            for cn_t in all_transition_amplitudes:
                psi = basis(N_plot, 0) * 0
                for n, cn in enumerate(cn_t):
                    psi += cn * basis(N_plot, n)
                psi_t.append(psi)

            a = destroy(N_plot)
            alpha_t = expect(a, psi_t)
            combined_alpha_t_real.append(alpha_t.real)
            combined_alpha_t_imag.append(alpha_t.imag)
            labels.append(f"pulse {pulse_idx} qubit {qubit_idx}")

    # return abs squared of last point
    final_real = combined_alpha_t_real[-1][-1]
    final_imag = combined_alpha_t_imag[-1][-1]
    abs_square = final_real**2 + final_imag**2
    return abs_square
    # all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]

    # N_plot = len(all_transition_amplitudes[0])

    # psi_t = []
    # for t,cn_t in enumerate(all_transition_amplitudes):
    #     psi = basis(N_plot, 0)*0
    #     for n,cn in enumerate(cn_t):
    #         psi += cn*basis(N_plot, n)
    #     psi_t.append(psi)

    # a = destroy(N_plot)
    # alpha_t = expect(a, psi_t)
    # # ensure numpy array
    # return np.atleast_1d(np.array(alpha_t))


def plot_all_amp_sim(amp_list, q_list, n=8, include_base=True):
    fig, ax = plt.subplots()


    # draw simulation points: color fixed per q, starred and bold
    for j, qv in enumerate(q_list):
        color = f"C{j % 10}"
        first = True
        for a in amp_list:
            a_check = a*100
            #round a_cgeck to the nearest integer to match filenames
            a_check = round(a_check)
            print(a_check)
            if(a_check == 10):
                a=0.1
                file = f'data_TID/const_pulse_a{a:.1f}/a{a:.1f}_q{qv}_n100/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n100_p192.json'
            elif(a_check == 12):
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n170/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n170_p192.json'
            elif(a_check == 14):
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n200/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n200_p192.json'
            elif(a_check == 16):
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n300/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n300_p192.json'
            elif(a_check == 18):
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n300/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n300_p192.json'
            elif(a_check == 6):
                a = 0.06
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n100/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n100_p192.json'
            elif(a_check == 7):
                a = 0.07
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n10/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n100_p192.json'
            elif(a_check == 8):
                a = 0.08
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n100/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n100_p192.json'
            else:
                file = f'data_TID/const_pulse_a{a:.2f}/a{a:.2f}_q{qv}_n10/amplitudes_all_pulses_amp{a:.6f}_q{qv}_n10_p192.json'
            
            print(file)
            print(qv, a)
            try:
                e_abs_square = load_alpha_t_sim(file)
            except FileNotFoundError:
                warn(f"file not found: {file}")
                continue
            # alpha_t = np.atleast_1d(np.array(alpha_t))
            # g_abs_square = np.abs(alpha_t)**2
            x = np.repeat(a, e_abs_square.size)

            label = f'q={qv}' if first else None
            # starred, larger and thicker marker edge for bold appearance
            ax.plot(x, e_abs_square, linestyle='-', marker='*', color=color,
                    markeredgecolor='k', markersize=10, markeredgewidth=1.5, label=label)
            first = False

    # optional: draw base/experimental points first (same axes)
    if include_base:
        plot_all_amp_antySym(amp_list, ax=ax, color='0.6', marker='o')
        # add a legend proxy for base points
        ax.plot([], [], linestyle='none', marker='o', color='0.6', markeredgecolor='k', markersize=6, label='base')

    ax.set_xlabel('Amplitude')
    ax.set_ylabel('|α|²     ', rotation=90)
    # ax.yaxis.set_rotation(90)
    ax.set_title('Transition Amplitudes for different drive amplitudes')
    ax.legend()
    plt.show()


# q_list = list(range(5, 12))
q_list = [6,7,8,9]
amp = np.arange(0.02, 0.17, 0.02)
amp = np.append([0.01], amp)
print(amp)
# plot_all_amp_antySym(amp)
plot_all_amp_sim(amp, q_list, n=4)
