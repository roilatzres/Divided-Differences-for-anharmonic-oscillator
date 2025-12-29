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

    g_abs_square = data['g_abs_square']
    # ensure we always return a 1D numpy array (handles scalar or list in json)
    return np.atleast_1d(np.array(g_abs_square))

#load all abs^2 data for all amplitudes in single plot with ampltiude x axis
def plot_all_amp_antySym(amp_list, ax=None, color='C0', marker='o'):
    # draw base (experimental) points on given axes (no show)
    if ax is None:
        _, ax = plt.subplots()
    for a in amp_list:
        file_g = f'data_orig/antisym_results/amp_{a:.4f}.json'
        print(file_g)
        try:
            g_abs_square = load_alpha_t(file_g)
        except FileNotFoundError:
            warn(f"file not found: {file_g}")
            continue
        g_abs_square = np.atleast_1d(np.array(g_abs_square))
        x = np.repeat(a, g_abs_square.size)
        ax.plot(x, g_abs_square, linestyle='None', marker=marker, color=color,
                markeredgecolor='k', markersize=6, label=None)
    return ax


def load_alpha_t_sim(file):
    with open(file, 'r') as f:
        data = json.load(f)

    all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]

    N_plot = len(all_transition_amplitudes[0])

    psi_t = []
    for t,cn_t in enumerate(all_transition_amplitudes):
        psi = basis(N_plot, 0)*0
        for n,cn in enumerate(cn_t):
            psi += cn*basis(N_plot, n)
        psi_t.append(psi)

    a = destroy(N_plot)
    alpha_t = expect(a, psi_t)
    # ensure numpy array
    return np.atleast_1d(np.array(alpha_t))


def plot_all_amp_sim(amp_list, q_list, n=8, include_base=True):
    fig, ax = plt.subplots()


    # draw simulation points: color fixed per q, starred and bold
    for j, qv in enumerate(q_list):
        color = f"C{j % 10}"
        first = True
        for a in amp_list:
            a = round(a,4)
            if ((a*1000) % 2) != 0:
                file = f'data_orig/main_base_a{a:.3f}/a{a:.3f}_q{qv}_n{n}/amplitudes_g_chg_amp{a:.6f}_q{qv}_n{n}.json'
            else:
                file = f'data_orig/main_base_a{a:.3f}/a{a:.2f}_q{qv}_n{n}/amplitudes_g_chg_amp{a:.6f}_q{qv}_n{n}.json'
            print(file)
            print(qv, a)
            try:
                alpha_t = load_alpha_t_sim(file)
            except FileNotFoundError:
                warn(f"file not found: {file}")
                continue
            alpha_t = np.atleast_1d(np.array(alpha_t))
            g_abs_square = np.abs(alpha_t)**2
            x = np.repeat(a, g_abs_square.size)

            label = f'q={qv}' if first else None
            # starred, larger and thicker marker edge for bold appearance
            ax.plot(x, g_abs_square, linestyle='None', marker='*', color=color,
                    markeredgecolor='k', markersize=10, markeredgewidth=1.5, label=label)
            first = False

    # optional: draw base/experimental points first (same axes)
    if include_base:
        plot_all_amp_antySym(amp_list, ax=ax, color='0.6', marker='o')
        # add a legend proxy for base points
        ax.plot([], [], linestyle='None', marker='o', color='0.6', markeredgecolor='k', markersize=6, label='base')

    ax.set_xlabel('Amplitude')
    ax.set_ylabel('|α|²     ', rotation=90)
    # ax.yaxis.set_rotation(90)
    ax.set_title('Transition Amplitudes for different drive amplitudes')
    ax.legend()
    plt.show()


# q_list = list(range(5, 12))
q_list = [6,8,10,12,14]
amp = np.linspace(0.015, 0.05, 8)
print(amp)
# plot_all_amp_antySym(amp)
plot_all_amp_sim(amp, q_list, n=4)
