import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.signal import gaussian, square
from warnings import warn
from qutip import *


amp_g_json = [ 'amplitudes_g_chg_amp.json',  'amplitudes_e_chg_amp1.json' ]
# amp_e_json = ['amplitudes_e_allQ.json' , 'amplitudes_e_chg_amp.json']

# for file in amp_e_json:
for file in amp_g_json:
    # load transition amplitudes
    with open(file, 'r') as f:
        data = json.load(f)

    all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]

    # print all_transition_amplitudes
    print(all_transition_amplitudes)

    N_plot = len(all_transition_amplitudes[0])
    print(N_plot)

    psi_t = []
    # create Qobj for each time frame in the transition amplitudes
    for t,cn_t in enumerate(all_transition_amplitudes):
        psi = basis(N_plot, 0)*0
        for n,cn in enumerate(cn_t):
            psi += cn*basis(N_plot, n)
        psi_t.append(psi)
        
    print(psi_t)

    a = destroy(N_plot)
    alpha_t = expect(a, psi_t)
    plt.figure()
    plt.plot(alpha_t.real, alpha_t.imag)
