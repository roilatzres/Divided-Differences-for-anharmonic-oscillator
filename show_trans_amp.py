import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.signal import gaussian, square
from warnings import warn
from qutip import *


# amp_json = [ 'amplitudes_e_chg_amp0.020000_q13.json',  'amplitudes_g_chg_amp0.020000_q13.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.010000_q13.json',  'amplitudes_g_chg_amp0.010000_q13.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.100000_q12.json',  'amplitudes_g_chg_amp0.100000_q12.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.100000_q11.json',  'amplitudes_g_chg_amp0.100000_q11.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.200000_q12.json',  'amplitudes_g_chg_amp0.200000_q12.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.020000_q12.json',  'amplitudes_g_chg_amp0.020000_q12.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.005000_q8.json',  'amplitudes_g_chg_amp0.005000_q8.json' ]
# amp_json = [ 'amplitudes_e_chg_amp0.005000_q10.json',  'amplitudes_g_chg_amp0.005000_q10.json' ]
####################################################################################################


# amp_json = [ 'amplitudes_g_chg_amp0.002000_q12.json',  'amplitudes_e_chg_amp0.002000_q12.json' ]
# amp_json = [ 'amplitudes_g_chg_amp0.001000_q10.json',  'amplitudes_e_chg_amp0.001000_q10.json' ]
# amp_json = [ 'amplitudes_g_chg_amp0.002000_q10.json',  'amplitudes_e_chg_amp0.002000_q10.json' ]
# amp_json = [ 'mid_amplitudes_e_chg_amp0.002000_q10.json',  'mid_amplitudes_g_chg_amp0.002000_q10.json' ]
# amp_json = [ 'mid_amplitudes_g_chg_amp0.002000_q10.json' ]
# amp_json = [ 'mid_amplitudes_e_chg_amp0.002000_q10.json' ]
# amp_json = [ 'mid_amplitudes_g_chg_amp0.002000_q10.json', 'amplitudes_g_chg_amp0.002000_q10.json' ]
# amp_json = [ 'mid_amplitudes_g_chg_amp-0.002000_q10.json', 'amplitudes_g_chg_amp-0.002000_q10.json' ]


# amp_json = [ 'mid_amplitudes_e_chg_amp0.002000_q10.json', 'base_amplitudes_e_chg_amp0.002000_q10.json', 'amplitudes_e_chg_amp0.002000_q10.json' ]
amp_json = [ 'mid_amplitudes_g_chg_amp-0.002000_q10.json', 
            '0_base_amplitudes_g_chg_amp-0.002000_q10.json',
            '1_base_amplitudes_g_chg_amp-0.002000_q10.json',
            '2_base_amplitudes_g_chg_amp-0.002000_q10.json',
            '3_base_amplitudes_g_chg_amp-0.002000_q10.json',
            '4_base_amplitudes_g_chg_amp-0.002000_q10.json',
            '5_base_amplitudes_g_chg_amp-0.002000_q10.json',
            'amplitudes_g_chg_amp-0.002000_q10.json' ]
amp_json = [ 'mid_amplitudes_g_chg_amp0.002000_q10.json', 
            '0_base_amplitudes_g_chg_amp0.002000_q10.json',
            '1_base_amplitudes_g_chg_amp0.002000_q10.json',
            '2_base_amplitudes_g_chg_amp0.002000_q10.json',
            '3_base_amplitudes_g_chg_amp0.002000_q10.json',
            '4_base_amplitudes_g_chg_amp0.002000_q10.json',
            '5_base_amplitudes_g_chg_amp0.002000_q10.json',
            'amplitudes_g_chg_amp0.002000_q10.json' ]
amp_json = [ 'mid_amplitudes_e_chg_amp0.020000_q18.json', 
            'amplitudes_e_chg_amp0.020000_q18.json' ,
            'amplitudes_g_chg_amp0.020000_q18.json' ,
            'amplitudes_g_chg_amp0.020000_q18.json' ]

########################################################################
# all_amp_g_json = []
# all_amp_e_json = []

# # create amp_json name files for different amplitude and different q
# amplitudes = [ 0.001000, 0.002000, 0.005000, 0.010000, 0.020000, 0.050000, 0.060000, 0.100000, 0.200000]

# for amp in amplitudes:
#     for q  in range(8, 14):
#         amp_e_json = 'amplitudes_e_chg_amp'+str(amp)+'q'+str(q)+'.json'
#         amp_g_json = 'amplitudes_g_chg_amp'+str(amp)+'q'+str(q)+'.json'
#         all_amp_e_json.append(amp_e_json)
#         all_amp_g_json.append(amp_g_json)
#########################################################################


# show only half of the amplitudes


for file in amp_json:
# for i, file in enumerate(amp_json):
# for file in all_amp_g_json:
    # load transition amplitudes

    with open(file, 'r') as f:
        data = json.load(f)

    all_transition_amplitudes = [[complex(d['real'], d['imag']) for d in sub_data] for sub_data in data]
    print(len(all_transition_amplitudes))
    
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
    plt.figure()
    plt.plot(alpha_t.real, alpha_t.imag)
    if '_e_' in file:
        plt.title('e')
    else:
        plt.title('g')

    # # output all alpha_t to a file
    # with open('alpha_t.json', 'w') as f:
    #     json.dump(alpha_t, f)

    # #print all alpha_t to txt file
    # with open('alpha_t.txt', 'w') as f:
    #     for item in alpha_t:
    #         f.write("%s\n" % item)
    
    #print last alpha_t
    print()
    print(file)
    print("original_alpha_t")
    print(alpha_t[-1])
    # print abs square of alpha_t
    print("abs_square_alpha_t")
    print(np.abs(alpha_t[-1])**2 )


    # print("new alpha_t")
    # print(alpha_t[-1])

    # print()