# %%

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 12:25:07 2021

@author: Eliya
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.signal import gaussian, square
from warnings import warn


def create_padded_gaussian(start, sigma, multiplier, tfinal, timestep = 1, is_cut_tail = True):
    before_zeros = np.zeros(round(start/timestep))
    gauss = gaussian(round((sigma*multiplier)/timestep), sigma/timestep)
    if is_cut_tail: gauss = gauss-np.min(gauss)
    after_zeros  = np.zeros(round((tfinal-start-sigma*multiplier)/timestep))
    padded_gauss = np.append(np.append(before_zeros, gauss), after_zeros)
    return padded_gauss

def create_pi_pulse(start, sigma, multiplier, phase, tfinal, timestep = 1, is_cut_tail = True):
    pulse = create_padded_gaussian(start, sigma, multiplier, tfinal, timestep, is_cut_tail)
    return np.exp(1j*phase) * np.pi / 2 * pulse / np.sum(pulse) / timestep

def create_pi2_pulse(start, sigma, multiplier, phase, tfinal, timestep = 1, is_cut_tail = True):
    pulse = create_padded_gaussian(start, sigma, multiplier, tfinal, timestep, is_cut_tail)
    return np.exp(1j*phase) * np.pi / 4 * pulse / np.sum(pulse) / timestep

def create_fast_condisp(amp, start, sigma, multiplier, detuning, times, omega_shift = 0, timestep = 1, is_cut_tail = True):
    gauss = create_padded_gaussian(start, sigma, multiplier, tfinal, timestep, is_cut_tail)
    waveform_upper = amp * gauss * np.exp(-1j * (times-start-sigma*multiplier/2+0.5) *   detuning)  / timestep / np.sum(gauss)
    waveform_lower = -amp * gauss * np.exp(-1j * (times-start-sigma*multiplier/2+0.5) * (-detuning)) / timestep / np.sum(gauss)
    condisp_pulse = waveform_upper + waveform_lower
    return condisp_pulse * np.exp(-1j*omega_shift*times)


#%%
fontsize=20
ticksize=14
is_progress_bar = True

#%% states
g = basis(2,0)
e = basis(2,1)
I = g + e

Nmax = 1000
a = destroy(Nmax)
vac = basis(Nmax,0)

#%%


kappa = 5e-6*2*np.pi

chi =  - 279e-6 * 2 * np.pi

T1 = 13e3
T2_ramsey = 12e3
T2 = 1/(1/2/T1+1/T2_ramsey)

H = chi * tensor(e.proj(), a.dag()*a)

HdriveQ = 1j  * tensor(qeye(2), a.dag()-a)
HdriveI = tensor(qeye(2), a.dag()+a)

HdriveY = tensor(sigmay(), qeye(Nmax))
HdriveX = tensor(sigmax(), qeye(Nmax))

#%% program

multiplier = 4 
sigma = 48 
# detuning =  chi/100
omega_shift = chi
# amp = 30000

phase_corr = 0
timestep = 1

# first condisp
condisp1_start_time = 0
condisp1_stop_time = multiplier * sigma

tfinal= condisp1_stop_time
detuning = 2 * np.pi / tfinal 
print(int((tfinal-timestep)/2))
times = np.linspace(0, int((tfinal-timestep)/2), int((tfinal/timestep)/2))
opts = Options(store_final_state=True, nsteps = len(times) * 100, max_step = 0.01)

is_cut_tail_antisym = False
amp = 0.001
# condisp_pulse1 = amp *  np.sin((times) * detuning) * np.exp(-1j*omega_shift*times) 
condisp_pulse1 = amp *  np.exp(-1j*omega_shift*times) 
condisp_pulse2 = - condisp_pulse1

# condisp_pulse1 = create_fast_condisp(amp=amp,
#                                      start = 0,
#                                      sigma = sigma,
#                                      multiplier = multiplier,
#                                      detuning = detuning,
#                                      omega_shift = omega_shift,
#                                      times = times, is_cut_tail = is_cut_tail_antisym)

mem_pulses = condisp_pulse1
mem_pulses2 = condisp_pulse2

fig,ax = plt.subplots(1,1)
ax.plot(times, (mem_pulses*np.exp(1j*chi*times)).real)
ax.plot(times, (mem_pulses*np.exp(1j*chi*times)).imag)
plt.title('Frame of e')
fig,ax = plt.subplots(1,1)
ax.plot(times, (mem_pulses).real)
ax.plot(times, (mem_pulses).imag)
plt.title('Frame of g')

#%% waveforms
init_state = tensor(g+e,vac).unit()

e_ops = []
e_ops.append(tensor(g.proj(), a))
e_ops.append(tensor(e.proj(), a))
c_ops = []

plt.pause(0.5)

res = mesolve([H, 
            [HdriveI, mem_pulses.real], [HdriveQ, mem_pulses.imag]], 
            init_state,
            times,
            c_ops = c_ops,
            e_ops = e_ops,
            progress_bar = is_progress_bar,
            options = opts
                )
res2 = mesolve([H, 
            [HdriveI, mem_pulses2.real], [HdriveQ, mem_pulses2.imag]], 
            res.final_state,
            times,
            c_ops = c_ops,
            e_ops = e_ops,
            progress_bar = is_progress_bar,
            options = opts
                )

#%%

xlim = 10
ylim = 10
num_of_pnts = 101
x = np.linspace(-xlim,xlim,num_of_pnts)
y = np.linspace(-ylim,ylim,num_of_pnts)

Nmax_to_plot = 4
focks_to_plot = range(0, Nmax_to_plot)

final_mem_state = res.final_state.ptrace(1).extract_states(focks_to_plot).tidyup()

final_mem_g_state = (tensor(g*g.dag(),qeye(Nmax)) * res.final_state ).ptrace(1).extract_states(focks_to_plot).tidyup()
final_mem_e_state = (tensor(e*e.dag(),qeye(Nmax)) * res.final_state ).ptrace(1).extract_states(focks_to_plot).tidyup()

# plot_wigner(final_mem_state)
plt.plot((res.expect[0].real)*2*np.sqrt(2),(res.expect[0].imag)*2*np.sqrt(2), label = 'g')
plt.plot((res.expect[1].real)*2*np.sqrt(2), (res.expect[1].imag)*2*np.sqrt(2), label = 'e')
plt.legend()

print("amp:", amp) 
print("show_trans_amp: ", amp*2)
# print('g:', res.expect[0][-1].real*2*np.sqrt(2), res.expect[0][-1].imag*2*np.sqrt(2))
# print('g:', res.expect[0][-1].real, res.expect[0][-1].imag)
print("for half time - ", (int(tfinal/2)-1))
print('g:', res.expect[0][int(tfinal/2)-1].real*2*np.sqrt(2), res.expect[0][int(tfinal/2)-1].imag*2*np.sqrt(2))
# print(res.expect[0][-1].real*2*np.sqrt(2), res.expect[0][-1].imag*2*np.sqrt(2))
print('e:', res.expect[1][int(tfinal/2)-1].real*2*np.sqrt(2), res.expect[1][int(tfinal/2)-1].imag*2*np.sqrt(2))
#print abs square of g and e
print("abs_square")
print('g:', np.abs(res.expect[0][int(tfinal/2)-1] *2*np.sqrt(2))**2)
print('e:', np.abs(res.expect[1][int(tfinal/2)-1] *2*np.sqrt(2))**2)    

# print("for half time - ", int(tfinal/2))
# print('g:', res.expect[0][int(tfinal/2)].real*2*np.sqrt(2), res.expect[0][int(tfinal/2)].imag*2*np.sqrt(2))
# # print(res.expect[0][-1].real*2*np.sqrt(2), res.expect[0][-1].imag*2*np.sqrt(2))
# print('e:', res.expect[1][int(tfinal/2)].real*2*np.sqrt(2), res.expect[1][int(tfinal/2)].imag*2*np.sqrt(2))
# #print abs square of g and e
# print("abs_square")
# print('g:', np.abs(res.expect[0][int(tfinal/2)] *2*np.sqrt(2))**2)
# print('e:', np.abs(res.expect[1][int(tfinal/2)] *2*np.sqrt(2))**2)    

# print("for half time - ", (int(tfinal/2) + 1))
# print('g:', res.expect[0][int(tfinal/2)+1].real*2*np.sqrt(2), res.expect[0][int(tfinal/2)+1].imag*2*np.sqrt(2))
# # print(res.expect[0][-1].real*2*np.sqrt(2), res.expect[0][-1].imag*2*np.sqrt(2))
# print('e:', res.expect[1][int(tfinal/2)+1].real*2*np.sqrt(2), res.expect[1][int(tfinal/2)+1].imag*2*np.sqrt(2))
# #print abs square of g and e
# print("abs_square")
# print('g:', np.abs(res.expect[0][int(tfinal/2)+1] *2*np.sqrt(2))**2)
# print('e:', np.abs(res.expect[1][int(tfinal/2)+1] *2*np.sqrt(2))**2)    

# print("without 2sqrt(2)")
# print('g:', res.expect[0][int(tfinal/2)].real, res.expect[0][int(tfinal/2)].imag)
# print('e:', res.expect[1][int(tfinal/2)].real, res.expect[1][int(tfinal/2)].imag)
# print("abs_square")
# print('g:', np.abs(res.expect[0][int(tfinal/2)])**2)
# print('e:', np.abs(res.expect[1][int(tfinal/2)])**2)    



# plot_wigner(final_mem_g_state) after mi us pulse
final_mem_state = res2.final_state.ptrace(1).extract_states(focks_to_plot).tidyup()

final_mem_g_state = (tensor(g*g.dag(),qeye(Nmax)) * res2.final_state ).ptrace(1).extract_states(focks_to_plot).tidyup()
final_mem_e_state = (tensor(e*e.dag(),qeye(Nmax)) * res2.final_state ).ptrace(1).extract_states(focks_to_plot).tidyup()

# plot_wigner(final_mem_state)
plt.plot((res2.expect[0].real)*2*np.sqrt(2),(res2.expect[0].imag)*2*np.sqrt(2), label = 'g')
plt.plot((res2.expect[1].real)*2*np.sqrt(2), (res2.expect[1].imag)*2*np.sqrt(2), label = 'e')
plt.legend()

print("amp:", amp) 
print("show_trans_amp: ", amp*2)
# print('g:', res2.expect[0][-1].real*2*np.sqrt(2), res2.expect[0][-1].imag*2*np.sqrt(2))
# print('g:', res2.expect[0][-1].real, res2.expect[0][-1].imag)
print("for half time - ", (int(tfinal/2)-1))
print('g:', res2.expect[0][int(tfinal/2)-1].real*2*np.sqrt(2), res2.expect[0][int(tfinal/2)-1].imag*2*np.sqrt(2))
# print(res2.expect[0][-1].real*2*np.sqrt(2), res2.expect[0][-1].imag*2*np.sqrt(2))
print('e:', res2.expect[1][int(tfinal/2)-1].real*2*np.sqrt(2), res2.expect[1][int(tfinal/2)-1].imag*2*np.sqrt(2))
#print abs square of g and e
print("abs_square")
print('g:', np.abs(res2.expect[0][int(tfinal/2)-1] *2*np.sqrt(2))**2)
print('e:', np.abs(res2.expect[1][int(tfinal/2)-1] *2*np.sqrt(2))**2)    

# print("for half time - ", int(tfinal/2))
# print('g:', res2.expect[0][int(tfinal/2)].real*2*np.sqrt(2), res2.expect[0][int(tfinal/2)].imag*2*np.sqrt(2))
# # print(res2.expect[0][-1].real*2*np.sqrt(2), res2.expect[0][-1].imag*2*np.sqrt(2))
# print('e:', res2.expect[1][int(tfinal/2)].real*2*np.sqrt(2), res2.expect[1][int(tfinal/2)].imag*2*np.sqrt(2))
# #print abs square of g and e
# print("abs_square")
# print('g:', np.abs(res2.expect[0][int(tfinal/2)] *2*np.sqrt(2))**2)
# print('e:', np.abs(res2.expect[1][int(tfinal/2)] *2*np.sqrt(2))**2)    

# print("for half time - ", (int(tfinal/2) + 1))
# print('g:', res2.expect[0][int(tfinal/2)+1].real*2*np.sqrt(2), res2.expect[0][int(tfinal/2)+1].imag*2*np.sqrt(2))
# # print(res2.expect[0][-1].real*2*np.sqrt(2), res2.expect[0][-1].imag*2*np.sqrt(2))
# print('e:', res2.expect[1][int(tfinal/2)+1].real*2*np.sqrt(2), res2.expect[1][int(tfinal/2)+1].imag*2*np.sqrt(2))
# #print abs square of g and e
# print("abs_square")
# print('g:', np.abs(res2.expect[0][int(tfinal/2)+1] *2*np.sqrt(2))**2)
# print('e:', np.abs(res2.expect[1][int(tfinal/2)+1] *2*np.sqrt(2))**2)    

print("without 2sqrt(2)")
print('g:', res2.expect[0][int(tfinal/2) - 1].real, res2.expect[0][int(tfinal/2) - 1].imag)
print('e:', res2.expect[1][int(tfinal/2 ) - 1].real, res2.expect[1][int(tfinal/2) - 1].imag)
print("abs_square")
print('g:', np.abs(res2.expect[0][int(tfinal/2) - 1])**2)
print('e:', np.abs(res2.expect[1][int(tfinal/2) - 1])**2)    
# plot_wigner(final_mem_g_state)
# plt.plot(res.expect[0].real*np.sqrt(2), res.expect[0].imag*np.sqrt(2))
# plot_wigner(final_mem_e_state)
# plt.plot(res.expect[1].real*np.sqrt(2), res.expect[1].imag*np.sqrt(2))

# %%
