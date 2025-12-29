# jc_drive_one_period.py
# Simulate a JC model: qubit + oscillator, with the oscillator driven resonantly
# by a sine waveform for exactly one period.

import numpy as np
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, qeye, destroy, num, sigmap, sigmam, sigmaz,
    mesolve, expect, Options
)

# -----------------------
# Parameters (ħ = 1 units)
# -----------------------
Ncav   = 3         # cavity Hilbert-space truncation
fr     = 7.0              # cavity frequency in Hz (sets the time unit)
wr     = 2*np.pi*fr       # cavity angular frequency
fq     = 5.0              # qubit transition frequency in Hz (can set = fr for resonance)
wq     = 2*np.pi*fq


Delta = -wr+wq
chi = -279e-6 * 2 * np.pi
Ec = 200e-3 * 2 * np.pi
# g      = np.sqrt(chi*Delta*(Delta-Ec)/Ec)
g      = np.sqrt(np.abs(chi*Delta))
Adr    = 0.01e-3 * 2 * np.pi         # drive amplitude (couples to a + a^\dagger)

# delta = 1e-3 * 2 * np.pi
delta = 1/192*2*np.pi

# shift = chi
shift = 0

# -----------------------
# Operators
# -----------------------
a   = destroy(Ncav)                       # cavity annihilation
I_c = qeye(Ncav)
I_q = qeye(2)

sm  = sigmam()
sp  = sigmap()
sz  = sigmaz()

# Tensor to joint space
a_t   = tensor(a, I_q)
adag_t= a_t.dag()
sm_t  = tensor(I_c, sm)
sp_t  = tensor(I_c, sp)
sz_t  = tensor(I_c, sz)

# Number operator and excited projector for observables
n_op  = tensor(num(Ncav), I_q)
Pe_op = tensor(I_c, basis(2,1)*basis(2,1).dag())  # |e><e|
Pg_op = tensor(I_c, basis(2,0)*basis(2,0).dag())  # |g><g|


# Time grid: exactly one period of the cavity
# -----------------------
T = 2 * np.pi / delta                 # one period of the drive/cavity in seconds
tlist = np.linspace(0.0, T, 10001)

# -----------------------
# Hamiltonian
# -----------------------
# Bare JC Hamiltonian (in lab frame)
H0 = chi*adag_t*a_t

def eps_t(t, args=None):
    # Drive for exactly one period of the cavity at amplitude Adr: Adr * sin(wr t)
    return Adr * np.sin(delta * t) * np.exp(-1j * chi * t*0)
    # return Adr * np.sin(delta * t) 
def eps_t_star(t, args=None):
    # Drive for exactly one period of the cavity at amplitude Adr: Adr * sin(wr t)
    return np.conjugate(Adr * np.sin(delta * t) * np.exp(-1j * chi * t*0))
    # return np.conjugate(Adr * np.sin(delta * t) )

# H = [H0,
H = [ 
     [a_t, eps_t_star], 
     [adag_t, eps_t], 
     [a_t * sp_t, g * np.exp(1j*(Delta-shift)*tlist)], #TODO: check conjugation
     [adag_t * sm_t, g * np.exp(-1j*(Delta-shift)*tlist)]] #TODO: check conjugation

# H = [H0, 
#     #  [a_t, eps_t_star], 
#     #  [adag_t, eps_t], 
#      [a_t * sp_t, g * np.exp(1j*(Delta-shift)*tlist)], #TODO: check conjugation
#      [adag_t * sm_t, g * np.exp(-1j*(Delta-shift)*tlist)]] #TODO: check conjugation
# -----------------------
# Initial state
# -----------------------
# Start in cavity vacuum and qubit ground: |0> ⊗ |g>
psi0 = tensor(basis(Ncav, 0), basis(2, 0))
psi1 = tensor(basis(Ncav, 0), basis(2, 1))

# -----------------------
# Collapse operators (optional)
# -----------------------
# -----------------------


# -----------------------
# Solve
# -----------------------
result0 = mesolve(H, psi0, tlist, c_ops=None, e_ops=[tensor(a,qeye(2)), Pe_op], progress_bar = True)
result1 = mesolve(H, psi1, tlist, c_ops=None, e_ops=[tensor(a,qeye(2)), Pe_op], progress_bar = True)

# alpha_res = result.expect[0] * np.exp(1j*wr*tlist)
alpha_res0 = result0.expect[0]
alpha_res1 = result1.expect[0]

# -----------------------
# Plot
# -----------------------

plt.figure()
plt.plot(alpha_res0.real, alpha_res0.imag, label='g', color='blue')
plt.plot(alpha_res0.real[0], alpha_res0.imag[0], 's')
plt.plot(alpha_res0.real[-1], alpha_res0.imag[-1], 'o')
plt.plot(alpha_res1.real, alpha_res1.imag, label='e', color='red')
plt.plot(alpha_res1.real[-1], alpha_res1.imag[-1], 'o')
plt.plot(alpha_res1.real[0], alpha_res1.imag[0], 's')
plt.legend()


plt.figure()
plt.plot(alpha_res0, label='g', color='blue')
plt.plot(alpha_res1, label='e', color='red')
plt.xlabel('Time (s)')
plt.ylabel('cavity amplitude')
plt.legend()

resultg0 = mesolve(H, psi0, tlist, c_ops=None, e_ops=[tensor(a,qeye(2)), Pg_op])
resultg1 = mesolve(H, psi1, tlist, c_ops=None, e_ops=[tensor(a,qeye(2)), Pg_op])

# alpha_res = result.expect[0] * np.exp(1j*wr*tlist)
alpha_res0g = resultg0.expect[0]
alpha_res1g = resultg1.expect[0]

# -----------------------
# Plot
# -----------------------

plt.figure()
plt.plot(alpha_res0g.real, alpha_res0g.imag, label='g', color='blue')
plt.plot(alpha_res0g.real[0], alpha_res0g.imag[0], 's')
plt.plot(alpha_res0g.real[-1], alpha_res0g.imag[-1], 'o')
plt.plot(alpha_res1g.real, alpha_res1g.imag, label='e', color='red')
plt.plot(alpha_res1g.real[-1], alpha_res1g.imag[-1], 'o')
plt.plot(alpha_res1g.real[0], alpha_res1g.imag[0], 's')
plt.legend()


plt.figure()
plt.plot(alpha_res0g, label='g', color='blue')
plt.plot(alpha_res1g, label='e', color='red')
plt.xlabel('Time (s)')
plt.ylabel('cavity amplitude')
plt.legend()
# -----------------------
# Notes:
# - Set fq = fr to place the qubit on resonance with the cavity.
# - Increase Ncav if the drive populates higher Fock states.
# - To include dissipation, set kappa, gamma1, gamma_phi > 0.
# - For a rotating-frame version, you can move to the drive frame and apply RWA.
# -----------------------