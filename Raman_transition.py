import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
from function import Gaussian_square
###### Basic index ######
##################################################################################################
#for fluxonium
EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
##################################################################################################
#for gaussian function
flat, sig, phi_0, time_step, Delta = 4.9556, 19.1211, 0, 0.02, 0.01
evals = fluxonium01.eigenvals(evals_count=N_q)# - fluxonium01.eigenvals(evals_count=N_q)[0]
tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)

##################################################################################################
# Hamiltonian 4*4 matrix
H = qt.Qobj(np.diag(evals))
# n_operator 4*4 matrix
n_opr = qt.Qobj(fluxonium01.matrixelement_table('n_operator', evals_count=N_q))
#begin state
psi0 = qt.basis(N_q, 2)
##################################################################################################
Rabi_03 = 0.02 * np.pi * 2
Rabi_23 = Rabi_03 * 0.4199
amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
##################################################################################################
e_ops_list = [qt.ket2dm(qt.basis(N_q, 0)), 
              qt.ket2dm(qt.basis(N_q, 1)), 
              qt.ket2dm(qt.basis(N_q, 2)), 
              qt.ket2dm(qt.basis(N_q, 3))]
results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)
# pi_pulse_index = np.argmax(results.expect[2])
# half_pi_pulse_index = np.argmin(np.abs(results.expect[2][:pi_pulse_index]- results.expect[2][pi_pulse_index] / 2))
# plt.plot(tlist, results.expect[2])
# plt.scatter(tlist[pi_pulse_index], results.expect[2][pi_pulse_index], color='black')
# plt.scatter(tlist[half_pi_pulse_index], results.expect[2][half_pi_pulse_index], color='red')
# plt.show()
qt.plot_expectation_values(results)
plt.show()
##################################################################################################
##################################################################################################