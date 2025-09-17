import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
from function import Gaussian_square
#################################################################################################
EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
flat, sig, phi_0, time_step, Delta, rabi_rate = 47.5836, 40.8315, 0, 0.02, 0.01, 0.4741
evals = fluxonium01.eigenvals(evals_count=N_q)# - fluxonium01.eigenvals(evals_count=N_q)[0]
tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
H = qt.Qobj(np.diag(evals))
n_opr = qt.Qobj(fluxonium01.matrixelement_table('n_operator', evals_count=N_q))
psi0 = qt.basis(N_q, 0)
e_ops_list = [qt.ket2dm(qt.basis(N_q, 0)), 
            qt.ket2dm(qt.basis(N_q, 1)), 
            qt.ket2dm(qt.basis(N_q, 2)), 
            qt.ket2dm(qt.basis(N_q, 3))]
Rabi_03 =  0.02 * np.pi * 2
Rabi_23 = Rabi_03 * rabi_rate
amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)
#################################################################################################
print(results.expect[0][-1], results.expect[2][-1], results.expect[3][-1])