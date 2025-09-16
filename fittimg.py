import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
from function import Gaussian_square
#################################################################################################
EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
flat, sig, phi_0, time_step, Delta = 50.0, 34.3, 0, 0.02, 0.01
evals = fluxonium01.eigenvals(evals_count=N_q)
print(evals)
tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
H = qt.Qobj(np.diag(evals))
n_opr = qt.Qobj(fluxonium01.matrixelement_table('n_operator', evals_count=N_q))
psi0 = qt.basis(N_q, 0)
e_ops_list = [qt.ket2dm(qt.basis(N_q, 0)), 
              qt.ket2dm(qt.basis(N_q, 1)), 
              qt.ket2dm(qt.basis(N_q, 2)), 
              qt.ket2dm(qt.basis(N_q, 3))]
