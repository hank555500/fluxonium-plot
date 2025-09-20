import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
from function import Gaussian_square
###### Basic index ######
##################################################################################################
#for fluxonium
EJ, EC, EL, phi_eA, N_q, in_state = 5.0, 1.5, 0.15, 0.0, 4, 0
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
e_ops_list = [qt.ket2dm(qt.basis(N_q, 0)), 
              qt.ket2dm(qt.basis(N_q, 1)), 
              qt.ket2dm(qt.basis(N_q, 2)), 
              qt.ket2dm(qt.basis(N_q, 3))]
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
# 
c_ops = []
##################################################################################################
# Hamiltonian Time evolution 4*4 matrix
Rabi_03 = 0.02 * np.pi * 2
Rabi_23 = Rabi_03 * 0.4199
def amp_03_func(t, args):
    idx = np.argmin(np.abs(tlist - t))
    return pulse_03[idx] * Rabi_03 / abs(n_opr[0,3])

def amp_23_func(t, args):
    idx = np.argmin(np.abs(tlist - t))
    return pulse_23[idx] * Rabi_23 / abs(n_opr[2,3])

H_evo = [H, [n_opr, amp_03_func], [n_opr, amp_23_func]]
##################################################################################################
props = qt.propagator(H_evo, tlist)
# Begin state
psi0 = qt.fock_dm(N_q, in_state) # density matrix
# print(np.abs(np.array(props[-1].as)))
final = props[-1] * psi0 * props[-1].dag()
for i in range(N_q):
    P_i = qt.ket2dm(qt.basis(N_q, i))
    pop = (P_i * final).tr().real
    print(f"Population at state {i}: {pop:.4f}")


