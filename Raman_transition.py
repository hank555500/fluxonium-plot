import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt


##################################################################################################
def Gaussian_square(freq, flat, sig, phi_0=0, num_sig=2, time_step=None, norm=True):
    """
    total duration = flat+sig*num_sig*2
    norm: normalize gaussian edge to 0 and 1
    unit: ns
    """
    if time_step==None:
        time_step = 0.1/freq
    t_list = np.linspace(0, flat+sig*num_sig*2, int((flat+sig*num_sig*2)/time_step)+1)
    pulse = np.ones(len(t_list))
    num_step = int(sig*num_sig/time_step)
    gauss = np.exp(-(t_list-t_list[0])**2/(2*sig**2))
    edge = gauss[:num_step]

    if norm:
        edge = (edge-edge.min()) / (edge.max()-edge.min())
    pulse[         :num_step] = edge[::-1]
    pulse[-num_step:        ] = edge

    return t_list, pulse*np.cos(2*np.pi*freq*t_list+phi_0)
##################################################################################################
###### Basic index ######
##################################################################################################
#for fluxonium
EJ, EC, EL, phi_eA, N_q = 5.0, 1.8, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
##################################################################################################
#for gaussian function
flat, sig, phi_0, time_step, Delta = 200, 34.3, 0, 0.02, 0.0 
evals = fluxonium01.eigenvals(evals_count=N_q)# - fluxonium01.eigenvals(evals_count=N_q)[0]
tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)

##################################################################################################
# Hamiltonian 4*4 matrix
H = qt.Qobj(np.diag(evals))
# n_operator 4*4 matrix
n_opr = qt.Qobj(fluxonium01.matrixelement_table('n_operator', evals_count=N_q))
#begin state
psi0 = qt.basis(N_q, 0)
##################################################################################################
Rabi_03, Rabi_23 = 0.01 * np.pi * 2, 0.01 * np.pi * 2
amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
##################################################################################################
results0 = qt.sesolve(H_evo, psi0, tlist, e_ops = qt.ket2dm(qt.basis(N_q, 0)))
results1 = qt.sesolve(H_evo, psi0, tlist, e_ops = qt.ket2dm(qt.basis(N_q, 1)))
results2 = qt.sesolve(H_evo, psi0, tlist, e_ops = qt.ket2dm(qt.basis(N_q, 2)))
results3 = qt.sesolve(H_evo, psi0, tlist, e_ops = qt.ket2dm(qt.basis(N_q, 3)))
pi_pulse_index = np.argmax(results2.expect[0])
half_pi_pulse_index = np.argmin(np.abs(results2.expect[0]- results2.expect[0][pi_pulse_index] / 2))
print(results0.expect[0][half_pi_pulse_index] + results1.expect[0][half_pi_pulse_index] + results2.expect[0][half_pi_pulse_index] + results3.expect[0][half_pi_pulse_index])
# plt.plot(tlist, results2.expect[0])
# plt.scatter(tlist[half_pi_pulse_index], results2.expect[0][half_pi_pulse_index], color='red')
# plt.show()
##################################################################################################
# e_ops_list = [qt.ket2dm(qt.basis(N_q, 0)), 
#               qt.ket2dm(qt.basis(N_q, 1)), 
#               qt.ket2dm(qt.basis(N_q, 2)), 
#               qt.ket2dm(qt.basis(N_q, 3))]
# results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)
# qt.plot_expectation_values(results)
# plt.show()
##################################################################################################
# for i in range(N_q):
#     results = qt.sesolve(H_evo, psi0, tlist, e_ops = qt.ket2dm(qt.basis(N_q, i)))
#     qt.plot_expectation_values(results)
#     plt.grid()
#     plt.show()
##################################################################################################
