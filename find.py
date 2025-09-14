import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
#################################################################################################
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
rabi_rate = np.linspace(0.0, 2, 100)
EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
flat, sig, phi_0, time_step, Delta = 200, 34.3, 0, 0.02, 0.0
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
mode = 2
state0pop, state1pop, state2pop, state3pop = [], [], [], []
##################################################################################################

# for i in rabi_rate:
#     Rabi_03 =  0.02 * np.pi * 2
#     Rabi_23 = Rabi_03 * i 
#     amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
#     H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
#     results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)
#     pi_pulse_index = np.argmax(results.expect[mode])
#     half_pi_pulse_index = np.argmin(np.abs(results.expect[mode][:pi_pulse_index] - results.expect[mode][pi_pulse_index] / 2))
#     state0pop.append(results.expect[0][half_pi_pulse_index])
#     state2pop.append(results.expect[2][half_pi_pulse_index])
#     state3pop.append(results.expect[3][half_pi_pulse_index])

#     # n = results.expect[0][half_pi_pulse_index] + results.expect[2][half_pi_pulse_index]
#     # if results.expect[0][half_pi_pulse_index] > 0.46 and results.expect[2][half_pi_pulse_index] > 0.46:
#     #     print(i, n, results.expect[3][half_pi_pulse_index])

# fig, (ax0, ax2, ax3) = plt.subplots(3, 1, sharex=True)
# ax0.plot(rabi_rate, state0pop)
# ax0.set_title('population od 0 stae')
# ax2.plot(rabi_rate, state2pop)
# ax3.plot(rabi_rate, state3pop)

# plt.show()
##################################################################################################
for i in rabi_rate :
    Rabi_03 =  0.02 * np.pi * 2
    Rabi_23 = Rabi_03 * i 
    amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
    H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
    results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)
    state0pop.append(np.array(results.expect[0]))
    state1pop.append(np.array(results.expect[1]))
    state2pop.append(np.array(results.expect[2]))
    state3pop.append(np.array(results.expect[3]))

plt.imshow(state0pop, aspect='auto', origin='lower',
           extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
           cmap='viridis')
plt.colorbar(label='Expectation Value')
plt.xlabel('Time')
plt.ylabel('Rabi Rate')
plt.title('state 0')
plt.show()

plt.imshow(state1pop, aspect='auto', origin='lower',
           extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
           cmap='viridis')
plt.colorbar(label='Expectation Value')
plt.xlabel('Time')
plt.ylabel('Rabi Rate')
plt.title('state 1')
plt.show()

plt.imshow(state2pop, aspect='auto', origin='lower',
           extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
           cmap='viridis')
plt.colorbar(label='Expectation Value')
plt.xlabel('Time')
plt.ylabel('Rabi Rate')
plt.title('state 2')
plt.show()

plt.imshow(state3pop, aspect='auto', origin='lower',
           extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
           cmap='viridis')
plt.colorbar(label='Expectation Value')
plt.xlabel('Time')
plt.ylabel('Rabi Rate')
plt.title('state 3')
plt.show()
