import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
from function import Gaussian_square
##################################################################################################
rabi_rate = np.linspace(0.0, 1, 101)
EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
flat, sig, phi_0, time_step, Delta = 50.0, 34.3, 0, 0.02, 0.01
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

fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
imshows = []
contour_levels = [0.5, 1.0]

X, Y = np.meshgrid(tlist, rabi_rate)

# State 0
im0 = axs[0, 0].imshow(state0pop, aspect='auto', origin='lower',
                    extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
                    cmap='cividis')

axs[0, 0].set_title('state 0')
fig.colorbar(im0, ax=axs[0, 0])

# State 1
im1 = axs[0, 1].imshow(state1pop, aspect='auto', origin='lower',
                    extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
                    cmap='cividis')

c0 = axs[0, 1].contour(X, Y, state0pop, levels=contour_levels, colors='w', linewidths=1.5)
# axs[0, 1].clabel(c0, fmt="%.1f", colors='w')

c2 = axs[0, 1].contour(X, Y, state2pop, levels=contour_levels, colors='yellow', linewidths=1.5)
# axs[0, 1].clabel(c2, fmt="%.1f", colors='yellow')

c3 = axs[0, 1].contour(X, Y, state3pop, levels=[0.01], colors='red', linewidths=1.5)
# axs[0, 1].clabel(c3, fmt="%.1f", colors='red')

axs[0, 1].set_title('state 1')
fig.colorbar(im1, ax=axs[0, 1])

# State 2
im2 = axs[1, 0].imshow(state2pop, aspect='auto', origin='lower',
                    extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
                    cmap='cividis')
axs[1, 0].set_title('state 2')
fig.colorbar(im2, ax=axs[1, 0])

# State 3
im3 = axs[1, 1].imshow(state3pop, aspect='auto', origin='lower',
                    extent=[tlist[0], tlist[-1], rabi_rate[0], rabi_rate[-1]],
                    cmap='cividis')

axs[1, 1].set_title('state 3')
fig.colorbar(im3, ax=axs[1, 1])

for ax in axs.flat:
    ax.set_xlabel('Time')
    ax.set_ylabel('Rabi Rate')

plt.tight_layout()
plt.show()