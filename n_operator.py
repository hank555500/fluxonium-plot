import scqubits as sc
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt

fluxonium = sc.Fluxonium(EJ=5, EC=1.5, EL=0.15, flux=0.0, cutoff=100)
flux_list = np.linspace(-0.5, 0.5, 100)

data = fluxonium.get_matelements_vs_paramvals(operator = 'n_operator' ,param_name = 'flux', param_vals = flux_list, evals_count = 6)

for i in range(5):
    for j in range(5):
        if j >= i:
            plt.plot(flux_list, np.abs(data.matrixelem_table[:, i, j]), label = f'{i},{j}')

plt.xlabel(r'$\Phi_{\mathrm{ext}}/\Phi_0$')
plt.ylabel(r'$\langle n \rangle$')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
###################################################################################################################
#內建製圖
# fig, ax = fluxonium.plot_matelem_vs_paramvals('n_operator', 'flux', flux_list, select_elems=4)
# plt.show()