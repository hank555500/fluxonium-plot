import scqubits as sc
import qutip as qu
import numpy as np 
import matplotlib.pyplot as plt

fluxonium = sc.Fluxonium(EJ = 5, EC = 1.5, EL = 0.150, flux = 0.0, cutoff = 100)
flux_list = np.linspace(-0.5, 0.5, 101)
num_of_eigst = 4

eigenergy, _ = fluxonium.eigensys(evals_count = num_of_eigst) #eigenergy will return 6 number of energy(GHz)1-D array.

# n_op_data = fluxonium.get_matelements_vs_paramvals(operator = 'n_operator' ,param_name = 'flux', param_vals = flux_list, evals_count = num_of_eigst)

for i in range(num_of_eigst):
    for j in range(num_of_eigst):
        if j > i and (j-i)%2 == 1 and j-i<=3:
            fij = np.abs(eigenergy[j]-eigenergy[i])
            plt.vlines(fij, 0, np.abs(fluxonium.matrixelement_table("n_operator", evals_count = num_of_eigst)[i, j]))
            plt.text(fij, np.abs(fluxonium.matrixelement_table("n_operator", evals_count = num_of_eigst)[i, j]), f'{i},{j}')

plt.xlabel(r'$\mathrm{transition\ frequency\ (GHz)}$')
plt.ylabel(r'$\langle n \rangle$')
plt.tight_layout()
plt.savefig(f'transition_frequency({num_of_eigst}esys).png')
# plt.show()

# print(fluxonium.hamiltonian(energy_esys=True))