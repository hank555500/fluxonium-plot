import scqubits as sc
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
##################################################################################################################
fluxonium = sc.Fluxonium(EJ = 5, EC = 1.5, EL = 0.15, flux = 0.5, cutoff = 100)
flux_list = np.linspace(-0.5, 0.5, 100)

spectrum_data = fluxonium.get_spectrum_vs_paramvals('flux', flux_list, evals_count=8, subtract_ground = True)
for idx in range(spectrum_data.energy_table.shape[1]):
    plt.plot(flux_list, spectrum_data.energy_table[:, idx], label=f'{idx}')

plt.xlabel(r'$\Phi_{\mathrm{ext}}/\Phi_0$')
plt.ylabel(r'$\mathrm{energy\ [GHz]}$')
plt.title('Fluxonium Energy Spectrum vs Flux')
plt.legend()
plt.show()
# #################################################################################################################
# #內建能階圖
# fluxonium.plot_wavefunction(esys=None, which=range(6), mode='real')
# plt.show()
# #內建製普譜線圖
fig, ax = fluxonium.plot_evals_vs_paramvals('flux', flux_list, evals_count=8, subtract_ground = True)
plt.show()