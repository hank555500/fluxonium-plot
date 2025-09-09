import scqubits as sc
import qutip as qu
import numpy as np 
import matplotlib.pyplot as plt

fluxonium = sc.Fluxonium(EJ = 5, EC = 1.5, EL = 0.15, flux = 0.5, cutoff = 100)
flux_list = np.linspace(-0.5, 0.5, 100)
num_of_eigst = 6

eigenergy = fluxonium.eigensys(evals_count = num_of_eigst)
n_op = fluxonium.get_matelements_vs_paramvals(operator = 'n_operator' ,param_name = 'flux', param_vals = flux_list, evals_count = num_of_eigst)

