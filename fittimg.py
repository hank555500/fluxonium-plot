import scqubits as scq
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt
from function import Gaussian_square
from scipy.optimize import minimize
#################################################################################################
EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 4 
fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)



def objective_function(params):
    flat_opt, sig_opt, rabi_rate_opt = params

    phi_0, time_step, Delta =0, 0.02, 0.01
    evals = fluxonium01.eigenvals(evals_count=N_q)
    H = qt.Qobj(np.diag(evals))
    n_opr = qt.Qobj(fluxonium01.matrixelement_table('n_operator', evals_count=N_q))
    psi0 = qt.basis(N_q, 0)
    Rabi_03 =  0.02 * np.pi * 2
    Rabi_23 = Rabi_03 * rabi_rate_opt 
    amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
    e_ops_list = [qt.ket2dm(qt.basis(N_q, 0)), qt.ket2dm(qt.basis(N_q, 1)), qt.ket2dm(qt.basis(N_q, 2)), qt.ket2dm(qt.basis(N_q, 3))]

    tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), flat_opt, sig_opt, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
    tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), flat_opt, sig_opt, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
    
    H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
    results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)

    populations = [results.expect[i][-1] for i in range(len(e_ops_list))]
    
    err = (populations[0] - 0.5) ** 2 + (populations[2] - 0.5) ** 2 + populations[1] ** 2 + populations[3] ** 2
    
    return err


initial_guess = [5.0, 10, 0.421]      # [flat, sig, rabi_rate]
bounds = [(0, 20), (0, 20), (0, 1)]   # 可根據問題調整

result = minimize(objective_function, initial_guess, method='Nelder-Mead', bounds = bounds,
                  options={'maxiter': 300, 'disp': True})

optimal_flat, optimal_sig, optimal_rabi_rate = result.x
print(f"最佳參數: flat = {optimal_flat:.4f}, sig = {optimal_sig:.4f}, rabi_rate = {optimal_rabi_rate:.4f}")
print(f"均方誤差: {result.fun:.8f}")
