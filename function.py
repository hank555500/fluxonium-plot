import numpy as np 
import scqubits as scq
import qutip as qt
#################################################################################################
def Gaussian_square(freq, flat, sig, phi_0=0, num_sig=2, time_step=None, norm=True):
    """
    total duration = flat+sig*num_sig*2
    norm: normalize gaussian edge to 0 and 1
    unit: ns
    """
    if time_step == None:
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
def objective_function(params):
    flat_opt, sig_opt, rabi_rate_opt = params
    EJ, EC, EL, phi_eA, N_q = 5.0, 1.5, 0.15, 0.0, 7 
    fluxonium01 = scq.Fluxonium(EC=EC, EL=EL, EJ=EJ, flux=phi_eA, cutoff=150, truncated_dim=N_q)
    phi_0, time_step, Delta =0, 0.02, 0.01
    evals = fluxonium01.eigenvals(evals_count=N_q)
    H = qt.Qobj(np.diag(evals))
    n_opr = qt.Qobj(fluxonium01.matrixelement_table('n_operator', evals_count=N_q))
    psi0 = qt.basis(N_q, 0)
    Rabi_03 =  0.02 * np.pi * 2
    Rabi_23 = Rabi_03 * rabi_rate_opt 
    amp_03, amp_23 = Rabi_03/abs(n_opr[0,3]), Rabi_23/abs(n_opr[2,3])
    e_ops_list = [qt.ket2dm(qt.basis(N_q, i)) for i in range(N_q)]

    tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), flat_opt, sig_opt, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
    tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), flat_opt, sig_opt, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
    
    H_evo = [H, [n_opr, amp_03 * pulse_03], [n_opr, amp_23 * pulse_23]]
    results = qt.sesolve(H_evo, psi0, tlist, e_ops = e_ops_list)

    populations = [results.expect[i][-1] for i in range(len(e_ops_list))]
    
    err = (populations[0] - 0.5) ** 2 + (populations[2] - 0.5) ** 2 + populations[1] ** 2 + populations[3] ** 2
    
    return err
##################################################################################################

