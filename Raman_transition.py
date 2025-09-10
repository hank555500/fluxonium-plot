import scqubits as sc
import qutip as qt
import numpy as np 
import matplotlib.pyplot as plt

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

flat, sig, phi_0, time_step = 200, 34.3, 0, 0.02
Rabi_0 = 0.010 * 2*np.pi
Delta  = 0.01414
tlist, pulse_03 = Gaussian_square((evals[3]-evals[0]-Delta)/(2*np.pi), 
                                  flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)
tlist, pulse_23 = Gaussian_square((evals[3]-evals[2]-Delta)/(2*np.pi), 
                                  flat, sig, phi_0=phi_0, num_sig=2, time_step=time_step, norm=True)