import numpy as np
import matplotlib.pyplot as plt
from qutip import *

omega = 1.0 * 2 * np.pi                          # Rabi頻率
tlist = np.linspace(0, 20, 1000)                 # 時間array
psi0 = basis(2, 0)                               # 初始狀態
H = 0.5 * omega * sigmax()                       # Hamiltonian
gamma = 0.001                                    # 自發衰減率
gamma_phi = 0.2                                  # 去相干率
##########################################################################################
#製圖
# fig, axes = plt.subplots(3, 3, figsize=(13, 10), sharex=True, sharey=True)
# ls = np.logspace(-6, 2, 9)                       #不同數量級
ls = np.linspace(0.25, 2, 9)                        #不同數值(0-0.8)

# for i, gamma in enumerate(ls):
#        row = i // 3
#        col = i % 3
#        ax = axes[row, col]
#        c_op = [np.sqrt(gamma) * sigmam(),                 # 自發衰減算符
#               np.sqrt(gamma_phi) * sigmaz()]              # 退相干算符
#        result = mesolve(H, psi0, tlist, c_op, [])         # 求解薛丁格方城
#        z = expect(sigmaz(), result.states)
#        ax.plot(tlist, z, 'k-' )
#        ax.set_title(f"$\\gamma$ = {gamma:.1f}")
#        ax.set_xlabel('time')
#        ax.set_ylabel('z')
#        ax.grid(True)

# plt.tight_layout()
# plt.show()

##########################################################################################
# c_op = [np.sqrt(gamma) * sigmam(),                 # 自發衰減算符
#        np.sqrt(gamma_phi) * sigmaz()]              # 退相干算符
# result = mesolve(H, psi0, tlist, c_op, [])         # 求解薛丁格方城
# # Bloch球坐標
# x = expect(sigmax(), result.states)
# y = expect(sigmay(), result.states)
# z = expect(sigmaz(), result.states)
# #show出球球
# b = Bloch()
# b.add_points([x, y, z], meth='l')
# b.show()
##########################################################################################
# # 繪圖
# plt.plot(tlist, z)
# plt.xlabel('time')
# plt.ylabel(r'$\langle \sigma_z \rangle$')
# plt.title('Rabi os')
# plt.grid(True)
# plt.show()
##########################################################################################
fig, axes = plt.subplots(3, 3, figsize=(13, 10),  subplot_kw={'projection': '3d'}, sharex=True, sharey=True)

for idx, gamma_phi in enumerate(ls):
    
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]
    c_op = [np.sqrt(gamma) * sigmam(),                 # 自發衰減算符
           np.sqrt(gamma_phi) * sigmaz()]              # 退相干算符
    result = mesolve(H, psi0, tlist, c_op, [])         # 求解薛丁格方城

    x = expect(sigmax(), result.states)
    y = expect(sigmay(), result.states)
    z = expect(sigmaz(), result.states)
    
    b = Bloch(axes=ax)
    b.add_points([x, y, z], meth='l')
    b.render() 
    ax.set_title(f"$\\gamma_\\phi$ = {gamma_phi:.6f}")


plt.tight_layout()
plt.show()