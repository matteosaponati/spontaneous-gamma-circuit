import numpy as np
import matplotlib.pyplot as plt
import sys

from utils.num_simulations import num_solution

"""
PARAMETERS
beta1, beta2, time, A
"""

A = np.zeros((2,2))
N = 100

beta1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.9449
beta2 = float(sys.argv[2]) if len(sys.argv) > 2 else -.9801
time = int(sys.argv[3]) if len(sys.argv) > 3 else 400
A[0,0] = float(sys.argv[4]) if len(sys.argv) > 4 else 1.
A[0,1] = float(sys.argv[5]) if len(sys.argv) > 5 else -1.1
A[1,0] = float(sys.argv[6]) if len(sys.argv) > 6 else 1.
A[1,1] = float(sys.argv[7]) if len(sys.argv) > 7 else -.84

"""
PLOT NUMERICAL SOLUTIONS and PHASE SPACE TRANSFORMATION
A: chosen affine transformation
beta1, beta2: chosen AR(2) parameters
time: length of the simulation
N: number of repetitions for averaging
savedir: directory for saving
"""

fig = plt.figure(figsize=(12,14), dpi=300)

"damped oscillations dynamics in time"
ar,E,I = num_solution(A,beta1,beta2,time=400)
plt.subplot(4,2,1)
plt.xlabel('timesteps')
plt.ylabel(r'$x_t$')
plt.plot(ar,color='xkcd:purple',linewidth=1.5)
plt.subplot(4,2,2)
plt.xlabel('timesteps')
plt.plot(E,color='firebrick',linewidth=1.5,label=r'$E_t$')
plt.plot(I,color='darkblue',linewidth=1.5,label=r'$I_t$')
plt.legend()

ar,E,I = num_solution(A,beta1,beta2,10000,N)
plt.subplot(4,2,5)
plt.xlabel(r'$x_{t-1}$')
plt.ylabel(r'$x_t$')
plt.hist2d(ar[:-1],ar[1:],bins=100,cmap='Purples')
plt.subplot(4,2,6)
plt.xlabel(r'$E_t$')
plt.ylabel(r'$I_t$')
plt.hist2d(E,I,bins=100,cmap='Purples')

"critical oscillations dynamics in time"
ar,E,I = num_solution(A,beta1,-1.,time)
plt.subplot(4,2,3)
plt.xlabel('timesteps')
plt.ylabel(r'$x_t$')
plt.plot(ar,color='xkcd:purple',linewidth=1.5)
plt.subplot(4,2,4)
plt.xlabel('timesteps')
plt.plot(E,color='firebrick',linewidth=1.5,label=r'$E_t$')
plt.plot(I,color='darkblue',linewidth=1.5,label=r'$I_t$')
plt.legend()

ar,E,I = num_solution(A,beta1,-1.,10000,N)
plt.subplot(4,2,7)
plt.xlabel(r'$x_{t-1}$')
plt.ylabel(r'$x_t$')
plt.hist2d(ar[:-1],ar[1:],bins=100,cmap='Purples')
plt.subplot(4,2,8)
plt.xlabel(r'$E_t$')
plt.ylabel(r'$I_t$')
plt.hist2d(E,I,bins=100,cmap='Purples')

fig.tight_layout(rect=[0, 0., 1, 0.98]) 
plt.savefig('plots/supp1_dynamics_A_{}.png'.format(A.flatten()), format='png', dpi=300)
plt.close()