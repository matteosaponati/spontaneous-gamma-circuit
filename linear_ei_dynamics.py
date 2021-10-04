import numpy as np
import matplotlib.pyplot as plt

"""
time: lenght of the simulation
A: transformation matrix
beta1, beta2: parameter of the AR(2) model
"""
def num_solution(time,A,beta1,beta2):
    ar = np.zeros(time)
    eps = np.random.randn(time)
    "numerical solution"
    for t in np.arange(2,time,1):
        ar[t] = beta1*ar[t-1] + beta2*ar[t-2] + eps[t]
    "affine transformation"
    E = A[0]*ar[1:] + A[1]*ar[:-1] 
    I = A[2]*ar[1:] + A[3]*ar[:-1] 
    return ar,E,I,eps

"""
compute weights, equation (20)
A: transformation matrix
beta1, beta2: parameter of the AR(2) model
"""
def compute_weights(A,beta1,beta2):
    weights=[]
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(A[3]*A[0]*beta1 + A[3]*A[1] - A[0]*A[2]*beta2))
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(-A[1]*A[0]*beta1 - A[1]*A[1] + A[0]*A[0]*beta2))
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(A[3]*A[2]*beta1 + A[3]*A[3] - A[2]*A[2]*beta2))
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(-A[1]*A[2]*beta1 - A[1]*A[3] + A[0]*A[2]*beta2))
    return weights, [r'$v_{ee}$',r'$v_{ei}$',r'$v_{ie}$',r'$v_{ii}$']

"""
averaging over many time series
to compute distribution in the phase space
"""
def phase_space_average(A,beta1,beta2,time,N):
    ar_m = np.zeros(time)
    E_m,I_m = np.zeros(time-1), np.zeros(time-1)
    
    for k in range(N):
        ar,E,I,_ = num_solution(time,A,beta1,beta2)
        ar_m, E_m, I_m += ar, E, I
    
    return ar_m/N,E_m/N,I_m/N

"""
FIG 1
A: chosen affine transformation
beta1, beta2: chosen AR(2) parameters
time: length of the simulation
N: number of repetitions for averaging
savedir: directory for saving
"""
def plot_dynamics(A,beta1,beta2,time,N,savedir):
    
    fig = plt.figure(figsize=(12,14), dpi=300)
    "damped oscillations dynamics in time"
    ar,E,I,_ = num_solution(time,A,beta1,beta2)
    plt.subplot(4,2,1)
    plt.xlabel('timesteps')
    plt.ylabel(r'$x_t$')
    plt.plot(ar,color='xkcd:purple',linewidth=1.5)
    plt.subplot(4,2,2)
    plt.xlabel('timesteps')
    plt.plot(E,color='firebrick',linewidth=1.5,label=r'$E_t$')
    plt.plot(I,color='darkblue',linewidth=1.5,label=r'$I_t$')
    plt.legend()
    "phase space"
    ar,E,I = phase_space_average(A,beta1,beta2,10000,N)
    print('done')
    plt.subplot(4,2,5)
    plt.xlabel(r'$x_{t-1}$')
    plt.ylabel(r'$x_t$')
    plt.hist2d(ar[:-1],ar[1:],bins=100,cmap='Purples')
    plt.subplot(4,2,6)
    plt.xlabel(r'$E_t$')
    plt.ylabel(r'$I_t$')
    plt.hist2d(E,I,bins=100,cmap='Purples')
    "critical oscillations dynamics in time"
    ar,E,I,_ = num_solution(time,A,beta1,-1.)
    plt.subplot(4,2,3)
    plt.xlabel('timesteps')
    plt.ylabel(r'$x_t$')
    plt.plot(ar,color='xkcd:purple',linewidth=1.5)
    plt.subplot(4,2,4)
    plt.xlabel('timesteps')
    plt.plot(E,color='firebrick',linewidth=1.5,label=r'$E_t$')
    plt.plot(I,color='darkblue',linewidth=1.5,label=r'$I_t$')
    plt.legend()
    "phase space"
    ar,E,I = phase_space_average(A,beta1,-1.,10000,N)
    print('done')
    plt.subplot(4,2,7)
    plt.xlabel(r'$x_{t-1}$')
    plt.ylabel(r'$x_t$')
    plt.hist2d(ar[:-1],ar[1:],bins=100,cmap='Purples')
    plt.subplot(4,2,8)
    plt.xlabel(r'$E_t$')
    plt.ylabel(r'$I_t$')
    plt.hist2d(E,I,bins=100,cmap='Purples')
    
    fig.tight_layout(rect=[0, 0., 1, 0.98]) 
    plt.savefig(savedir+'/dynamics_A_{}.png'.format(A), format='png', dpi=300)
    plt.close()
    del ar,E,I
    return
