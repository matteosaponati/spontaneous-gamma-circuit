## LINEAR E-I SYSTEM

import sys
import numpy as np

## set your main directory
main_dir =''
sys.path.append(main_dir)

import linear_ei_parspace as plt_parspace
import linear_ei_dynamics as plt_dyn

#%%
## PARAMETERS
## beta1, beta2: parameters of the AR(2) model
## A: coefficients of the affine transformation
## a[0] = a_1, a[1] = a_2, a[2] = a_3, a[3] = a_4 
beta1 = 1.9449
beta2 = -.9801
A = [1,-1.1,1,-.84]

#%%
## NUMERICAL SOLUTION FUNCTION
## time: lenght of the simulation
## A: transformation matrix
## beta1, beta2: parameter of the AR(2) model
def num_solution(time,A,beta1,beta2):
    ## initialize variables
    ar = np.zeros(time)
    eps = np.random.randn(time)
    ## numerical solution
    for t in np.arange(2,time,1):
        ar[t] = beta1*ar[t-1] + beta2*ar[t-2] + eps[t]
    ## affine transformation
    E = A[0]*ar[1:] + A[1]*ar[:-1] 
    I = A[2]*ar[1:] + A[3]*ar[:-1] 
    return ar,E,I,eps

## COMPUTE WEIGHTS
## A: transformation matrix
## beta1, beta2: AR(2) parameter
def compute_weights(A,beta1,beta2):
    weights=[]
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(A[3]*A[0]*beta1 + A[3]*A[1] - A[0]*A[2]*beta2))
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(-A[1]*A[0]*beta1 - A[1]*A[1] + A[0]*A[0]*beta2))
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(A[3]*A[2]*beta1 + A[3]*A[3] - A[2]*A[2]*beta2))
    weights.append((1/(A[0]*A[3]-A[1]*A[2]))*(-A[1]*A[2]*beta1 - A[1]*A[3] + A[0]*A[2]*beta2))
    return weights, [r'$v_{ee}$',r'$v_{ei}$',r'$v_{ie}$',r'$v_{ii}$']

#%%
## NUMERICAL SOLUTIONS
## numerically solve the dynamics of the AR(2) model and compute
## the affine transformation in the (E,I) space
## time: length of the simulation

time = 1000
ar,E,I,eps = num_solution(time,A,beta1,beta2)

#%% 
## FIG 1
## PLOT NUMERICAL SOLUTIONS and PHASE SPACE TRANSFORMATION
## A: chosen affine transformation
## beta1, beta2: chosen AR(2) parameters
## time: length of the simulation
## N: number of repetitions for averaging
## savedir: directory for saving

time = 400
N = 10000
plt_dyn.plot_dynamics(A,beta1,beta2,time,N,main_dir)

#%%
## FIG 2
## WEIGHTS IN THE (beta1,beta2) PARAMETER SPACE
## A: chosen affine transformation
## b1_range: minimum and maximum value for the beta1 span [b1_min,b1_max]
## b2_range: minimum and maximum value for the beta2 span [b2_min,b2_max]
## savedir: directory for saving
b1_range, b2_range = [-2,2], [-1,1]
span = 800

plt_parspace.weights_b1b2(A,b1_range,b2_range,span,main_dir)
#plt_parspace.weights_b1b2_proxy(A,main_dir)

 #%%
## FIG 3
## PHASE SHIFT AND E-I BALANCE IN THE (a2,a4) SPACE
## beta1, beta2 : chosen AR(2) parameters
## A: chosen affine transformation
## a2_range: minimum and maximum value for the a2 span [a2_min,a2_max]
## a4_range: minimum and maximum value for the a4 span [a4_min,a4_max]
## savedir: directory for saving
a2_range, a4_range = [-4,4], [-4,4]
span = 800

plt_parspace.balance_a2a4(beta1,beta2,A,a2_range,a4_range,span,main_dir)
plt_parspace.weights_a2a4(beta1,beta2,A,a2_range,a4_range,span,main_dir)
plt_parspace.phase_a2a4(beta1,beta2,A,a2_range,a4_range,span,main_dir)