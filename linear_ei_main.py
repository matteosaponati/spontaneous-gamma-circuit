"LINEAR E-I SYSTEM"
import sys
main_dir = ''
sys.path.append(main_dir)

import linear_ei_parspace as parspace
import linear_ei_dynamics as dyn

"""
NUMERICAL SOLUTIONS
find the solution numerically for the AR(2) model
compute affine transformation in the (E,I) space
"""
beta1 = 1.9449
beta2 = -.9801
A = [1,-1.1,1,-.84]
time = 1000
ar,E,I,eps = dyn.num_solution(time,A,beta1,beta2)

"------------------------"

"""
Supplementary Material: FIG 1 
PLOT NUMERICAL SOLUTIONS and PHASE SPACE TRANSFORMATION
A: chosen affine transformation
beta1, beta2: chosen AR(2) parameters
time: length of the simulation
N: number of repetitions for averaging
savedir: directory for saving
"""
time = 400
N = 100
dyn.plot_dynamics(A,beta1,beta2,time,N,main_dir)

"""
Supplementary Material FIG 2
WEIGHTS IN THE (beta1,beta2) PARAMETER SPACE
A: chosen affine transformation
b1_range: minimum and maximum value for the beta1 span [b1_min,b1_max]
b2_range: minimum and maximum value for the beta2 span [b2_min,b2_max]
savedir: directory for saving
"""
b1_range, b2_range = [-2,2], [-1,1]
span = 800
parspace.weights_b1b2(A,b1_range,b2_range,span,main_dir)

"""
Supplementary Material FIG 3
PHASE SHIFT AND E-I BALANCE IN THE (a2,a4) SPACE
beta1, beta2 : chosen AR(2) parameters
A: chosen affine transformation
a2_range: minimum and maximum value for the a2 span [a2_min,a2_max]
a4_range: minimum and maximum value for the a4 span [a4_min,a4_max]
savedir: directory for saving
"""
a2_range, a4_range = [-4,4], [-4,4]
span = 800
parspace.balance_a2a4(beta1,beta2,A,a2_range,a4_range,span,main_dir)
parspace.weights_a2a4(beta1,beta2,A,a2_range,a4_range,span,main_dir)
parspace.phase_a2a4(beta1,beta2,A,a2_range,a4_range,span,main_dir)
