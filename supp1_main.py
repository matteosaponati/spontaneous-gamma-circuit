import linear_ei_funs as funs

beta1 = 1.9449
beta2 = -.9801
A = [1,-1.1,1,-.84]
time = 400
N = 100

"""
PLOT NUMERICAL SOLUTIONS and PHASE SPACE TRANSFORMATION
A: chosen affine transformation
beta1, beta2: chosen AR(2) parameters
time: length of the simulation
N: number of repetitions for averaging
savedir: directory for saving
"""

funs.plot_dynamics(A,beta1,beta2,time,N)
