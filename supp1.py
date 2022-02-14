import funs as funs
import sys
import ast

"""
PARAMETERS
beta1, beta2, A, time (total number of timestep for num simulation)
"""
beta1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.9449
beta2 = float(sys.argv[2]) if len(sys.argv) > 2 else -.9801
A = ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else [1,-1.1,1,-.84]
time = int(sys.argv[4]) if len(sys.argv) > 4 else 400

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
