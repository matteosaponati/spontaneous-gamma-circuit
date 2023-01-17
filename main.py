import numpy as np
import matplotlib.pyplot as plt
import sys

from utils.num_simulations import num_solution

"""
PARAMETERS

beta1, beta2, time, A
"""
beta1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.9449
beta2 = float(sys.argv[2]) if len(sys.argv) > 2 else -.9801
time = int(sys.argv[3]) if len(sys.argv) > 3 else 400
A[0,0] = float(sys.argv[4]) if len(sys.argv) > 4 else 1.
A[0,1] = float(sys.argv[5]) if len(sys.argv) > 5 else -1.1
A[1,0] = float(sys.argv[6]) if len(sys.argv) > 6 else 1.
A[1,1] = float(sys.argv[7]) if len(sys.argv) > 7 else -.84

"""
NUMERICAL SOLUTIONS

find the solution numerically for the AR(2) model
compute affine transformation in the (E,I) space
"""
ar,E,I = funs.num_solution(beta1,beta2,A,time)

'save'
np.save('results/num_solution_ei',[E,I])
np.save('results/num_solution_ar',ar)

