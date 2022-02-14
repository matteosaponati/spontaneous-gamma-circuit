import numpy as np
import sys
import ast
import funs as funs

"""
PARAMETERS
beta1, beta2, A, time (total number of timestep for num simulation)
"""
beta1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.9449
beta2 = float(sys.argv[2]) if len(sys.argv) > 2 else -.9801
A = ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else [1,-1.1,1,-.84]
time = int(sys.argv[4]) if len(sys.argv) > 4 else 1000

"""
NUMERICAL SOLUTIONS
find the solution numerically for the AR(2) model
compute affine transformation in the (E,I) space
"""
ar,E,I,eps = funs.num_solution(beta1,beta2,A,time)

'save'
np.save('num_solution_ei',[E,I])
np.save('num_solution_ar',ar)

