import numpy as np
import linear_ei_funs as funs

beta1 = 1.9449
beta2 = -.9801
A = [1,-1.1,1,-.84]
time = 1000

"""
NUMERICAL SOLUTIONS
find the solution numerically for the AR(2) model
compute affine transformation in the (E,I) space
"""
ar,E,I,eps = funs.num_solution(time,A,beta1,beta2)

'save'
np.save('num_solution_ei',[E,I])
np.save('num_solution_ar',[ar])

