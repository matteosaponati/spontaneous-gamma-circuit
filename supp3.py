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

a2_range, a4_range = [-4,4], [-4,4]
span = 800

"""
PHASE SHIFT AND E-I BALANCE IN THE (a2,a4) SPACE
beta1, beta2 : chosen AR(2) parameters
A: chosen affine transformation
a2_range: minimum and maximum value for the a2 span [a2_min,a2_max]
a4_range: minimum and maximum value for the a4 span [a4_min,a4_max]
savedir: directory for saving
"""
print('EI balance in the (a2,a4) space')
funs.balance_a2a4(beta1,beta2,A,a2_range,a4_range,span)
print('weights in the (a2,a4) space')
funs.weights_a2a4(beta1,beta2,A,a2_range,a4_range,span)
print('phase shift in the (a2,a4) space')
funs.phase_a2a4(beta1,beta2,A,a2_range,a4_range,span)