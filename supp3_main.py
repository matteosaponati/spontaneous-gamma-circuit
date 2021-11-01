import linear_ei_funs as funs

beta1 = 1.9449
beta2 = -.9801
A = [1,-1.1,1,-.84]
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