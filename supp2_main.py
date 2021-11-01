import linear_ei_funs as funs

A = [1,-1.1,1,-.84]
b1_range, b2_range = [-2,2], [-1,1]
span = 800

"""
WEIGHTS IN THE (beta1,beta2) PARAMETER SPACE
A: chosen affine transformation
b1_range: minimum and maximum value for the beta1 span [b1_min,b1_max]
b2_range: minimum and maximum value for the beta2 span [b2_min,b2_max]
savedir: directory for saving
"""
funs.weights_b1b2(A,b1_range,b2_range,span)