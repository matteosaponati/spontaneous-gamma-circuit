import numpy as np
import sys
from itertools import product

import matplotlib.pyplot as plt
import matplotlib.colors as colors
class MidpointNormalize(colors.Normalize):
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)
	def __call__(self, value, clip=None):
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
from utils.num_simulations import num_solution, compute_weights

"""
PARAMETERS
beta1, beta2, A, time (total number of timestep for num simulation)
"""

A = np.zeros((2,2))

beta1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.9449
beta2 = float(sys.argv[2]) if len(sys.argv) > 2 else -.9801
time = int(sys.argv[3]) if len(sys.argv) > 3 else 400
A[0,0] = float(sys.argv[4]) if len(sys.argv) > 4 else 1.
A[0,1] = float(sys.argv[5]) if len(sys.argv) > 5 else -1.1
A[1,0] = float(sys.argv[6]) if len(sys.argv) > 6 else 1.
A[1,1] = float(sys.argv[7]) if len(sys.argv) > 7 else -.84

a2_range, a4_range = [-4,4], [-4,4]
span = 100

'---------------------------------------------------------------------------------'
'plots'

'Weights for different combinations of (a2,a4) values'

a2 = np.linspace(a2_range[0],a2_range[1],span)
a4 = np.linspace(a4_range[1],a4_range[0],span)

V_tot = np.zeros((span,span,2,2))
for x,y in product(range(span),range(span)):
    A[0,1], A[1,1] = a2[x], a4[y]
    V_tot[x,y,:,:] = compute_weights(A, beta1, beta2)
    
V_plots = [V_tot[:,:,0,0],V_tot[:,:,0,1],V_tot[:,:,1,0],V_tot[:,:,1,1]]
for k in range(4): V_plots[k][abs(V_plots[k]) >1e10]=np.nan

fig = plt.figure(figsize=(10,8), dpi=300)
weights_name = [r'$v_{ee}$',r'$v_{ei}$',r'$v_{ie}$',r'$v_{ii}$']
for k in range(4):
    plt.subplot(2,2,k+1)
    plt.gca().set_title(weights_name[k])
    plt.ylabel(r'$a_4$')
    plt.xlabel(r'$a_2$')
    plt.xticks(np.linspace(0,span,6),np.round(np.linspace(a2_range[0],a2_range[1],6),2).tolist())
    plt.yticks(np.linspace(0,span,6),np.round(np.linspace(a4_range[1],a4_range[0],6),2).tolist())
    
    plt.axhline(y=(np.abs(a4)).argmin(),linewidth=.5, linestyle='dashed',color='k')
    plt.axvline(x=(np.abs(a2)).argmin(),linewidth=.5, linestyle='dashed',color='k')
    plt.xlim(0,span)
    plt.ylim(span,0)
    
    plt.contour(np.where(abs(V_plots[k].T)<1,1,0),colors='xkcd:navy blue',linewidths=.5,alpha=.7)
    plt.imshow(V_plots[k].T,aspect='auto',cmap='PRGn',norm=MidpointNormalize(midpoint=0), vmin=-20, vmax=20,alpha=.9)
    plt.colorbar()
fig.tight_layout(rect=[0, 0.01, 1, 0.96]) 
plt.savefig('plots/supp_3_weights_b1_{}_b2_{}_A_{}.png'.format(beta1,beta2,A), format='png', dpi=300)
plt.close()

'E-I balance for different combinations of (a2,a4) values'

"""
create_mask_weights: region of coupling weights smaller than 1
create_mask_det: region of det = +/-1
"""
def create_mask_weights(weights):
    mask = np.zeros_like(weights[0].T)
    for idx in np.ndenumerate(weights[0]):
        if  abs(weights[1].T[idx[0]]) <1 and abs(weights[2].T[idx[0]]) <1 : 
            mask[idx[0]] = 1         
    return mask
def create_mask_det(A):
    a2_span = np.tile(np.linspace(a2_range[0],a2_range[1],span),(span,1))
    a4_span = np.tile(np.linspace(a4_range[1],a4_range[0],span),(span,1)).transpose()
    det = (A[0,0]*a4_span - a2_span*A[1,0])
    return np.where(abs(det)<1,1,0)

'Weights for different combinations of (a2,a4) values'
fig = plt.figure(figsize=(6,6), dpi=300)
plt.ylabel(r'$a_4$')
plt.xlabel(r'$a_2$')
plt.xticks(np.linspace(0,span,6),np.round(np.linspace(a2_range[0],a2_range[1],6),2).tolist())
plt.yticks(np.linspace(0,span,6),np.round(np.linspace(a4_range[1],a4_range[0],6),2).tolist())

plt.axhline(y=(np.abs(a4)).argmin(),linewidth=.5, linestyle='dashed',color='k')
plt.axvline(x=(np.abs(a2)).argmin(),linewidth=.5, linestyle='dashed',color='k')
plt.xlim(0,span)
plt.ylim(span,0)

plt.contour(create_mask_weights(V_plots),colors='xkcd:navy blue',linewidths=.5,alpha=.7)
plt.contour(create_mask_det(A),colors='xkcd:mauve',linewidths=.5,alpha=.7)
plt.imshow(V_plots[2].T+V_plots[1].T,aspect='auto',cmap='PRGn',norm=MidpointNormalize(midpoint=0),alpha=.9)
plt.colorbar()

fig.tight_layout(rect=[0, 0.01, 1, 0.98]) 
plt.savefig('plots/supp_3_EI_balance_b1_{}_b2_{}_A_{}.png'.format(beta1,beta2,A), format='png', dpi=300)
plt.close()

"""
Phase shift for different combinations of (a2,a4) values

for every (a_2,a_4) simulate the dynamics of the E-I system and compute the mean angle between I and the AR(2) model, the mean angle
between E and the AR(2) model and then compute the difference between the two estimated angles (default time = 1000)
"""

a2_span = np.tile(np.linspace(a2_range[0],a2_range[1],span),(span,1))
a4_span = np.tile(np.linspace(a4_range[1],a4_range[0],span),(span,1)).transpose()

phase_shift = np.zeros_like(V_plots[0].T)
for idx in np.ndenumerate(a2_span):
    A[0,1],A[1,1] = a2_span[idx[0]],a4_span[idx[0]]
    ar,E,I = num_solution(A,beta1,beta2,1000)
    phase_shift[idx[0]] = np.mean(np.arctan2(I[5:],ar[5:-1])) - np.mean(np.arctan2(E[5:],ar[5:-1]))

fig = plt.figure(figsize=(6,6), dpi=300)
plt.ylabel(r'$a_4$')
plt.xlabel(r'$a_2$')
plt.xticks(np.linspace(0,span,6),np.round(np.linspace(a2_range[0],a2_range[1],6),2).tolist())
plt.yticks(np.linspace(0,span,6),np.round(np.linspace(a4_range[1],a4_range[0],6),2).tolist())

plt.axhline(y=span/2-(np.abs(a4_span[:,0])).argmin(),linewidth=.5, linestyle='dashed',color='xkcd:light grey')
plt.axvline(x=(np.abs(a2_span[0,:])).argmin(),linewidth=.5, linestyle='dashed',color='xkcd:light grey')
plt.xlim(0,span)
plt.ylim(span,0)

plt.contour(create_mask_weights(V_plots),colors='xkcd:navy blue',linewidths=.5,alpha=.7)
plt.contour(create_mask_det(A),colors='xkcd:mauve',linewidths=.5,alpha=.7)
plt.imshow(phase_shift,aspect='auto',cmap='PRGn',norm=MidpointNormalize(midpoint=0),alpha=.9)
plt.plot(a2_span[0,:],a4_span[:,0],linewidth=1,color='xkcd:light grey')
plt.colorbar()

fig.tight_layout(rect=[0, 0.01, 1, 0.98]) 
plt.savefig('plots/supp_3_phase_shift_b1_{}_b2_{}_A_{}.png'.format(beta1,beta2,A), format='png', dpi=300)
plt.close()
