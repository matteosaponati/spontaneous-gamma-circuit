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

from utils.num_simulations import compute_weights

A = np.zeros((2,2))

A[0,0] = float(sys.argv[1]) if len(sys.argv) > 1 else 1.
A[0,1] = float(sys.argv[2]) if len(sys.argv) > 2 else -1.1
A[1,0] = float(sys.argv[3]) if len(sys.argv) > 3 else 1.
A[1,1] = float(sys.argv[4]) if len(sys.argv) > 4 else -.84

'define parameter sweep list'
b1_range, b2_range = [-2,2], [-1,1]
span = 800
beta1 = np.linspace(b1_range[0],b1_range[1],span)
beta2 = np.linspace(b2_range[1],b2_range[0],span)

V_tot = np.zeros((2,2,span,span))
for x,y in product(range(span),range(span)):
    V_tot[:,:,x,y] = compute_weights(A, beta1[x], beta2[y])

'-----------------------'
'plot'
fig = plt.figure(figsize=(10,8), dpi=300)

V_plots = [V_tot[0,0,:,:],V_tot[0,1,:,:],V_tot[1,0,:,:],V_tot[1,1,:,:]]

for k in range(4):

    plt.subplot(2,2,k+1)
    plt.ylabel(r'$\beta_2$')
    plt.xlabel(r'$\beta_1$')
    
    plt.xticks(np.linspace(0,span,6),np.round(np.linspace(b1_range[0],b1_range[1],6),2).tolist())
    plt.yticks(np.linspace(0,span,6),np.round(np.linspace(b2_range[1],b2_range[0],6),2).tolist())
    
    plt.axhline(y=span/2,linewidth=.5, linestyle='dashed',color='xkcd:light grey')
    plt.axvline(x=span/2,linewidth=.5, linestyle='dashed',color='xkcd:light grey')
    
    plt.xlim(0,span)
    plt.ylim(span,0)
    
    plt.imshow(V_plots[k].T,aspect='auto',cmap='PRGn',norm=MidpointNormalize(midpoint=0))
    plt.colorbar()
    
    'show region of parameter space based on the dynamical proprieties of the AR(2) model (see Hamilton 1994)'
    plt.plot(np.arange(0,span,1),(2*np.arange(-span/2,span/2,1)),linewidth=.5,color='k')
    plt.plot(np.arange(0,span,1),-(2*np.arange(-span/2,span/2,1)),linewidth=.5,color='k')
    plt.plot(np.arange(0,span,1),span/2+(1/(span/2))*((np.arange(0,span,1)-span/2)**2),linewidth=.5,color='k') 
    
fig.tight_layout(rect=[0, 0.01, 1, 0.98]) 
plt.savefig('plots/supp_1_weights_A_{}.png'.format(A), format='png', dpi=300)
plt.close()   