import numpy as np

'------------------------'
def num_solution(A,beta1,beta2,time,N=1):
    """
    time: lenght of the simulation
    A: transformation matrix
    beta1, beta2: parameter of the AR(2) model
    N: number of simulations
    """
    
    ar = np.zeros((time,N))
    eps = np.random.randn(time,N)
    
    for t in np.arange(2,time,1):
        ar[t,:] = beta1*ar[t-1,:] + beta2*ar[t-2,:] + eps[t,:]
        
    [E,I] = np.einsum('ij,jkl->ikl',A,np.stack([ar[1:,:],ar[:-1,:]],axis=0))
    
    return ar.mean(axis=1),E.mean(axis=1),I.mean(axis=1)

'------------------------'
def compute_weights(A,beta1,beta2):
    """
    compute weights, (equation 17)
    A: transformation matrix
    beta1, beta2: parameter of the AR(2) model
    """
    
    B = np.array([[beta1,beta2],[1,0]])
    if np.linalg.det(A) != 0:
        V = A @ B @ np.linalg.inv(A)
    else: V = np.zeros((2,2))
    
    return V

'------------------------'
