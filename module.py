"""Utility to evaluate the expected number of new events"""
import math
import numpy as np
def smoothedGT(Phi,n,t,sd_samples):
    '''Compute the smoothed Good-Toulmin estimator.
    Base on arXiv:1511.07428.
    Input list:
        Phi: histogram of histogram (or list of prevalences)
        n: sample size
        t: fold-changes of the new sample size wrt n
        sd_samples: number of samples from the smoothing distro
    '''
    U_list = []
    if t <= 1 :
        coeff = []
        for i in range(len(Phi)):
            coeff.append(-(-t)**(i+1))
        U_GT = np.dot(coeff,Phi) # the expected number of new species with a new sampling of the same size (t=1)
        U_list.append(U_GT)
    if t>1:
        k = math.ceil( 0.5*math.log(1.0*n*t**2/(t-1),3) )
        q = 2.0/(t+2)
        L_list = np.random.binomial(k,q,size=(sd_samples,1))
        for L in L_list:
            coeff = []
            for i in range(L[0]):
                coeff.append(-(-t)**(i+1))
            if len(coeff)-len(Phi) > 0:
                Phi = np.lib.pad(Phi,(0,len(coeff)-len(Phi)),'constant',constant_values=0)
            U_GT = np.dot(coeff,Phi[:L[0]]) # the expected number of new species with a new sampling of the same size (t=1)
            U_list.append(U_GT)
    return np.mean(U_list)
