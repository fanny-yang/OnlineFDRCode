# Bonferroni

import numpy as np
from numpy import sqrt, log, exp, mean, cumsum, sum, zeros, ones, argsort, argmin, argmax, array, maximum
        
class BONF_proc_batch: 

    tmp = range(1, 10000)

    gamma_vec =  np.true_divide(np.log(np.maximum(tmp,np.ones(len(tmp))*2)), np.multiply(tmp, np.exp(np.sqrt(np.log(np.maximum(np.ones(len(tmp)), tmp))))))
    gamma_vec = gamma_vec / np.float(sum(gamma_vec)) 

    def __init__(self, alpha0, numhyp):

        self.alpha0 = alpha0
        self.alpha = np.zeros(numhyp + 1)
        self.alpha[0:2] = [0, self.gamma_vec[0]*self.alpha0]  # vector of alpha_js, first can be ignored
        self.wealth_vec = np.zeros(numhyp + 1)
        self.wealth_vec[0] = self.alpha0
    
    def run_fdr(self, pvec):
              
        numhyp = len(pvec)
        last_rej = 0
        rej = np.zeros(numhyp + 1)
        
        
        for k in range(0, numhyp):

            # Get rejection
            this_alpha = self.alpha[k + 1] # make sure first one doesn't do bullshit
            rej[k + 1] = (pvec[k] < this_alpha)

            self.wealth_vec[k + 1] = self.wealth_vec[k] - this_alpha
            # Calc new alpha
            #next_alpha = self.alpha0/(2*(k+2)**2)
            next_alpha = self.gamma_vec[k + 1]*self.alpha0
            if k < numhyp - 1:
                self.alpha[k + 2] = next_alpha


        rej = rej[1:]
        self.alpha = self.alpha[1:]

        return rej
