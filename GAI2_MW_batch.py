# This is  LORD GAI++ including memory and weight

import numpy as np
from numpy import sqrt, log, exp, mean, cumsum, sum, zeros, ones, argsort, argmin, argmax, array, maximum

class GAI2_MW_proc_batch: 

    tmp = range(1, 10000)

    # discount factor \gamma_m
    gamma_vec =  np.true_divide(log(np.maximum(tmp, ones(len(tmp))*2)), np.multiply(tmp, exp(sqrt(log(np.maximum(ones(len(tmp)), tmp))))))
    gamma_vec = gamma_vec / np.float(sum(gamma_vec))    

    def __init__(self, alpha0, numhyp, startfac, pr_w, pen_w, mempar):
        pen_min = min(pen_w)
        self.alpha0 = alpha0
        self.w0 = min(pen_min, startfac)*self.alpha0
        self.wealth_vec = np.zeros(numhyp + 1)
        self.wealth_vec[0] = self.w0
        self.alpha = np.zeros(numhyp + 1)
        self.alpha[0:2] = [0, self.gamma_vec[0]*self.w0]  # vector of alpha_js
        self.last_rej = 0 # save last time of rejection
        self.mempar = mempar
        self.pr_w = pr_w
        self.pen_w = pen_w


    # Can add different thresh function
    def thresh_func(self, x):
        return x

    def run_fdr(self, pvec):

        numhyp = len(pvec)
        last_rej = 0
        first = 0 # Whether we are past the first rejection
        flag = 1
        rej = np.zeros(numhyp+1)
        phi = np.zeros(numhyp+1)
        psi = np.zeros(numhyp+1)
        pr_w = self.pr_w
        pen_w = self.pen_w
        # Have to change that since prior weights need to be adjusted 
        r = [np.true_divide(self.thresh_func(pen_w[i]), pr_w[i]) for i in range(numhyp)]

        #penw, prw, pvec without k+1 (just shifted for purpose of having a wealth[0] (rest starts at 1)
        # here t=k+1 (t in paper)
        last_rej = []
        psi_rej = []

        for k in range(0, numhyp):                

            if self.wealth_vec[k] > 0:
                this_alpha = self.alpha[k+1] # make sure first one doesn't do bullshit

                # Calc b, phi
                b_k = self.alpha0 - (1-first)*np.true_divide(self.w0, pen_w[k])
                phi[k + 1] = min(this_alpha, self.mempar*self.wealth_vec[k] + (1-self.mempar)*flag*self.w0)
                # Adjust prior weight depending on penw, phi, b - calc prw, r
                max_weight = (phi[k+1]*self.thresh_func(pen_w[k]))/((1-b_k)*this_alpha)
                pr_w[k] = min(pr_w[k], max_weight)
                r[k] = np.true_divide(self.thresh_func(pen_w[k]), pr_w[k])

                # Calc psi
                # The max is to get rid of numerical issues when computing max_weight
                psi[k + 1] = max(min(phi[k + 1] + pen_w[k]*b_k, np.true_divide(phi[k + 1], this_alpha)*r[k] - pen_w[k] + pen_w[k]*b_k),0)
   
            
                # Rejection decision
                rej[k + 1] = (pvec[k] < np.true_divide(this_alpha,r[k]))

                if (rej[k + 1] == 1):
                    if (first == 0):
                        first = 1
                    last_rej = np.append(last_rej, k + 1).astype(int)
                    psi_rej = np.append(psi_rej, psi[k + 1]).astype(float)
                
                # Update wealth
                wealth = self.mempar*self.wealth_vec[k] + (1-self.mempar)*flag*self.w0 - phi[k + 1] + rej[k + 1]*psi[k + 1]

                # Calc new alpha
                
                if len(last_rej) > 0:
                    # first_gam = self.mempar**(k+1-last_rej[0])*self.gamma_vec[k + 1 - last_rej[0]]
                    #t_taoj = ((k+1)*np.ones(len(last_rej[0:-1]),dtype=int) - last_rej[0:-1])
                    # sum_gam = sum(np.multiply(self.mempar**t_taoj, self.gamma_vec[t_taoj])) - first_gam + self.gamma_vec[k+1 - last_rej[-1]]
                    # next_alpha = self.gamma_vec[k+1]*self.wealth_vec[0] + (self.alpha0 - self.w0/float(pen_w[last_rej[0]]))*first_gam + self.alpha0*sum_gam

                    #gam_vec = np.append(np.multiply(self.mempar**t_taoj, self.gamma_vec[t_taoj]), self.gamma_vec[k + 1 - last_rej[-1]])

                    t_taoj = ((k+1)*np.ones(len(last_rej),dtype=int) - last_rej)
                    gam_vec = np.multiply(self.mempar**t_taoj, self.gamma_vec[t_taoj])

                    next_alpha = self.gamma_vec[k + 1]*self.wealth_vec[0] + np.dot(psi_rej, gam_vec)

                else:
                    sum_gam = 0
                    next_alpha = self.gamma_vec[k+1]*self.wealth_vec[0]

                # next_alpha = self.mempar**(k+1-last_rej)*self.gamma_vec[k+1 - last_rej]*self.wealth_vec[last_rej]
                #next_alpha = self.gamma_vec[k+1 - last_rej]*self.wealth_vec[last_rej]
            else:
                break

            self.wealth_vec[k + 1] = wealth
            if k < numhyp-1:
                self.alpha[k+2] = next_alpha
            
            # After past the first reject, set flag to 0
            if (first == 1):
                flag = 0

        # Cut off the first zero
        rej = rej[1:]
        self.alpha = self.alpha[1:]

        return rej
        
            
    
    
