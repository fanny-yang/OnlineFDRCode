import numpy as np
#import ipdb

class LORD2_proc_batch: 

    tmp = range(1, 10000)
    # discount factor \gamma_m
    gamma_vec =  np.true_divide(np.log(np.maximum(tmp,np.ones(len(tmp))*2)), np.multiply(tmp, np.exp(np.sqrt(np.log(np.maximum(np.ones(len(tmp)), tmp))))))
    gamma_vec = gamma_vec / np.float(sum(gamma_vec))    

    def __init__(self, alpha0, numhyp, startfac):
        self.alpha0 = alpha0
        self.w0 = startfac*self.alpha0
        self.wealth_vec = np.zeros(numhyp + 1)
        self.wealth_vec[0] = self.w0
        self.alpha = np.zeros(numhyp + 1)
        self.alpha[0:2] = [0, self.gamma_vec[0]*self.w0]  # vector of alpha_js, first can be ignored
        self.b0 = self.alpha0 - self.w0
    
    def run_fdr(self, pvec):

        numhyp = len(pvec)
        # last rejections vector, starting as empty
        last_rej = []
        
        rej = np.zeros(numhyp + 1)

        for k in range(0, numhyp):    

            if self.wealth_vec[k] > 0:
                # Get rejection
                this_alpha = self.alpha[k + 1] # make sure first one doesn't do bullshit
                rej[k + 1] = (pvec[k] < this_alpha)

                # Calc wealth
                if (rej[k + 1] == 1):
                    last_rej = np.append(last_rej, k + 1).astype(int)
                # Update wealth
                wealth = self.wealth_vec[k] - this_alpha + rej[k + 1]*self.b0
                self.wealth_vec[k + 1] = wealth

                # Calc new alpha
                if len(last_rej) > 0:

                    sum_gam = sum(self.gamma_vec[((k+1)*np.ones(len(last_rej),dtype=int) - last_rej)])
                    # sum_gam = self.gamma_vec[[int(i) for i in ((k+1)*np.ones(len(last_rej)) - last_rej)]][0]
                    
                    next_alpha = self.gamma_vec[k+1]*self.wealth_vec[0] + self.b0*sum_gam

                else:
                    sum_gam = 0
                    next_alpha = self.gamma_vec[k+1]*self.wealth_vec[0]
                
                if k < numhyp - 1:
                    self.alpha[k + 2] = next_alpha    
            else: 
                break
        rej = rej[1:]
        self.alpha = self.alpha[1:]
        return rej
