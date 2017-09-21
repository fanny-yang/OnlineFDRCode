import numpy as np

from scipy.stats import norm
from scipy.stats import bernoulli

class rowexp_new_batch:

    def __init__(self, NUMHYP, numdraws, alt_vec, mu0, mu_alt_vec, sigma = 1): 
        self.numhyp = NUMHYP
        self.alt_vec = alt_vec
        self.mu0 = mu0
        self.mu_vec = mu0*np.ones(NUMHYP) + np.multiply(alt_vec, mu_alt_vec)
        self.pvec = np.zeros(NUMHYP)
        self.numdraws = numdraws

        '''
        Function drawing p-values: Mixture of two Gaussians
        '''
    def gauss_two_mix(self, sigma = 1, rndsd = 0):
        np.random.seed(rndsd)
        Z = np.zeros([self.numdraws, self.numhyp])
        
        # Draw Z values
        for i in range(self.numdraws):
            Z[i] = self.mu_vec + np.random.randn(self.numhyp)*sigma # draw gaussian acc. to hypothesis, if sigma are all same  
        Z_av = [np.average([Z[k][i] for k in range(self.numdraws)]) for i in range(self.numhyp)]
        # Compute p-values and save
        self.pvec = [2*norm.cdf(-np.true_divide(abs(z), np.true_divide(sigma, np.sqrt(self.numdraws)))) for z in Z_av]


    def bernoulli_draws(self, rndsd = 0):
        np.random.seed(rndsd)
        # Draw a couple of Bernoullis and compute of successes
        binom_succ = [np.random.binomial(self.numdraws,self.mus[self.alt_vec[i]],1) for i in range(self.numhyp)]
        # Compute p val
        self.pvec = [stats.binom_test(binom_succ, self.numdraws, self.mus[self.alt_vec[i]]) for i in range(self.numhyp)]
