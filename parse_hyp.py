# FDR framework
import numpy as np
from numpy import sqrt, log, exp, mean, cumsum, sum, zeros, ones, argsort, argmin, argmax, array, maximum, concatenate
from numpy.random import randn, rand
np.set_printoptions(precision = 4)

import os
import time
#import ipdb
from datetime import datetime
from generate_hyp import*

def get_hyp(pi_max, pi1c, num_hyp):
        
    # Read hyp from file
    filename_pre = "H_PM%.1f_PIC%d_NH%d_" % (pi_max, pi1c, num_hyp)
    hypo_filename = [filename for filename in os.listdir('./expsettings') if filename.startswith(filename_pre)]
    if len(hypo_filename) > 0:
        # Just take the first sample
        hyp_mat = np.loadtxt('./expsettings/%s' % hypo_filename[0])    
    else:
        print("Hyp file doesn't exist, thus generating the file now ...")
        # Generate 100 draws of num_hyp hypotheses with given pi_1 setting
        hyp_mat = generate_hyp(pi1c, pi_max, num_hyp, 100)

    # Choose some Hypvector could choose a different sample
    Hypo = hyp_mat[0]
    
    return Hypo

