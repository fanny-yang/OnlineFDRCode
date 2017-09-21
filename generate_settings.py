# Import Python libraries
import numpy as np
from numpy import sqrt, log, exp, mean, cumsum, sum, zeros, ones, argsort, argmin, argmax, array, maximum, concatenate
from numpy.random import randn, rand
np.set_printoptions(precision = 4)
import os
import scipy.optimize as optim
from scipy.stats import norm
from scipy.stats import bernoulli
import time
from datetime import datetime

from toimport import *


def get_hyp(pi_max, pi1c, num_hyp):
        
    # Read hyp from file
    filename_pre = "H_PM%.2f_PIC%d_NH%d_" % (pi_max, pi1c, num_hyp)
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


def generate_hyp(pi1c, pi_max, max_hyp, samples):
    
    pi_filename = 'PIC.dat'
    
    # ---- Get pi1 progression ---- #
    if pi_max == 0:
        #%%% If want a homogeneous progression of pi1
        f = open('./expsettings/%s' % pi_filename)
        pi_list = f.read().split('\n')

        # Convert to number array
        cache1 = pi_list[pi1c].split(' ')
        pi1_vec = np.array([float(cache1[i]) for i in range(len(cache1))])

        # Determine vector of lengths for constant pieces (with length hyp-step)
        hyp_steps = len(set(pi1_vec)) # number of piecewise constant
        len_per = max_hyp / hyp_steps # lengths of pcw constant
        len_last = max_hyp - hyp_steps*len_per + len_per 
        length_vec = np.concatenate((len_per*np.ones(hyp_steps-1), np.array([len_last])))
    else: 
        #%%% If want constant pi1
        pi1_vec = np.ones(max_hyp)*pi_max

        # Caculate lengths of constant pieces using pi1_vec and max_hyp
        hyp_steps = 1
        length_vec = [max_hyp]
        
    hyp_mat = np.zeros([samples, max_hyp])    
    
    # ---- Sample hypotheses vectors using the pi1 progression ------ #
    for i in range(samples):
        
        Hyp = np.array([])
        for j in range(hyp_steps):        
            Hyp = np.concatenate((Hyp, bernoulli.rvs(pi1_vec[j], size = length_vec[j])))
        
        hyp_mat[i] = Hyp
            
    
    # ----- Save sample hypotheses vectors ----- #
    dirname = './expsettings'
    filename = "H_PM%.2f_PIC%d_NH%d_" % (pi_max, pi1c, max_hyp)
    saveres(dirname, filename, hyp_mat)
    return hyp_mat
    
# def create_hyp(pi1c, pirange, hyprange = np.arange(1000, 10000,1000)):
#     for max_hyp in hyprange:
#         for pi_max in pirange:
#             generate_hyp(pi1c, pi_max, max_hyp, 20)

def get_absvec(abs_style, abs_prob, num_hyp):

    # Read abs vec from file
    filename_pre = "A_S%d_P%.2f_NH%d_" % (abs_style, abs_prob, num_hyp)
    abs_filename = [filename for filename in os.listdir('./expsettings') if filename.startswith(filename_pre)]
    if len(abs_filename) > 0:
        # Just take the first sample
        abs_mat = np.loadtxt('./expsettings/%s' % abs_filename[0])   
    else:
        print("Absvec file doesn't exist, thus generating the file now ...")
        # Generate 100 draws of num_hyp hypotheses with given pi_1 setting
        abs_mat = generate_abs(abs_style, abs_prob, num_hyp,100)

    # Choose some Hypvector could choose a different sample

    abs_vec = abs_mat[0]
    
    return abs_vec    

def generate_absvec(abs_style, abs_prob, num_hyp):
    pi_filename = 'PIC.dat'
    # abs_style is essentially pi1c for hyp
    pi1c = abs_style
    absvec = np.zeros(length)

    # ---- Get pi1 progression ---- #
    if abs_prob == 0:
        #%%% If want a homogeneous progression of pi1
        f = open('./expsettings/%s' % pi_filename)
        pi_list = f.read().split('\n')

        # Convert to number array
        cache1 = pi_list[pi1c].split(' ')
        pi1_vec = np.array([float(cache1[i]) for i in range(len(cache1))])

        # Determine vector of lengths for constant pieces (with length hyp-step)
        abs_steps = len(set(pi1_vec)) # number of piecewise constant
        len_per = max_hyp / abs_steps # lengths of pcw constant
        len_last = max_hyp - abs_steps*len_per + len_per 
        length_vec = np.concatenate((len_per*np.ones(abs_steps-1), np.array([len_last])))
        
    abs_mat = np.zeros([samples, max_hyp])    
    
    # ---- Sample hypotheses vectors using the pi1 progression ------ #
    for i in range(samples):
        
        Abs = np.array([])
        for j in range(hyp_steps):        
            Abs = np.concatenate((Abs, bernoulli.rvs(pi1_vec[j], size = length_vec[j])))
        
        abs_mat[i] = Abs
            
    
    # ----- Save sample hypotheses vectors ----- #
    dirname = './expsettings'
    filename = "A_P%.2f_S%d_NH%d_" % (abs_prob, abs_style, num_hyp)
    saveres(dirname, filename, abs_mat)
    return abs_mat

    # # Random
    # if abs_style == 0:
    #     # Draw Bernoulli
    #     abs_vec = [np.random.binomial(1, abs_prob) for i in range(num_hyp)]
    #  elif abs_style == 1:
    #      # Draw with varying pi1c

def create_pen(penw_style, penw_const, prw_vec,  NUMHYP):
    # Constant weights
    if penw_style == 1:
        penw_vec = np.ones(NUMHYP)*penw_const
    # Linear Decreasing
    elif penw_style == 2:
        penw_vec = np.linspace(penw_const, 1, num = NUMHYP)
    # Linear increasing
    elif penw_style == 3:
        penw_vec = np.linspace(1, penw_const, num = NUMHYP)
    # Correlated
    elif penw_style == 4:
        penw_vec = prw_vec
    return penw_vec

def create_pr(prw_style, prw_const, m_corr, Hypo, NUMHYP):
    prw_vec = Hypo.astype(float)
    diff = m_corr - 1
    zero_indices= np.where(np.array(Hypo == 0))[0]

    # Correlated
    if prw_style == 1:
        # If Hyp=1 set 2
        prw_vec = prw_vec*m_corr
        # In the additive model
        prw_vec[zero_indices] = (1-diff)*np.ones(len(zero_indices))
        # # In the ratio model
        # prw_vec[zero_indices] = np.true_divide(np.ones(len(zero_indices)),m_corr)

        
    # Anti correlated
    elif prw_style == 2:
        # In the additive model
        prw_vec = (1-diff)*np.ones(len(zero_indices)) 
        # # In the ratio model
        # prw_vec = np.true_divide(np.ones(NUMHYP),m_corr)
        prw_vec[zero_indices] = np.ones(len(zero_indices))*m_corr
        
    # Constant 
    elif prw_style == 3:
        prw_vec = np.ones(NUMHYP)*prw_const
    return prw_vec
        
        
        
