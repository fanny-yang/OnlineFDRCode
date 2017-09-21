# Import Python libraries
import numpy as np
np.set_printoptions(precision = 4)

import os
#import matplotlib.pyplot as plt
import scipy.optimize as optim
from scipy.stats import norm
from scipy.stats import bernoulli

# Import FDR procedures
from AlphaInvest_batch import*
from AlphaInvestPlus_batch import*
from LORD2_batch import*
from GAIPlus2_batch import*
#from LORD_batch import*
#from GAIPlus_batch import*
#from GAI_MW_batch import*
from Bonf_batch import*
from GAI2_MW_batch import*
from GAI2_MW_AD_batch import*

# Import utilities
from rowexp_new_batch import*
from toimport import*
from settings_util import*
  
################ Running entire framework  ####################

def run_single(NUMRUN, NUMHYP, NUMDRAWS, mu_gap, penw_style,  penw_const, prw_style,  prw_const, m_corr, pi1c, pi_max, mempar, alpha0, mod_choice, FDR, sigma = 1, verbose = False, TimeString = False, rndseed = 0, startfac = 0.1, abs_style = 0, abs_eps = 0.1, abs_prob = 0, wealth_eps = 0.05): 
# m_corr needs to be an array
# model choice 1: two Gaussian mixture; 2: Bernoulli
# Set pi_max = 0 if you want to use pi1c, set pi1c = 0 if wanna use pi_max
# Penw_style (penalty weight): 1: Constant, 2: Linearly decreasing, 3: Linearly increasing, 4: Correlated
# Prw_style (prior weight): 1: Correlated with m_corr, 2: Anti correlated with m_corr, 3: Constant; 4: random
# m_corr needs to be > 1
    
    if rndseed == 0:
        TimeString = True

    ##------------- Setting hypotheses, penalty and prior weights -------------## 
    if TimeString:
        time_str = datetime.today().strftime("%m%d%y_%H%M")
    else:
        time_str = '0'

    ##### Parse Hypo vector 
    if pi_max > 0:
        pi1c = 0
    elif pi1c > 0:
        pi_max = 0 # Just to make sure pi1c is being used when desired
    Hypo = get_hyp(pi_max, pi1c, NUMHYP)
    # Convert to int
    Hypo = Hypo.astype(int)
    num_alt = np.sum(Hypo)

    if 'mem' not in proc_list[FDR]:
        mempar_this = 1 
        m_corr = 1
    else:
        mempar_this = mempar
        m_corr = m_corr[0]

    #####  Parse abstinence vector - either there is enforced abstinence by this vector (or adaptive if wealth is low)
    abs_vec = np.zeros(NUMHYP)

    ##### Generate prior and penalty weights (if correlated, do so) or parse from file
    
    prw_vec = create_pr(prw_style, prw_const, m_corr, Hypo, NUMHYP)
    penw_vec = create_pen(penw_style, penw_const, prw_vec, NUMHYP)

       
    #### Set file and dirnames ##########
    dir_name = './dat'
    pr_filename = 'PR_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.1f_NR%d_%s' % (mu_gap, sigma, FDR, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_this, NUMHYP, NUMDRAWS, pi_max, NUMRUN, time_str)
    ad_filename = 'AD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_NR%d_%s' % (mu_gap, sigma, FDR, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_this, NUMHYP, NUMDRAWS, pi_max, NUMRUN, time_str)
    td_filename = 'TD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_NR%d_%s' % (mu_gap, sigma, FDR, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_this, NUMHYP, NUMDRAWS, pi_max, NUMRUN, time_str)
    
    # ----------- initialize result vectors and mats ------------- ##
    pval_mat = np.zeros([NUMHYP, NUMRUN])
    rej_mat = np.zeros([NUMHYP, NUMRUN])
    alpha_mat = np.zeros([NUMHYP, NUMRUN])
    wealth_mat = np.zeros([NUMHYP, NUMRUN])
    falrej_vec = np.zeros(NUMRUN)
    correj_vec = np.zeros(NUMRUN)
    totrej_vec = np.zeros(NUMRUN)
    FDR_mat = np.zeros([NUMHYP, NUMRUN])
    TDR_mat = np.zeros([NUMHYP, NUMRUN])
    falrej_mat = np.zeros([NUMHYP, NUMRUN])
    correj_mat = np.zeros([NUMHYP, NUMRUN])
    memFDR_mat = np.zeros([NUMHYP, NUMRUN])

    # ----------------- Run experiments (with same mus) ------------------ # 
    for l in range(NUMRUN):

        sigma_mu = np.sqrt(2*np.log(NUMHYP))

        #%%%%%%%%  Initialize theta_j, FDR and experiments %%%%%%%%#
        # Some random seed
        if (rndseed == 1):
            rndsd = l+50
        else:
            rndsd = None

        # Create a vector of gaps (i.i.d. from N(0,\sigma^2)
        if mu_gap == 0:
            gap = np.random.randn(NUMHYP)*sigma_mu
        else:
            gap = np.random.randn(NUMHYP)*mu_gap


        # Initialize FDR
        if FDR == 0:
            proc = BONF_proc_batch(alpha0, NUMHYP)
        elif FDR == 1:
            proc = GAI2_proc_batch(alpha0, NUMHYP, startfac)
        elif FDR == 2:
            proc = LORD2_proc_batch(alpha0, NUMHYP, startfac)
        elif FDR == 3:
            proc = GAI2_MW_proc_batch(alpha0, NUMHYP, startfac, prw_vec, penw_vec, mempar)
        elif FDR == 4:
            proc = GAI2_MW_AD_proc_batch(alpha0, NUMHYP, startfac, prw_vec, penw_vec, mempar, abs_vec, abs_eps, wealth_eps)
            

        # Initialize roexp: Alternative hypothesis is mean = mu_alt = gap (- or + is symmetric)
        mu_alt = gap 
        # mu_alt = 0 # If you want same mu for alternative as well!
        this_exp = rowexp_new_batch(NUMHYP, NUMDRAWS, Hypo, 0, gap)

        #%%%%%%%%% Run experiments: Get sample and p-values etc. %%%%%%%%%%%%%
        # Run random experiments with same random seed for all FDR procedures
        if mod_choice == 1:
            this_exp.gauss_two_mix(sigma, rndsd)
        elif mod_choice == 2:
            this_exp.bernoulli_draws(rndsd)
        pval_mat[:, l] = this_exp.pvec

        #%%%%%%%%%% Run FDR, get rejection and next alpha %%%%%%%%%%%%
        rej_mat[:, l] = proc.run_fdr(this_exp.pvec)
        alpha_mat[:, l] = proc.alpha    # Note this starts at alpha[1]
        wealth_mat[:, l] = proc.wealth_vec[0:NUMHYP] # Note this starts at W[0]

        #%%%%%%%%%%  Save results %%%%%%%%%%%%%%
        # total measures
        falrej_singlerun = np.array(rej_mat[:,l])*np.array(1-Hypo)
        correj_singlerun = np.array(rej_mat[:,l])*np.array(Hypo)
        totrej_singlerun = np.array(rej_mat[:,l])
        falrej_vec[l] = np.sum(falrej_singlerun)
        correj_vec[l] = np.sum(correj_singlerun)
        totrej_vec[l] = np.sum(totrej_singlerun)
        falrej_mat[:, l] = falrej_singlerun
        
        
        # Compute time dependent memFDR, FDR
        # Saves memFDR of the same parameter
        memFDR_mat[:, l] = compute_memFDR(totrej_singlerun, falrej_singlerun, penw_vec, mempar)
        FDR_mat[:, l] = compute_memFDR(totrej_singlerun, falrej_singlerun, penw_vec, 1)


        
    # -----------------  Compute average quantities we care about ------------- #
    FDR_vec = np.true_divide(falrej_vec, [max(totrej_vec[i],1) for i in range(len(totrej_vec))])
    TDR_vec = np.true_divide(correj_vec, num_alt)
    FDR_vec = [FDR_mat[NUMHYP-1][l] for l in range(NUMRUN)]
    memFDR_vec = [memFDR_mat[NUMHYP-1][l] for l in range(NUMRUN)]

    if verbose == 1:
        print("done with computation")

    # Save data
    saveres(dir_name, td_filename, np.r_[FDR_mat, memFDR_mat, rej_mat, falrej_mat,  wealth_mat, pval_mat, alpha_mat])
    saveres(dir_name, ad_filename, [TDR_vec, FDR_vec, memFDR_vec])
    
