### Import Python libraries
import numpy as np
from numpy import sqrt, log, exp, mean, cumsum, sum, zeros, ones, argsort, argmin, argmax, array, maximum, concatenate
from numpy.random import randn, rand
np.set_printoptions(precision = 4)

import os
#import matplotlib.pyplot as plt
import time
from datetime import datetime
import sys


### Import utilities for plotting
from plotting import*
from settings_util import*
from toimport import*

# Plot_styles:
# 0: memFDR over time
# 1: memFDR over time with alpha death
# 2: Power/FDR over pi
# 3: Power/FDR over pi with weights


def plot_results(plot_style, plot_numrun, whichrun, FDRrange, pirange, mempar_plot, hyprange, mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr_range, pi1c, mempar, NUMHYP, NUMDRAWS = 1, startfac = 0.1, alpha0 = 0.05):
# hyprange: range(2000,3000,1) or 0 (if you want entire range)
# pirange[]: needs to be [0] if no constant pi1 chosen

    plot_dirname = './plots'
    numrun = 100000

    #%%%%%%%%%%%%%%%%%%%%  PLOTS vs. Hyp (time)  %%%%%%%%%%%%%%%%%%%%%%

    # Time variant plots
    # Here there are hyp_steps > 1
    # Find all possible pi_max
    
    ######## Plots vs. Hyp diff FDR, fixed run, fixed mempar #########
    if plot_style == 0 or plot_style == 1:

        numFDR = len(FDRrange)
        
        pi_max = pirange[0]
        m_corr = m_corr_range[0]

        # ----------- LOAD DATA --------
        FDR_mat = [None]*len(FDRrange)
        memFDR_mat = [None]*len(FDRrange)
        wealth_mat = [None]*len(FDRrange)
        memTDR_mat = [None]*len(FDRrange)
        TDR_mat = [None]*len(FDRrange)
        alpha_mat = [None]*len(FDRrange)

        for FDR_j, FDR in enumerate(FDRrange):
            if 'mem' not in proc_list[FDR]:
                mempar_r = 1
                m_corr_r = 1
            else:
                mempar_r = mempar
                m_corr_r = m_corr


            filename_pre = 'TD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_' % (mu_gap, sigma, FDR, penw_style, penw_const, prw_style, prw_const, m_corr_r, pi1c, mempar_r, NUMHYP, NUMDRAWS, pi_max)
            all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]
            
            if all_filenames == []:
                print("No files found!")
                print(filename_pre)
                sys.exit()
            
            # Load results
            result_mat = np.loadtxt('./dat/%s' % all_filenames[0])
            FDR_vec = result_mat[0:NUMHYP, whichrun]
            rej_vec = result_mat[2*NUMHYP:3*NUMHYP, whichrun]
            falrej_vec = result_mat[3*NUMHYP:4*NUMHYP, whichrun]
            wealth_vec = result_mat[4*NUMHYP:5*NUMHYP, whichrun]
            alpha_vec = result_mat[6*NUMHYP:7*NUMHYP, whichrun]

            # Get true Hypo vector
            if pi1c > 0:
                pi_max = 0 # Just to make sure pi1c is beieng used when desired
            Hypo = get_hyp(pi_max, pi1c, NUMHYP)
            Hypo = Hypo.astype(int)
            
            # Generate prior and penalty weights
            prw_vec = create_pr(prw_style, prw_const, m_corr, Hypo, NUMHYP)
            penw_vec = create_pen(penw_style, penw_const, prw_vec, NUMHYP)

            # Compute memFDR with different parameters
            memFDR_vec = compute_memFDR(rej_vec, falrej_vec, penw_vec, mempar_plot)
            FDR_vec = compute_memFDR(rej_vec, falrej_vec, penw_vec, 1)
            TDR_vec = compute_memTDR(rej_vec, Hypo, penw_vec, 1)
            memTDR_vec = compute_memTDR(rej_vec, Hypo, penw_vec, mempar_plot)

            # Save to matrix
            FDR_mat[FDR_j] = FDR_vec
            memFDR_mat[FDR_j] = memFDR_vec
            wealth_mat[FDR_j] = wealth_vec
            TDR_mat[FDR_j] = TDR_vec
            memTDR_mat[FDR_j] = memTDR_vec
            alpha_mat[FDR_j] = alpha_vec
        
        # -------- PLOT ---------------
        # Set x axis
        if len(hyprange) == 1:
            xs = range(NUMHYP)
            hyplen = NUMHYP
        else:
            # Cut the matrices
            xs = hyprange
            hyplen = len(hyprange)
            FDR_mat = np.array(FDR_mat)[:,0:len(hyprange)]
            memFDR_mat = np.array(memFDR_mat)[:,0:len(hyprange)]
            wealth_mat = np.array(wealth_mat)[:,0:len(hyprange)]
            memTDR_mat = np.array(memTDR_mat)[:,0:len(hyprange)]
            TDR_mat = np.array(TDR_mat)[:,0:len(hyprange)]
            alpha_mat = np.array(alpha_mat)[:,0:len(hyprange)]
        
        legends_list = np.array(proc_list).take(FDRrange)[0:numFDR]
        
        
        ##### FDR vs. HYP #####
        
        if plot_style == 1:
            leg_col = 1
        else:
            leg_col = 2

        #### FDP vs HYP ###
        #        filename = 'FDRvsHP_Plot%d_MG%.1f_Si%.1f_PWE%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_HR%d_R%d' %  (plot_style, mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar, NUMHYP, NUMDRAWS, pi_max, hyplen, whichrun)

#        plot_curves_mat(xs, FDR_mat, legends_list, plot_dirname, filename,  'Hypothesis index', 'FDR($J$)', 0, leg_col = leg_col) 

        ##### memFDR (strictly memFDP) vs. HYP ####
        filename = 'memFDP%.2fvsHP_Plot%d_MG%.1f_Si%.1f_PWE%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_HR%d_R%d' %  (mempar_plot, plot_style,mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar, NUMHYP, NUMDRAWS, pi_max, hyplen, whichrun)

        plot_curves_mat(xs, memFDR_mat, legends_list, plot_dirname, filename,  'Hypothesis index', 'mem-FDP($J$) with $\delta = $%.2f' % mempar_plot, 0, leg_col = leg_col) 


        #### Wealth vs. HYP ####
        filename = 'WealthvsHP_Plot%d_MG%.1f_Si%.1f_PWE%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_HR%d_R%d' %  (plot_style, mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar, NUMHYP, NUMDRAWS, pi_max, hyplen, whichrun)
        plot_curves_mat(xs, wealth_mat, legends_list, plot_dirname, filename,  'Hypothesis index', 'Wealth($J$)', 0, leg_col = leg_col) 

        #### Power vs. HYP ####
#        filename = 'PowervsHP_Plot%d_MG%.1f_Si%.1f_PWE%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_HR%d_R%d' %  (plot_style, mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar, NUMHYP, NUMDRAWS, pi_max, hyplen, whichrun)

#        plot_curves_mat(xs, TDR_mat, legends_list, plot_dirname, filename,  'Hypothesis index', 'Power($J$)', 0, leg_col = leg_col) 

        #### memTDR vs. HYP ####
        filename = 'memPower%.2fvsHP_Plot%d_MG%.1f_Si%.1f_PWE%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_HR%d_R%d' %  (mempar_plot, plot_style, mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar, NUMHYP, NUMDRAWS, pi_max, hyplen, whichrun)

        plot_curves_mat(xs, memTDR_mat, legends_list, plot_dirname, filename,  'Hypothesis index', 'mem-Power($J$)', 0, leg_col = leg_col) 

        #### alpha vs. HYP ####
        filename = 'alphavsHP%.2fvsHP_Plot%d_MG%.1f_Si%.1f_PWE%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_HR%d_R%d' %  (mempar_plot, plot_style, mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar, NUMHYP, NUMDRAWS, pi_max, hyplen, whichrun)

        plot_curves_mat(xs, alpha_mat, legends_list, plot_dirname, filename,  'Hypothesis index', '$alpha(J)$', 0, leg_col = leg_col) 


    #%%%%%%%%%%%%%%%%%%%  PLOTS vs. pi1 (without weights) %%%%%%%%%%%%%%%%%%%%%%%%%%


    elif plot_style == 2:

        m_corr = m_corr_range[0]
        numFDR = len(FDRrange)
               
        # ---------- LOAD DATA --------------
        for FDR_j, FDR in enumerate(FDRrange):
            if 'mem' not in proc_list[FDR]:
                mempar_r = 1
            else:
                mempar_r = mempar

            filename_pre = 'AD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_' % (mu_gap, sigma, FDR, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_r, NUMHYP, NUMDRAWS)
            all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]

            if all_filenames == []:
                print("No file found!")
                print(filename_pre)
                sys.exit()
            
            # Get different pis
            pos_PM_start = [all_filenames[i].index('PM') for i in range(len(all_filenames))]
            pos_PM_end = [all_filenames[i].index('_NR') for i in range(len(all_filenames))]
            PM_vec = [float(all_filenames[i][pos_PM_start[i] + 2:pos_PM_end[i]]) for i in range(len(all_filenames))]

            order = np.argsort(PM_vec)
            PM_list = pirange

            # Initialize result matrices
            if FDR_j == 0:
                TDR_av = np.zeros([len(FDRrange), len(PM_list)])
                TDR_std = np.zeros([len(FDRrange), len(PM_list)])
                FDR_av = np.zeros([len(FDRrange), len(PM_list)])
                FDR_std = np.zeros([len(FDRrange), len(PM_list)])             
                memFDR_av = np.zeros([len(FDRrange), len(PM_list)])
                memFDR_std = np.zeros([len(FDRrange), len(PM_list)])

            # Merge everything with the same NA and NH
            for k, PM in enumerate(PM_list):
                indices = np.where(np.array(PM_vec) == PM)[0]
                result_mat = []
                # Load resultmats and append 
                for j, idx in enumerate(indices):
                    result_mat_cache = np.loadtxt('./dat/%s' % all_filenames[idx])
                    if (j == 0):
                        result_mat = result_mat_cache
                    else:
                        result_mat = np.c_[result_mat, result_mat_cache]
                
                if result_mat == []:
                    ipdb.set_trace()
                    break
                numrun = len(result_mat[0])
                # Get first vector for TDR
                TDR_vec = result_mat[0]
                TDR_av[FDR_j][k] = np.average(TDR_vec)
                TDR_std[FDR_j][k] = np.true_divide(np.std(TDR_vec),np.sqrt(numrun))
                # FDR
                FDR_vec = result_mat[1]
                FDR_av[FDR_j][k] = np.average(FDR_vec)
                FDR_std[FDR_j][k] = np.true_divide(np.std(FDR_vec), np.sqrt(numrun))
                # memFDR
                memFDR_vec = result_mat[2]
                memFDR_av[FDR_j][k] = np.average(memFDR_vec)
                memFDR_std[FDR_j][k] = np.true_divide(np.std(memFDR_vec), np.sqrt(numrun)) 

        
        # -------- PLOT ---------------
       
        xs = PM_list
        x_label = '$\pi_1$'
        legends_list = np.array(proc_list).take(FDRrange)[0:numFDR]
        target_vec = alpha0*np.ones(len(xs))


        ##### FDR vs. pi #####

        
        filename = 'FDRvsPI_MG%.1f_Si%.1f_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d' %  (mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_r, NUMHYP, NUMDRAWS)
        plot_errors_mat(xs, FDR_av, FDR_std, legends_list, plot_dirname, filename, x_label, 'FDR', plus = False)

        ##### TDR vs. pi ####
        filename = 'PowervsPI_MG%.1f_Si%.1f_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d' %  (mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_r, NUMHYP, NUMDRAWS)
        plot_errors_mat(xs, TDR_av, TDR_std, legends_list, plot_dirname, filename, x_label, 'Power', plus = False)


        #%%%%%%%%%%%%%%%%%%%  PLOTS vs. pi1 for different weights %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elif plot_style == 3:
        
        FDR = FDRrange[0]
        numMC = len(m_corr_range)
        mempar_r = mempar
       
        # ---------- LOAD DATA --------------
        for mc_j, m_corr in enumerate(m_corr_range):
            
            filename_pre = 'AD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_' % (mu_gap, sigma, FDR, penw_style, penw_const, prw_style, prw_const, m_corr, pi1c, mempar_r, NUMHYP, NUMDRAWS)
            all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]

            if all_filenames == []:
                print("No file found!")
                print(filename_pre)
                sys.exit()
            
            # Get different pis
            pos_PM_start = [all_filenames[i].index('PM') for i in range(len(all_filenames))]
            pos_PM_end = [all_filenames[i].index('_NR') for i in range(len(all_filenames))]
            PM_vec = [float(all_filenames[i][pos_PM_start[i] + 2:pos_PM_end[i]]) for i in range(len(all_filenames))]

            order = np.argsort(PM_vec)
            PM_list = sorted(set(np.array(PM_vec)[order]))

            # Initialize result matrices
            if mc_j == 0:
                TDR_av = np.zeros([numMC, len(PM_list)])
                TDR_std = np.zeros([numMC, len(PM_list)])
                FDR_av = np.zeros([numMC, len(PM_list)])
                FDR_std = np.zeros([numMC, len(PM_list)])             
                memFDR_av = np.zeros([numMC, len(PM_list)])
                memFDR_std = np.zeros([numMC, len(PM_list)])

            # Merge everything with the same NA and NH
            for k, PM in enumerate(PM_list):
                indices = np.where(np.array(PM_vec) == PM)[0]
                result_mat = []
                # Load resultmats and append 
                for j, idx in enumerate(indices):
                    result_mat_cache = np.loadtxt('./dat/%s' % all_filenames[idx])
                    if (j == 0):
                        result_mat = result_mat_cache
                    else:
                        result_mat = np.c_[result_mat, result_mat_cache]
                
                numrun = len(result_mat[0])
                # Get first vector for TDR
                TDR_vec = result_mat[0]
                TDR_av[mc_j][k] = np.average(TDR_vec)
                TDR_std[mc_j][k] = np.true_divide(np.std(TDR_vec),np.sqrt(numrun))
                # FDR
                FDR_vec = result_mat[1]
                FDR_av[mc_j][k] = np.average(FDR_vec)
                FDR_std[mc_j][k] = np.true_divide(np.std(FDR_vec), np.sqrt(numrun))
                # memFDR
                memFDR_vec = result_mat[2]
                memFDR_av[mc_j][k] = np.average(memFDR_vec)
                memFDR_std[mc_j][k] = np.true_divide(np.std(memFDR_vec), np.sqrt(numrun)) 


        
        # -------- PLOT ---------------
        xs = PM_list
        x_label = '$\pi_1$'
        
        # Create legend
        legends_list = []
        for m_j, m_corr in enumerate(m_corr_range):
            legends_list.append('a = %.2f' % (m_corr-1))
                            
        ##### FDR vs. pi #####

        filename = 'FDRvsPIvsMC_MG%.1f_Si%.1f_PEW%d_PEWC%d_PRW%d_PRC%d_PIC%d_MP%.2f_NH%d_ND%d' %  (mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, pi1c, mempar_r, NUMHYP, NUMDRAWS)
        plot_errors_mat(xs, FDR_av, FDR_std, legends_list, plot_dirname, filename, x_label, 'FDR')

        ##### TDR vs. pi ####
        filename = 'PowervsPIvsMC_MG%.1f_Si%.1f_PEW%d_PEWC%d_PRW%d_PRC%d_PIC%d_MP%.2f_NH%d_ND%d' %  (mu_gap, sigma, penw_style, penw_const, prw_style, prw_const, pi1c, mempar_r, NUMHYP, NUMDRAWS)
        plot_errors_mat(xs, TDR_av, TDR_std, legends_list, plot_dirname, filename, x_label, 'Power')


