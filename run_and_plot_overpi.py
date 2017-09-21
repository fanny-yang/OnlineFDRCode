import logging, argparse
import numpy as np
from exp_FDR_batch_new import*
from plot_batch_results import*
from toimport import *


def main():

    if not os.path.exists('./dat'):
        os.makedirs('./dat')

    #########%%%%%%  SET PARAMETERS FOR RUNNING EXPERIMENT %%%%%%%##########

    FDRrange = str2list(args.FDRrange)
    if not hasattr(args, 'mempar_control'):
        mempar_control = args.mempar
    else:
        mempar_control = args.mempar_control
    pirange = str2list(args.pirange, 'float')

    #########%%%%% SET PARAMETERS FOR RUNNING PLOTTING  ##################### 

    hyprange = [0]
    plot_numrun = 1 # Plot how many of the trials that are run

    ########%%%%%%%%%%%%%%%%% RUN EXPERIMENT %%%%%%%%########################
    
    for pi_max in pirange:
        if args.weights == 1:
            FDR = FDRrange[0]
            m_corr_range = str2list(args.m_corr, 'float')
            for m_corr in m_corr_range:

                # Prevent from running if data already exists
                filename_pre = 'AD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_NR%d' % (args.mu_gap, 1, FDR, 1, 1, 1, 0, m_corr, 0, mempar_control, args.num_hyp, 1, pi_max, args.num_runs)
                all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]

                # Run if data doesn't exist yet
                if all_filenames == []:
                    print("Running experiment for FDR %d correlation  %.2f and pi %.1f" % (FDR, m_corr, pi_max)) 
                    run_single(args.num_runs, args.num_hyp, 1, args.mu_gap, 1,  1, 1,  0, [m_corr], 0, pi_max, mempar_control, args.alpha0, args.mod_choice, FDR, sigma = 1, verbose = False)
                else:
                    print("Experiments are already run for correlation %.2f and pi %.1f" % (m_corr, pi_max))

        else:    
            # Run single FDR 
            for FDR in FDRrange:

                if 'mem' not in proc_list[FDR]: 
                    mempar_control = 1 
                    m_corr = 1
                else:
                    m_corr = float(args.m_corr)

                # Prevent from running if data already exists
                filename_pre = 'AD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_NR%d' % (args.mu_gap, 1, FDR, args.penw_style, args.penw_const, args.prw_style, args.prw_const, m_corr, 0, mempar_control, args.num_hyp, 1, pi_max, args.num_runs)
                all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]

                # Run experiment if data doesn't exist yet
                if all_filenames == []:
                    print("Running experiment for FDR procedure %s and pi %.1f" % (proc_list[FDR], pi_max)) 
                    run_single(args.num_runs, args.num_hyp, 1, args.mu_gap, args.penw_style,  args.penw_const, args.prw_style,  args.prw_const, [m_corr], 0, pi_max, mempar_control, args.alpha0, args.mod_choice, FDR, sigma = 1, verbose = False)
                else:
                    print("Experiments for FDR procedure %s are already run %.1f" % (proc_list[FDR], pi_max))

    # Plot different measures over hypotheses for different FDR
    if args.plot == 1:
        print("Now plotting ... ")
        if args.weights == 1:
            plot_results(3, 5, 0, [FDR], pirange, mempar_control, hyprange, args.mu_gap, 1, 1, 1, 1, 0, m_corr_range, 0, args.mempar, args.num_hyp)
        else:
            plot_results(args.plot_style, 5, 0, FDRrange, pirange, mempar_control, hyprange, args.mu_gap, 1, args.penw_style, args.penw_const, args.prw_style, args.prw_const, [m_corr], 0, args.mempar, args.num_hyp)
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--FDRrange', type=str, default = '1,2,0')
    parser.add_argument('--mempar',type=float, default = 1)
    parser.add_argument('--num-runs', type=int, default = 200)
    parser.add_argument('--num-hyp', type=int, default = 1000)
    parser.add_argument('--plot-style', type = int, default = 2)
    parser.add_argument('--prw-style', type=int, default = 3)
    parser.add_argument('--prw-const', type=int, default = 1)
    parser.add_argument('--m-corr', type=str, default = '1')
    parser.add_argument('--penw-style', type=int, default = 1)
    parser.add_argument('--penw-const', type=int, default = 1)
    parser.add_argument('--alpha0', type=float, default = 0.05)
    parser.add_argument('--mu-gap', type=float, default = 0)
    parser.add_argument('--mod-choice', type=int, default = 1)
    parser.add_argument('--pirange', type=str, default = '0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9')
    parser.add_argument('--weights', type=int, default=0)
    parser.add_argument('--plot', type=int, default=1)
    args = parser.parse_args()
    logging.info(args)
    main()
