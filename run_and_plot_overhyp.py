import logging, argparse
import numpy as np

# Import own utilities
from exp_FDR_batch_new import*
from plot_batch_results import*
from toimport import *

def main():

    #########%%%%%%  SET PARAMETERS FOR RUNNING EXPERIMENT %%%%%%%##########

    FDRrange = str2list(args.FDRrange)
    if args.mempar_control is None:
        mempar_control = args.mempar
    else:
        mempar_control = args.mempar_control

    #########%%%%% SET PARAMETERS FOR RUNNING PLOTTING  ##################### 

    hyprange = [0]
    plot_numrun = 1 # Plot how many of the trials that are run
    whichrun = 1 # random really
    if args.num_runs < whichrun:
        logging.info("Can't choose a run that was not conducted")
        sys.exit()
    pirange = [args.pi_max]
    if args.pi_max != 0:
        args.pi1c = 0

    # Run single FDR and check over hypothesis
    for FDR in FDRrange:

        # Prevent from running if data already exists
        if 'mem' not in proc_list[FDR]: 
            mempar_control_this = 1 
        else:
            mempar_control_this = mempar_control

        filename_pre = 'TD_MG%.1f_Si%.1f_FDR%d_PEW%d_PEWC%d_PRW%d_PRC%d_MC%.4f_PIC%d_MP%.2f_NH%d_ND%d_PM%.2f_NR%d' % (args.mu_gap, 1, FDR, args.penw_style, args.penw_const, args.prw_style, args.prw_const, args.m_corr, args.pi1c, mempar_control_this, args.num_hyp, 1, args.pi_max, args.num_runs)
        if not os.path.exists('./dat'):
            os.makedirs('./dat')

        all_filenames = [filename for filename in os.listdir('./dat') if filename.startswith(filename_pre)]
        if all_filenames == []:
            print("Running experiment for FDR procedure %s" % proc_list[FDR]) 
            run_single(args.num_runs, args.num_hyp, 1, args.mu_gap, args.penw_style,  args.penw_const, args.prw_style,  args.prw_const, [args.m_corr], args.pi1c, args.pi_max, mempar_control_this, args.alpha0, args.mod_choice, FDR, sigma = 1, verbose = False, rndseed = args.randseed)
        else:
            print("Experiments are already run")
        
    # Plot different measures over hypotheses for different FDR
    print("Now plotting ... ")
    plot_results(args.plot_style, 5, whichrun, FDRrange, pirange, args.mempar, hyprange, args.mu_gap, 1, args.penw_style, args.penw_const, args.prw_style, args.prw_const, [args.m_corr], args.pi1c, mempar_control, args.num_hyp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--FDRrange', type=str, default = '3,1')
    parser.add_argument('--mempar',type=float, default = 0.99)
    parser.add_argument('--num-runs', type=int, default = 10)
    parser.add_argument('--num-hyp', type=int, default = 2000)
    parser.add_argument('--plot-style', type = int, default = 0)
    parser.add_argument('--prw-style', type=int, default = 3)
    parser.add_argument('--prw-const', type=int, default = 1)
    parser.add_argument('--m-corr', type=int, default = 1)
    parser.add_argument('--pi1c', type=int, default = 2)
    parser.add_argument('--pi-max', type=float, default = 0)
    parser.add_argument('--penw-style', type=int, default = 1)
    parser.add_argument('--penw-const', type=int, default = 1)
    parser.add_argument('--alpha0', type=float, default = 0.05)
    parser.add_argument('--mu-gap', type=float, default = 0)
    parser.add_argument('--mod-choice', type=int, default = 1)
    parser.add_argument('--mempar-control', type=float)
    parser.add_argument('--randseed', type=int, default=1)
    args = parser.parse_args()
    logging.info(args)
    main()
