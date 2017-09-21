--------------------------- README --------------------

HowTo reproduce the plots in the paper (with the settings mentioned in the paper), plots are saved as .pdf files in the folder "plots", data in the folder "dat" saved as .dat files

1. Plot of FDP and Power of LORD++, LORD and Bonferroni procedures to illustrate higher power for GAI++ compared to simple GAI procedures

python run_and_plot_overpi.py

Plots will be saved as FDRvsPI_MG0.0_Si1.0_PEW1_PEWC1_PRW3 and PowervsPI_MG0.0_Si1.0_PEW1_PEWC1_PRW3

2. Plot FDR and Power with correlated prior weights for LORD++ for different correlation coefficients

python run_and_plot_overpi.py --m-corr 0.95,1.0,1.05 --weights 1 --FDRrange 3

Plots will be saved as FDRvsPIvsMC... and PowervsPIvsMC...

3. Plot memFDR, memPower, wealth and alpha against time/hypothesis for memFDR++

python run_and_plot_overhyp.py --FDRrange 3,1

Plots will be saved as memFDR0.99vsHP_Plot0..., memTDR0.99vsHP_Plot0..., WealthvsHP_Plot0..., alphavsHP_Plot0 etc. 

4. Alpha-death memFDR, power and wealth plots

python run_and_plot_overhyp.py --FDRrange 3,4 --plot-style 1 --pi-max 0.01 

Plots will be saved as memFDR0.99vsHP_Plot1..., memTDR0.99vsHP_Plot1..., WealthvsHP_Plot1..., alphavsHP_Plot1 etc.

--------------------------------------------------------

Control and vary settings of experiments and plots by passing arguments for the above scripts as follows:

Relevant parameters to set for run_and_plot_xxx:
--FDRrange (string of the form 'x,y,z' where x,y,z are integers): Which FDR procedure to choose, with 
	0: Bonferroni, 1: LORD++, 2: LORD, 3: mem-LORD++, 4: mem-LORD++ (abstinent)
--mempar (float): Memory parameter \delta
--num-runs (integer): Number of independent random draws of one experimental setting (over which you want to average)
--num-hyp (integer): Number of hypotheses per experiment
--plot_style (integer): 0: memory quantities vs. number of hypothesis, 2: Power/FDR vs. pi1 for different FDR, 3: Power/FDR vs. pi1 for different prior weights 
--prw-style (integer): 1: Correlated with m_corr, 2: Anti correlated with m_corr, 3: Constant; 4: random!
--prw-const (integer): Constant for option prw_style = 3 (only relevant when this is set)
--m-corr (string of the form 'x,y,z' where x,y,z, are floats): Correlation coefficient of weight with true indicator whether hypothesis is null.
--penw-style (integer): 1: Constant, 2: Linearly decreasing, 3: Linearly increasing, 4: Correlated
--penw-const (integer)
--alpha0 (float): desired false discovery rate
--mu-gap (float): Gap between null and alternative mean, usually just \sqrt{2 log (# Hyp)}
--mod-choice (int): Type of probability distribution for each hypothesis 1: Gaussian, 2: Bernoulli draws
--plot (binary): 0: Only run experiments 1: Run experiments (only if not run yet) and plot
used only for ...._overpi.py:
--pirange (int): Range of pi1 for which you want to run and plot experiments for different weights
--weights (binary): Whether we want weight experiments plots
used only for ..._overhyp.py:
--pi-max (float): fixed pi1 (set 0 if pi1c should be used)
--pi1c (integer): Index of the choice of progression of pi1 (evenly spread over time) as in the corresponding line number in expsettings/PIC.dat. Arbitrary progressions can be added.


Relevant files used in run_and_plot_xxx:

exp_FDR_batch_new.py: Contains function run_single() to run NUMRUN of experiments with FDR procedure with a fixed single setting
rowexp_new_batch.py: Contains functions gauss_two_mix() and bernoulli_draws() which draws p-values for each hypothesiso
plot_batch_results.py: Contains plot_results() to plot plot_style 

Utilities: 
Put them all in one settings file. 
settings_util.py: 
	generate_hyp: creates pi1 vector, draws multiple samples of hypothesis vectors and saves in a file 
	get_hyp: Loads hypothesis vectors from file and outputs one
	create_pen: Create vector with penalty weights according to some style penw_style (see options above)
	create_pr: Create vector with prior weights according to some prw_style (see options above)



