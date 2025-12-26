TrackMomentum

Contact: Shane Horner
Email: sh22f@fsu.edu

OVERVIEW:

These macros serve as machinery to analyze the momentum resolution of sPHENIX
using K-Shorts by perturbing simulation pT. They plot mass resolution as a function of
pT so that the perturbed simulations can be compared to the data. The pT resolution has
two independent components - the multiple scattering term, and the detector resolution
term. This machinery runs over a grid of possible combinations of each term, and generates
a contour plot of the combinations that most closely replicate the data.

MY WORKING DIRECTORY:
/sphenix/user/shanehorner/macros/detectors/sPHENIX/unbinned
Can be used as an example of how to implement this tool effectively. 

TOOLS:

roo_DeltaPAnalysis:

	Runs over KshortReconstruction output files using kshort_variables. Perturbs the x
and y components of each pion in each reconstructed pair, then re-reconstructs each pair
with the perturbed vectors, and re-calculates invariant mass and pT. Each perturbed mass
and respective pT is saved into a RooDataSet as a point.


multirun_analyze_kshort_mom:

	Runs over Tony Frawley's kshort data files to create and save a RooDataSet of Mass
vs pT using kshort_variables. 


single_plot_analyze_kshort_mom:

	Groups the outputs from a specified RooDataSet into pT slices. Plots and saves a
histogram of reconstructed mass, and fits the mass width for each slice using a
gaussian+exponential function in RooFit. Plots and saves a TGraphErrors containing each pT
slice and its corresponding mass width. Use this to get the TGraphErrors for the data. It
can also be used to inspect individual perturbed simulations more thoroughly.
	

fast_bulk_plot_analyze_kshort_mom:

	This is a condensed version of single_plot_analyze_kshort_mom designed to run over
the entire grid and save a TGraphErrors of mass width as a function of pT for each
perturbation.


chi_square_kshort_2d:

	Generates chi-square surface (TGraph2D) by comparing the simulations to the data.
Extracts contours and finds the grid point which minimize the chi-square value. Saves the
contours and TGraph2D as plots. Outputs the lowest chi-sq value in the terminal and the
corresponding file to use for compare_kshort_simvdata.

compare_kshort_simvdata:

	Sanity check - use this to overlay TGraphErrors of the data with your perturbed
TGraphErrors of choice. For example, chi_square_kshort_2d will output the filename which
matches the data closest, so you can compare the data with the file it gives you. 


plot_sim_width:

	(Extra tool) Extracts the widths in slices of the raw simulation data without any
perturbations, and fits a TF1 of the multiple scattering and detector resolution
parameters to delta pt vs pt of the simulations.


SCRIPTS:

runit_kshort_muddy:

	Feeds input argument to roo_DeltaPAnalysis, which determines the perturbation
parameters for that element of the grid. Then, runs fast_bulk_plot_analyze_kshort_mom with
the same input argument. 


condor_muddy_kshort_mom:

	Queue the number of elements in the grid (ngrid in roo_DeltaPAnalysis)^2 to run
runit_kshort_muddy over the entire grid.























	
	
	

