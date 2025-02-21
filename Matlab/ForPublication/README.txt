README for Matlab files for 

'Biochemically sensing cellular growth rate facilitates its 
robust optimal adaptation to changing conditions'

by R. Planqué, J. Hulshof and F.J. Bruggeman

	Copyright 2025, 
	Robert Planqué, Dept of Mathematics, Vrije Universiteit Amsterdam



To recreate Figure 4 in the main text, run in Matlab

>> qorac(400,1,[2,9,0.1],1);

To explore different scenarios, update the parameters in qorac.m, 
and/or the definition of the metabolic pathway (at the end of
that script). 

You will then need to update the control laws. To do that, run

>> scan_x0

In this script, the nutrient environment is varied, and 
for each choice the optimal allocation is computed. The output files

alphaopts1.mat
muopts1.mat
x0opts.mat

contain the optimal ribosomal allocation (chi in the paper), optimal growth rates
(lambda in the paper), and corresponding environmental nutrient concentrations, resp. 

In qorac() you can also run the pathway with fixed allocation, using

>> qorac(400,[0.3 0.3 0.4],2,1);

which uses the allocation vector [0.3 0.3 0.4], and nutrient concentration = 2. 
