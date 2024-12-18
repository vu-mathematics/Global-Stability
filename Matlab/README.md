# README for qORAC matlab files

	Copyright 2017, 
	Robert Planqué, Dept of Mathematics, Vrije Universiteit Amsterdam

	Recreates the images in 
	"Maintaining maximal metabolic flux by gene expression control" by 
	Planque, Hulshof, Teusink and Bruggeman, Cell Systems.
	
	The code is available at
	http://www.few.vu.nl/~rplanque/Research/qORAC/


First, within the directory in which you save the matlab files,
create two directories, "pics" and "Symbolic_functions". 
The first contains output figures, the second contains
the matlab symbolic functions with which the programs do most computations.
These functions are specific for each pathway, and need
to be generated each time parameters, kinetics or pathway structure are changed.

Note that the code files are not generic code for arbitrary pathways. 
Each file specifies the qORAC framework for a specific pathway. 
In future releases we plan to make this code more generic, 
so that it is easier to play around with different pathways. 

### Figure 2

To recreate the enzyme synthesis relations in **Figure 2B** for the network shown in Figure 2A, run

<code> >> sensor\_optimum\_enz\_CN\_newton</code>

The corresponding image is <code> CN\_1\_input\_output.jpg </code>

To recreate the subplots in **Figure 2C**, run

<code> >> [t,x,x0] = daes_CN(1,1,1,50);</code>

This gives several figures, called <code>CN\_1\_....jpg</code>. 

To rerun the script without computing the symbolic functions again, do

<code> >> [t,x,x0] = daes_CN(1,1,0,50);</code>

### Figure 3

To recreate **Figure 3A**, use

<code> >> daes\_extra\_param</code>

This example uses metabolites 2, 3 and 6 as sensors. 

For **Figure 3B**, we choose different sensors (namely 2, 4 and 6). These do follow the optimum (so the sensors are valid), but for which the combined pathway + qORAC does not steer to optimum, do 

<code> >> daes\_extra\_param\_wrong</code>

### Figure S1

For the complicated double branched pathway, with two allosteric interactions, Figure S1 in the Supplement, run

<code>[t,x,x0] = daes\_double\_branched\_sym;</code>

### Figure S2

This is the minimal ICs plot in the Supplement. 

<code> >> [t,x,v] = daes\_CN\_minimal\_ICs(1,1,0,50);</code>


### Figure S3

The code to make Supplementary Figure S3 is

<code> >> daes\_linearchain\_reversal</code>

You may have to run the program more than once to calculate through the critical thermodynamic equilibrium point (despite the fact that the equations have been changed appropriately to make the Implicit Function Theorem valid in this point). In each new run, a different initial condition is chosen. A good run has time running from 0 to 300. After running the code once, you can run it faster by not computing the symbolic functions again, using

<code> >> daes\_linearchain\_reversal(1)</code>
