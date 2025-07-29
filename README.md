# thermo-diffusive_sfd
This code is extension of previous thermo-diffusive code. In this the periodic  oscillations of flame arising due to diffussive thermal imbalance or Lewis number effects can be damped to obtain 
a good approximation of steady-state or stationary solution. In this way a fast solution to steady state can be obtained. The algorithm used is the Selective Frequency Damping algorithm described in the following paper:
Bastien E. Jordi, Colin J. Cotter, Spencer J. Sherwin; An adaptive selective frequency damping method. Physics of Fluids 1 September 2015; 27 (9): 094104. https://doi.org/10.1063/1.4932107
This approximation can be refined using the Newton solver.

In order to run, A list of input paramters are required and is described in inputs2d_info.txt
Also a sample input file inputs2d.txt is also included.
to compile use 
gfortran -o 2d tdcode_2d.f90
./compile.sh
to run use
./2d<inputs2d.txt
inlet.dat contains the mass fraction of fuel and oxidiser at the inlet and is generated using a mathematica code using the formula for premxideness
a sample field file output_inc.dat is required to start the code. 
The output is stored in output_filter.dat
a norm.dat file is registered to track the convergence
