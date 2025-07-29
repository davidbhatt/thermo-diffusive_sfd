# thermo-diffusive_sfd
This code is extension of previous thermo-diffusive code. In this the oscillations of flame arising due to diffussive thermal imbalance or Lewsi number effects can be damped to obtain 
a good approximation of steady-state or stationary solution. In this way a fast solution to steady state can be obtained.
This apprimation can be efined using the Newton solver.
  !More details can be found at following publication:                                                                  
  David S. Bhatt, Daniel Rodríguez, Linear analysis of thermo-diffusive instability from edge flames to fully premixed laminar flames with a wide range of Damköhler number and Lewis number greater than unity. (submitted to Combustion and Flame)

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
