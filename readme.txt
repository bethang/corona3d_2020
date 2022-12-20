corona3d_2020 aims to be a flexible planetary exosphere modeling program, able to handle
multiple particle types with a variety of initial distribution options and background
species profiles. See the sample 'corona3d_2020.cfg' file under the 'src' folder for
details on the currently available configuration parameter options.

**********************
** QUICKSTART GUIDE **
**********************

This software has been tested on both Linux and Mac systems using the GCC compilier.
To get started, first download the source code to the desired directory from a terminal:

     cd <desired directory>
     git clone https://github.com/rdelliott76/corona3d_2020.git

Now cd into the src directory and compile the application:

     cd corona3d_2020/src/
     make

Note that you may receive a compilation error due to missing the Eigen library. If this
happens, try this on Ubuntu:

     sudo apt install libeigen3-dev

On Mac:

     brew install eigen

You may need to look up how to install the Eigen library on your particular system if
neither of these commands works. Once installed, retry make, and it should work now. The relevant
'include' line may also need to be commented in in Particle.hpp.

Once installed, modify the corona3d_2020.cfg file to have the desired parameters. To
run the program, execute the following command from within the src directory:

     ./corona3d_2020



*************************
** SOURCE CODE DETAILS **
*************************

Following is a (somewhat) detailed explanation of the source code files and program logic
to aid future developers who would like to contribute to this software. There are also
fairly extensive comments throughout all of the source code files to aid in understanding.
A general knowledge of object-oriented programming concepts in C++ will be immensely helpful.

** GENERAL PROGRAM LOGIC OVERVIEW **

The first noteable bit of code that gets executed upon running the program is in the
Common_Functions.cpp file, specifically the get_seed() function. This bit of code seeds the
random number generator. It first looks for a file named "rng_seed" within the same directory
as the executable. This file should contain a single long long integer if you would like to run
a replicateable simulation for testing purposes. If the "rng_seed" file does not exist, then
the current system clock time is used as the seed. Execution now passes to main.cpp.

The main.cpp source file begins by reading in the main configuration file. This file must be
named "corona3d_2020.cfg" and must be in the same directory as the executeable. After reading
in all of the needed paramerters, main.cpp then constructs all of the necessary objects to
run the simulation. These objects are constructed and initialized in the following order: Planet,
array of Particles (actually, an array of shared pointers to Particle child classes, usually O or H),
Distribution (actually, a shared pointer to one of the Distribution child classes, usually
Distribution_Hot_H or Distribution_Hot_O), Background_Species, Atmosphere.

Once all the necessary objects are constructed, the run_simulation method within Atmosphere.cpp is called.
This routine executes the main loop of the program. More specifically, it consists of an outer for loop
that iterates for the number of time steps specified in the main configuration file. Within this loop,
there is a nested for loop that cycles through all of the active particles in the simulation. For each time
step, for each active particle, several things happen:  First, any console output or particle trace data is
output, if necessary. Then, the radial density and EDF counts are updated. Next, the particle's position is
updated according to its current velocity vector by timestepping the 3-dimensional equation of motion using
the value of dt supplied in the main configuration file. Then, the particle's velocity vector is updated using
the acceleration due to gravity at its current altitude. Next, a collision probability check is performed using
the background species' densities at the current altitude along with the available energy-dependent collision
cross sections supplied in the configuration files. If it is determined that a collision occurred, then the
particle's velocity vector is again updated using convervation of momentum (only elastic collisions are considered
currently). Finally, a deactivation check is done using three deacitivation criteria.

***************************************************************
** PLOTTING SCRIPTS FOR "HCO+ DISSOCIATIVE RECOMBINATION:    **
** A SIGNIFICANT DRIVER OF NONTHERMAL HYDROGEN LOSS AT MARS" **
***************************************************************
The plotting scripts are in src/plotting_scripts and must be run from inside that directory.
They can be run via the Terminal using python3 by running 'python3 name_of_script.py' or 'python name_of_script.py'.
The files as saved here produce the figures included in 'HCO+ dissociative recombination: A significant driver of nonthermal hydrogen loss at Mars.'
Further details for each script can be found below.

1) Figure 1, panels (a)-(c): Use plot_Figure1abc.py.
   Plots density and production rate profiles from input density and temperature csv files.
   Input files are located in src/inputs/Mars/.

   Function inputs: SA_list = list of solar activity condition strings (strings are either 'LSA' or 'HSA').
   Calls functions 'find_bg_dens' and 'find_col_dens_above' from Common_Plotting_Functions.py.
   Calls function 'find_rate_coefficient' from same file (function at bottom of plot_Figure1abc.py).

   find_bg_dens(SA, planet).
   Inputs: SA = solar activity string (either 'LSA' or 'HSA'); planet: planet string (either 'Mars' or 'Venus').
   Reads background density input files and returns 'bg_z' (list of altitudes for density data (cm)) and
   'CO2_density', 'O_density', 'CO_density', and 'N2_density' (lists of densities for each background species (cm-3)).

   find_col_dens_above(bg_z,CO2_density,O_density,CO_density,N2_density,z_new,bin_size,SA,planet).
   Inputs: 'bg_z' = list of altitudes for background density data from input file (cm).
           'CO2_density', 'O_density', 'CO_density', and 'N2_density' = lists of densities for each background species (cm-3).
	   'z_new' = list containing new altitude grid (cm).
	   'bin_size' = integer or float, distance between uniformly-separated altitude values in z_new (cm).
	   'SA' = solar activity string ('LSA' or 'HSA').
	   'planet' = planet string ('Mars' or 'Venus').
   Returns list of floats for column density above each altitude in z_new, of same length as z_new.

   find_rate_coefficient(mechanism,T_type,Tn_new,Te_new,Ti_new,z).
   Inputs: 'mechanism' string for reaction.
   	   'T_type' = string for temperature required for rate coefficient (options are 'Tn' (neutral), Ti (ion), or Te (electron)).
   	   'Tn_new' = if T_type is Tn, this is a list of neutral temperatures interpolated for altitude grid z_new.
	   'Te_new' = if T_type is Te, this is a list of electron temperatures interpolated for altitude grid z_new.
	   'Ti_new' = if T_type is Ti, this is a list of ion temperatures interpolated for altitude grid z_new.
	   'z' = altitude at which function calculates rate coefficient.
   Returns rate coefficient for selected mechanism at altitude z.

2) Figure 1, lower right panel: Use plot_particle_positions.py.
   This script utilizes output data from src/model_output_data to plot the positions of the test particles at the
   beginning and end of the simulation.

   It currently uses files from a specific low solar activity output directory, for which the particle positions were tracked.
   The timestep at which the positions are plotted can be changed in lines 21 and 50.

3) Figure 2: Use plot_Figure2.py.
   This script uses output data from src/model_output_data/LSA/1/ and src/model_output_data/HSA/1/ to plot the densities
   and energy distributions of the hot H particles from one low and one high solar activity model.
   The chosen output files can be specified in L29-30.

   The script plot_Figure2.py calls 'density_integration.py,' which is used to calculate line-of-sight densities for
   brightness calculations.
   
   Within density_integration.py:
   nadir_density(sc_alt,filename).
   Inputs: sc_alt = spacecraft altitude (km); filename = output file name with 1d density data.
   Returns column density (cm-2) and brightness (Rayleighs) along line-of-sight looking down towards the planet.

   limb_density(filename).
   Inputs: filename = output file name with density data along line-of-sight tangential to planet (see Supplementary Material).
   Returns column density (cm-2) and brightness (Rayleighs) along limb line-of-sight (cm-2).
 
4) Figure 3: Use plot_Figure3.py.
   This script plots the escape probability curves by calling the python script 'escape_prob_curves.py.'
   It uses data from model_output_data/escape_probabilities_LSA.csv and model_output_data/escape_probabilities_HSA.csv
   to calculate and plot the fitted escape probability curves from the escape probability model output data.

   escape_probabilitity_curves(return_or_plot).
   Input: 'return_or_plot' = string (either 'Return' or 'Plot').
   Using 'Return' returns lists of escape rates varying with altitude and lists of corresponding altitudes.
   Using 'Plot' plots escape probability profiles, as in Figure 3.

Figures are saved in src/plotting_scripts/figures/

Supplementary Tables:

- The values in Table S1 correspond to the 'Total loss rate' output, printed to the Terminal after each run and saved in model_output_data/../../output/console.out.

- The values in Table S2 are escape fractions that were calculated and manually input to model_output_data/escape_probabilities_LSA.csv and model_output_data/escape_probabilities_HSA.csv after each run.





