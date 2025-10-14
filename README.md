# PBS3
For the third assignment of Particle Based Simulations given at the Eindhoven University of Technology, it was tasked to implement code to implement code to make a simulation technique called Dissipative Particle Dynamics. The following file describes the made files by the students for the simulation technique.

# Running the simulation
In the 'main.c' file, the main loop code is presented. On line 44 and 45, there are definitions that can be commented to save on computation time. Each of the definitions have different functions:
    - HISTOGRAM: This command includes the calculation for the velocity, density and volume fraction histograms and saving it to the correct output files. 
    - NUMPART_CALC: This command controls some of the parameters that can quickly be changed by the user to obtain data quickly. The data that can be changed starts from line 76 in the main.c' file. 

# Forces.c/.h
The forces (conservative, dissipative, random and bond forces) were worked out in their respective part and validated which is explained in the 'Particle Based Simulations - Group A3_06.pdf' file. 

# grf.c/.h
This file contains the code that computes the radial distribution function by using a histogram and printing the data to a file.

# grf.c/.h
Similar to the file described above, this file contains the code that computes the radial distribution function by using a histogram and printing the data to a file. This has been added to a different file to make sure no confusion between the students could arise what is calculated where.

# parametercalc.c/.h
In this file, the volume fraction of binary is calculated and stored in a histogram for later data visualization. In addition, the chi-interaction parameter is calculated from the obtained histogram.

# vector_functions.c/.h
This file contains basic vector calculations that is used for force calculations.

# plotting.ipynb
This file plots the data obtained by running the simulation. If parameters need to be altered for the visualization, the simulation needs to be run again to abtain the most accurate results.