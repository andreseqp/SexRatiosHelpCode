# Queen-worker conflict can drive the evolution of social polymorphism and split sex ratios

This repository provides the code to simulate the the co-evolution of helping
behaviour, sex-ratios and fratricidal behaviour. The model is
coded in c++ language. Visualization of the output is done in R. In this readme 
we describe how to run the simulations and how to plot their outcome. 

## Model simulations

As stated before, code for the simulation model is written in C++ language. 
The file contained at the base folder: random.h contain useful functions and 
random number generators used in the stochastic simulations. 
The source file of the simulation model is fratricide.cpp.
With these files executables can be compile for the simulation with standard c++
compilers. At the beginning of the main function in the cpp file parameters can 
be change prior to compilation. Execution of the program will produce the 
simulation files necessary plot the outcome. The program runs by default 10 
replicates for the parameters chosen. 

## Visualization

Figure 2 in the manuscript shows dynamics of the deterministic model (numerical
simulations), and stochastic simulations. The stochastic simulations are 
represented by choosing a representative run. But, the details of the different runs 
differ from each other. Thus, to get a good fit between numerical an stochastic 
simulations some parameters will have to be tweaked in the numerical 
simulations too. In particular variable *pertTime* must be set the time in which 
a transition from solitary to social life happen in the stochastic simulations. 
A part from this change, running the file *Figure2.R* should produce figure 2 of 
the manuscript. 

Figure 3, 4 and 5 should be produced by running the corresponding files 
*Figure3.R*, *Figure4.R* and *Figure5.R*.  
