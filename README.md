# Quasifree Dynamics for a Two Band Charge Transport System

This project calculates observables concerned with charge transport in a quasi-free two-band system: one *insulating* band and one *conducting* band. 

The project requires numpy and modules from the Python standard library.

## Running the Program

This program has the following optional arguments:

### Parameter Values

- L (default value: 100; int or list of ints) - the length of the system (including resevoirs) is then 2L+1
- l (default value: 2; int or list of ints)) - the length of the insulatling band (= the conducting band) is 2l+1
- epsilon (default value: 0.7; float or list of floats) - the strength of hopping between adjacent sites in the conducting band.
- kappa (default value: 0.85; float or list of floats) - the strength of hopping between adjacent sites in the resevoir.
- theta (default value: 0.85; float or list of floats) the strength of hopping between the system and the resevoir.
- mu_0 (default value: 0; float or list of floats)
- mu_12 (default value: 0.5; float or list of floats)
- mu_11 (default value: 0.5; float or list of floats) 
- gamma (default value: 0; float or list of floats) - the strength of hopping between the bands
- eta (default value: 0.1; float or list of floats) - the voltage
- beta (default value: 300; float or list of floats) - the temperature of the system in Kelvin

### Observables

 If the following values are set to True, evolution.py calculates:
 
- C (default value: False; Boolean) - Current
- CT (default value: False; Boolean) - Charge transport at individual sites
- CS (default value: False; Boolean) - Charge squared (used to calculate variance)
- QF (default value: False; Boolean) - Quantum fluctuation or current variance.

### Options

- AutoSR (default value: False; Boolean) - Calculates the mean of the observables by finding the stationary regime automatically (by increasing L if necessary)
- dirname (default value: False; string) - Creates a directory called dirname where the results are saved. If dirname is not given, evolution.py stores the results in a directory called *date*_quasifree_calcs where date is today's date. If the directory already exists, the user will be prompted to enter another directory name.

Upon running the program 

- With AutoSR set to False: The user will be prompted to add the name of a numpy file containing a list of times (which should be given in picoseconds) at which to calculate the chosen observables. A directory name dirname will be created along with .npy files containing the results of the specified observables. A metadata file is also generated containing the corressponding parameters.
- With AutoSR set to True: The stationary regime will be found automatically. A secondary metadata file with the times and value of L for which which the observables were calculated is created along with the times at which the observables were calculated. The program with calculate the mean of the observable in the stationary regime. The sequence of means will be recorded in a single numpy file.
- If PD (particle density) or CT (charge transport) are set to true the user will be prompted to enter a list of sites at which to calculate the these observables.

## Example usage

### Example 1    
`I=qf_calcs.Evolution(Ls=250, ls=[2,3,4,5], C=True, QF=True, AutoSR=False,dirname='testSR')` 

`please input the filename of a .npy file with desired times`

`times.npy`
 
 ### Example 2    
`I=qf_calcs.Evolution(etas=0.1, PD=True, AutoSR=False,dirname='testSR')` 

`please input a list of sites with sites at which to calculate the particle density`

`[-5,-2,0,2,5]`
