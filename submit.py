import time
from TOFFEE import TOFFEE_class

main = TOFFEE_class()
time_1 = time.perf_counter()



#### USER innputs

# select method for producing covariance data
main.covariance_method = 'NJOY'

# Do you want to submit the MCNP file automatically
main.run_mcnp = False
### options:
    #True 
    #False
    
# Name of your mcnp input file
main.mcnp_file_name = 'ksen_test.inp'

# Your desired energy group structure
#SCALE 56 group
main.energy_bin_structure = [1.00E-11,4.00E-09,1.00E-08,2.53E-08,4.00E-08,5.00E-08,6.00E-08,8.00E-08,1.00E-07,1.50E-07,2.00E-07,2.50E-07,3.25E-07,3.50E-07,3.75E-07,4.50E-07,6.25E-07,1.01E-06,1.08E-06,1.13E-06,5.00E-06,6.25E-06,6.50E-06,6.88E-06,7.00E-06,2.05E-05,2.12E-05,2.18E-05,3.60E-05,3.71E-05,6.50E-05,6.75E-05,1.01E-04,1.05E-04,1.16E-04,1.18E-04,1.88E-04,1.92E-04,2.25E-03,3.74E-03,1.70E-02,2.00E-02,5.00E-02,2.00E-01,2.70E-01,3.30E-01,4.70E-01,6.00E-01,7.50E-01,8.61E-01,1.20E+00,1.50E+00,1.85E+00,3.00E+00,4.30E+00,6.43E+00,2.00E+01]

#SCALE 44 group
#main.energy_bin_structure = [1.00E-11,3.00E-09,7.50E-09,1.00E-08,2.53E-08,3.00E-08,4.00E-08,5.00E-08,7.00E-08,1.00E-07,1.50E-07,2.00E-07,2.25E-07,2.50E-07,2.75E-07,3.25E-07,3.50E-07,3.75E-07,4.00E-07,6.25E-07,1.00E-06,1.77E-06,3.00E-06,4.75E-06,6.00E-06,8.10E-06,1.00E-05,3.00E-05,1.00E-04,5.50E-04,3.00E-03,1.70E-02,2.50E-02,1.00E-01,4.00E-01,9.00E-01,1.40E+00,1.85E+00,2.35E+00,2.48E+00,3.00E+00,4.80E+00,6.43E+00,8.19E+00,2.00E+01]



# List of reactions you want evaluated
main.reactions_to_evaluate = ['rxn=2', 'rxn=4', 'rxn=16', 'rxn=18', 'rxn=102', 'rxn=103', 'rxn=104', 'rxn=105', 'rxn=106', 'rxn=107']

# Names for the reactions as seen in MCNP (reactions must be in the same location in the list)
main.reaction_names = ['elastic', 'inelastic','n,2n', 'fission','n,gamma','n,p', 'n,d', 'n,t', 'n,3he', 'n,alpha']


# Percent change in density of the nuclides being perturbed (ONLY for PERT method)
main.density_change = 1.0



main.run_mcnp_sensitivity()


print('done in ' + str(round(time.perf_counter()-time_1,3)) + ' seconds')
