# TOol For Fast Error Estimation (TOFFEE)
TOFFEE is a Python3 tool that automates cross-section uncertainty analysis using the tool available in MCNP6.3 and NJOY. The program utilizes the PERT and KSEN methods to generate sensitivities, uses NJOY to generate covariance data, and uses NUMPY to do the sandwich rule to calculate uncertainites.

# How to install and test TOFFEE
1. Download and extract TOFFEE_package.zip into a folder.
2. Download the ENDF6 format nuclear data library that you will be using in your MCNP evaluation. This code was designed using ENDF-B/VIII.0 data.This library should be within a folder that is named `endf_neutron_libraries` by default. An example for hydrogen-1 `endf_neutron_libraries/n-001_H_001.endf` which must be within the folder you extract the TOFFEE_package.  
3. I have included two example bash scripts that can be used to submit the new generated input file/files. These files should be editied to match your job submission system. If you submit MCNP files with the command line, you can look at these scripts to determine the files you should evaluate and the name of the output files. For KCODE files, 1 file will be evaluated. For SDEF files, several files will need to be evaluated within the `working_dir` folder.
4. To start TOFFEE, use this command `python3 submit.py` to evaluate the covariance library and input file.
5. Run the generated input file `perturbed_mcnp.inp` with the output `perturbed_mcnp.out`.
6. Once the simualtion is complete, use this command `python3 analyze.py` to generate a text file with all variance data, plots of all covariance matrices and sensitivity vectors, and a plot with top contributors to the total response variance.


# How to use TOFFEE and available options
1. Add your MCNP file for to the folder. The file `ksen_test.inp` is an example input file. The only addition to your input file will be two comment lines that identify the beginning and ending of your materials. The comments are `c Materials` and `c End Materials`, respectively.
2. In the file `submit.py`, you can edit the name of the input file in line 21.
3. If you would like a different energy discretization (Default is SCALE 56 group structure), change the vector in line 25.
4. If you would like a different list of reaction, change the vectors in line 33 and 36.
5. The plotting options can be turned off in lines 25, 27, 32, and 35 in the file `analyze.py` for the plots you would like the code to make. 
