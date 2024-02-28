# TOol For Fast Error Estimation (TOFFEE)
TOFFEE is a Python3 tool that automates cross-section uncertainty analysis using the tool available in MCNP6.3 and NJOY. The program utilizes the PERT and KSEN methods to generate sensitivities, uses [NJOY21](https://github.com/njoy/NJOY21) to generate covariance data, and uses NUMPY to do the sandwich rule to propogate uncertainites.

# How to install and test TOFFEE
1. Download and extract TOFFEE_package.zip into a folder.
2. Download the ENDF6 format nuclear data library that you will be using in your MCNP evaluation. This code was designed using [ENDF-B/VIII.0](https://www.nndc.bnl.gov/endf-b8.0/download.html) data.This library should be within a folder that is named `endf_neutron_libraries` by default. An example for hydrogen-1 `endf_neutron_libraries/n-001_H_001.endf` which must be within the folder you extract the TOFFEE_package.  
3. I have included two example bash scripts that can be used to submit the new generated input file/files. These files should be editied to match your job submission system. If you submit MCNP files with the command line, you can look at these scripts to determine the files you should evaluate and the name of the output files. For KCODE files, 1 file will be evaluated. For SDEF files, several files will need to be evaluated within the `working_dir` folder.
4. The `NJOY_run_file.txt` is the file that evaluates NJOY. This bash script will needed to be editied to reference the location of your NJOY install.
5. To start TOFFEE, use this command `python3 submit.py` to evaluate the covariance library and input file.
6. Run the generated input file `perturbed_mcnp.inp` with the output `perturbed_mcnp.out`.
7. Once the simualtion is complete, use this command `python3 analyze.py` to generate a text file with all variance data, plots of all covariance matrices and sensitivity vectors, and a plot with top contributors to the total response variance.


# How to use TOFFEE for KCODE evaluations
1. Add your MCNP file for to the folder. The file `ksen_test.inp` is an example input file. The only addition to your input file will be two comment lines that identify the beginning and ending of your materials. The comments are `c Materials` and `c End Materials`, respectively.
2. In the file `submit.py`, you can edit the name of the input file in line 21.
3. Follow the proccess from the test problem from steps 5-7 to evaluate your MCNP file. 

# How to use TOFFEE for SDEF evaluations
1.  Add your MCNP file for to the folder. The file `ksen_test.inp` is an example input file. The only addition to your input file will be two comment lines that identify the beginning and ending of your materials. The comments are `c Materials` and `c End Materials`, respectively.
2. In the file `submit.py`, you can edit the name of the input file in line 21.
3. `TOFFEE.py` will need to be altered to specify the tally identifiers. `TOFFEE.py` is set to read an F4 tally that is on cell 1 without a multiplier. All of these parameters can be specified to find any tally that is affected by PERT cards. TOFFEE looks for the "1tally" output within the MCNP output files. The cell identifier (and multiplier for reaction rates) are used to identifiy when TOFFEE can read the sensitivity coefficients. The variables `main.tally_identifier`, `main.location_identifier`, `main.multiplier_identifier`, and `main.tally_multiplier` need to be altered in the file `submit.py` to be for your specific tally.
4. Use this command `python3 submit.py` to evaluate the covariance library and input file.
5. Evaluate all files in the `working_dir` directory.
6. Use this command `python3 analyze.py` to analyze the data.

# Available options for data analysis
3. If you would like a different energy discretization ([Default is SCALE 56 group structure](https://scale-manual.ornl.gov/XSLib.html#the-56-group-library)), change the vector in line 25 in the file `submit.py`.
4. If you would like a different list of reaction, change the vectors in line 33 and 36 in the file `submit.py`.
5. The plotting options can be turned off in lines 25, 27, 32, and 35 in the file `analyze.py` for the plots you would like the code to make. 
