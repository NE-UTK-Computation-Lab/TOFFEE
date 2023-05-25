# TOol For Fast Error Estimation (TOFFEE)
TOFFEE is a python script that automates cross-section uncertainty analysis using the tool available in MCNP6.2 and NJOY. The program utilizes the PERT and KSEN methods to generate sensitivities, uses NJOY to generate covariance data, and uses NUMPY to do the sandwich rule to calculate uncertainites.

# How to install and use TOFFEE
1. Download TOFFEE.py, element_letters.csv, ERRORR_reader, and njoy_template.txt into a folder.
2. The programs `submit_to_cluster`(line 265 in TOFFEE.py) and `submit_njoy_run`(line 114 in ERRORR_reader.py) must be altered to submit the mcnp and njoy runs for your specific computation machine.
3. Alter `self.endf_folder` in line 23 of ERRORR_reader to the ENDF6 formated neutron cross section library of your choise
4. Add your MCNP file for to the folder.
5. Change the name of the file "base_sdef_temp.txt" to the name of your MCNP file in lines 1334 and 1343.
6. If your MCNP file is source driven, the script is ready to run.

To run a kcode file in TOFFEE.py, change 'pert' to 'ksen' in line 1334. In line 1343 change `submit_pert` to `submit_ksen` and change `covariance,0.01,True` to `covariance,True`
