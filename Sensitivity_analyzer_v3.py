import os
import time
import copy
import numpy as np
import csv
#########Inputs/setting print text
# This is the list of reactions as they are named
#list_of_rxn=['rxn 1', 'rxn 2', 'rxn 4', 'rxn 18', 'rxn 102' ]
# This is the name given to each reaction that is called. This list should mirror the following list
reaction_names=['Total', 'Elastic', 'Total inelastic', 'Total Fission','(n-gamma)']
# This is a list of the reactions called in a format MCNP can read. The number is the MT number of the reaction desired.
reaction_calls=['rxn=1', 'rxn=2', 'rxn=4', 'rxn=18', 'rxn=102' ]

#This list is a set of naming pairs for the elements present elements found in MCNP material cards and a common terminology for the element
material_names=[['1001.70c','H1'],['8016.70c','O16'],['92234.70c','U234'],['92235.70c','U235'],['92236.70c','U236'],['92238.70c','U238']]

#energy_bins=[1.00E-10,5.00E-10,7.50E-10,0.000000001,1.20E-09,1.50E-09,0.000000002,2.50E-09,0.000000003,0.000000004,0.000000005,7.50E-09,0.00000001,2.53E-08,0.00000003,0.00000004,0.00000005,0.00000006,0.00000007,0.00000008,0.00000009,0.0000001,0.000000125,0.00000015,0.000000175,0.0000002,0.000000225,0.00000025,0.000000275,0.0000003,0.000000325,0.00000035,0.000000375,0.0000004,0.00000045,0.0000005,0.00000055,0.0000006,0.000000625,0.00000065,0.0000007,0.00000075,0.0000008,0.00000085,0.0000009,0.000000925,0.00000095,0.000000975,0.000001,0.00000101,0.00000102,0.00000103,0.00000104,0.00000105,0.00000106,0.00000107,0.00000108,0.00000109,0.0000011,0.00000111,0.00000112,0.00000113,0.00000114,0.00000115,0.000001175,0.0000012,0.000001225,0.00000125,0.0000013,0.00000135,0.0000014,0.00000145,0.0000015,0.00000159,0.00000168,0.00000177,0.00000186,0.00000194,0.000002,0.00000212,0.00000221,0.0000023,0.00000238,0.00000247,0.00000257,0.00000267,0.00000277,0.00000287,0.00000297,0.000003,0.0000031,0.0000032,0.0000035,0.00000373,0.0000041,0.0000047,0.000005,0.0000054,0.000006,0.00000625,0.0000065,0.00000675,0.000006875,0.000007,0.00000715,0.0000081,0.0000091,0.00001,0.0000115,0.0000119,0.0000129,0.0000144,0.000016,0.000017,0.0000185,0.0000194,0.00002,0.0000205,0.0000212,0.00002175,0.0000225,0.000025,0.0000275,0.00003,0.00003125,0.00003175,0.00003325,0.00003375,0.000035,0.0000355,0.000036,0.000037,0.00003713,0.00003727,0.00003763,0.000038,0.0000391,0.0000396,0.000041,0.0000424,0.000044,0.0000452,0.0000483,0.0000506,0.0000534,0.000058,0.000061,0.000063,0.000065,0.0000675,0.000072,0.000076,0.00008,0.0000817,0.00009,0.000097,0.0001012,0.000105,0.000108,0.000113,0.000116,0.0001175,0.000119,0.000122,0.000143,0.00017,0.00018,0.0001877,0.0001885,0.0001915,0.000193,0.000202,0.0002074,0.0002095,0.00022,0.00024,0.000285,0.000305,0.00055,0.00067,0.000683,0.00095,0.00115,0.0015,0.00155,0.0018,0.0022,0.00225,0.0025,0.003,0.00374,0.0039,0.0057,0.00803,0.0095,0.013,0.017,0.02,0.03,0.045,0.05,0.052,0.06,0.073,0.075,0.082,0.085,0.1,0.1283,0.149,0.2,0.27,0.33,0.4,0.42,0.44,0.47,0.492,0.55,0.573,0.6,0.67,0.679,0.75,0.82,0.8611,0.875,0.9,0.92,1.01,1.1,1.2,1.25,1.317,1.356,1.4,1.5,1.85,2.354,2.479,3,4.304,4.8,6.434,8.187,10,12.84,13.84,14.55,15.68,17.33,20]
# This is the value of the energy bin in MeV
#These energy bins are the ones used in the covariance matrix used to find the variance due to the cross section
energy_bins_temp=[1.00E-5,4.00E-03,1.00E-02,2.53E-02,4.00E-02,5.00E-02,6.00E-02,8.00E-02,1.00E-01,1.50E-01,2.00E-01,2.50E-01,3.25E-01,3.50E-01,3.75E-01,4.50E-01,6.25E-01,1.01E+00,1.08E+00,1.13E+00,5.00E+00,6.25E+00,6.50E+00,6.88E+00,7.00E+00,2.05E+01,2.12E+01,2.18E+01,3.60E+01,3.71E+01,6.50E+01,6.75E+01,1.01E+02,1.05E+02,1.16E+02,1.18E+02,1.88E+02,1.92E+02,2.25E+03,3.74E+03,1.70E+04,2.00E+04,5.00E+04,2.00E+05,2.70E+05,3.30E+05,4.70E+05,6.00E+05,7.50E+05,8.61E+05,1.20E+06,1.50E+06,1.85E+06,3.00E+06,4.30E+06,6.43E+06,2.00E+07]
#Converts energy bin values to eV
energy_bins=[round(a/10**6,12) for a in energy_bins_temp]





##########create perturbation strings
##inputs
#template_name is the name of the mcnp file the sensitivity analysis will be used on
#density_change is the value of percent of macroscopic cross section perturbed. A good value to use is 0.01. This is used for ksen problems, but still should be given an arbitrary value.
#reactions is a list of reactions that a preturbation is desired. This refers to the variable reaction_calls used created in the input section
#erg_bins is a list of the energy bin system used for the perturbations. This variable is energy_bins, as seen in the inputs section, for my case.
#run_type is a variable that describes the type of pertuebation cards that will be created. Three options are currently available:'kpert','ksen', and 'pert' 
def make_pert_cards(template_name,density_change,reactions,erg_bins,run_type):
    #opening the mcnp file to read the materials and elements of each material
    template=open(template_name,'r')
    #creating variables to build later/ setting defaults
    list_of_materials=[]
    elements_out=[]
    in_materials=False
    j=0
    #####This loop will read the MCNP input file to create a list of materials and elements of that material
    #This loop steps through the MCNP input line by line
    for line in template:
        #This function creates a list of the elements of the current line
        line_split=line.split()
        #This line looks for an identifier within the MCNP code the tells where the materials list begins
        if "c Materials" in line:
            in_materials=True
            continue
        # this determines if we are in the list of materials
        # if we are, then we begin to create a list of the materials within the system.
        if in_materials==True:
            #This if is true if we are at the first element of a material. This initiallizes the temparary variables:
            #list_of_elements-the element names in MCNP form. EXP:1001.70c is Hydrogen-1 with the 70c cross section library values used
            #element_material- the atom fraction given in the MCNP input material card. The fraction must be in the atom fraction form.    
            if list(line_split[0])[0]=='m':
                #initialization of the first material card within the system
                if j==0:
                    list_of_elements=[]
                    element_material=[]
                    # the number of the m card given i.e.M21 in MCNP
                    mat_num=line_split[0]
                    list_of_elements.append(line_split[1])
                    element_material.append(line_split[2])
                    j+=1
                    #All other materials in the system start with this section
                else:
                    mat_num_list=list(mat_num)
                    temp_num=''
                    for j in range(len(mat_num_list)):
                        if j>0:
                            temp_num+=str(mat_num_list[j])
                    temp_list=[temp_num,list_of_elements,element_material]
                    list_of_materials.append(temp_list)
                    list_of_elements=[]
                    element_material=[]
                    mat_num=line_split[0]
                    list_of_elements.append(line_split[1])
                    element_material.append(line_split[2])
                    j+=1
            #This line determines if we are at the end of the materials list to end the building of the material library            
            elif'c End Materials' in line:
                mat_num_list=list(mat_num)
                temp_num=''
                for i in range(len(mat_num_list)):
                    if i>0:
                        temp_num+=str(mat_num_list[i])
                temp_list=[temp_num,list_of_elements,element_material]
                list_of_materials.append(temp_list)
                list_of_elements=[]
                in_materials=False
            #This line skips commented lines
            elif line_split[0]=='c':
                continue
            #This line adds all additional materials within the system. This condition is most common once in materials
            else:
                list_of_elements.append(line_split[0])
                element_material.append(line_split[1])
    template.close()
    #####This section makes a list of the cells where each of the materials are located. It also saves the normal density of the material
    #This section opens the MCNP input
    template=open(template_name,'r')
    in_cells=False
    cell_materials=[]
    cell_materials_check=[]
    #initializes elements of a list for each material within the system
    for materials in list_of_materials:
        cell_materials.append('')
        cell_materials_check.append(False)
    #adds the cell numbers that contain each material
    for line in template:
        # This looks for a call in the MCNP template that tells when the cells section starts
        if 'c Cells' in line:
            in_cells= True
            continue
        # Once in the cells the loop looks for different material names to determine if the cell has material x in it
        if in_cells== True:
            line_split=line.split()
            for i in range(len(list_of_materials)):
                if line_split[1]==list_of_materials[i][0]:
                    if cell_materials_check[i]==False:
                        cell_materials[i]+=line_split[0]
                        cell_materials_check[i]=True
                        list_of_materials[i].append(line_split[2])
                    else:
                        cell_materials[i]+=','+line_split[0]
        # This is a call that looks for a line that determines the end of the cells section of the MCNP input file
        if 'c End Cells'in line:
            in_cells=False
    #####This section creates the new added materials with a perturbation in each isotope. Only used for PERT and KPERT runs.
    list_of_new_materials=[]
    for x in range(len(list_of_materials)):
        list_temp=''
        element_temp=[]
        for i in range(len(list_of_materials[x][1])):
            for j in range(len(list_of_materials[x][1])):
                #Creates the first line of the materials that do not have the perturbation on this element
                if j==0 and j!=i:
                    list_temp+='m'+str(list_of_materials[x][0])+str(i+1)+'    '+list_of_materials[x][1][j]+' '+list_of_materials[x][2][j]+'\n'
                #Creates the first line of the materials that do have the perturbation on this element
                elif j==0 and j==i:
                    list_temp+='m'+str(list_of_materials[x][0])+str(i+1)+'    '+list_of_materials[x][1][j]+' '+str(round(float(list_of_materials[x][2][j])*(1+float(density_change)),9))+'\n'
                    element_temp.append(list_of_materials[x][1][j])
                #Creates the other lines of the materials that do the perturbation on this element
                elif j!=0 and j==i:
                    list_temp+='       '+list_of_materials[x][1][j]+' '+str(round(float(list_of_materials[x][2][j])*(1+float(density_change)),9))+'\n'
                    element_temp.append(list_of_materials[x][1][j])
                #Creates the other lines of the materials that do not have the perturbation on this element
                else:
                    list_temp+='       '+list_of_materials[x][1][j]+' '+list_of_materials[x][2][j]+'\n'
        elements_out.append(element_temp)    
        list_of_new_materials.append(list_temp)
    #####This section creates the added MCNP calls that are used to find the sensititvities
    #This section is for creating the kpert cards need for keff sensitivities
    if run_type=='kpert':
        list_of_pert=[]
        # creates a set of cards for every material
        for x in range(len(list_of_materials)):
            temp_list_pert_mat=[]
            # Every isotope within material x
            for i in range(len(list_of_materials[x][1])):
                temp_list_pert=[]
                # Every reaction of interest. Note: The reactions chosen appear for all isotopes even if they do not have that particular reaction. EXAMPlE: If fission reaction is desired for the system with water and uranium, a H-1 fission sensitivity will be calculated, but will return void.
                for j in range(len(reactions)):
                    temp=[]
                    # a kpert card is need for each energy bin range
                    # cell refers to the cells of interest for the material being perturbed, rho is the perturbed desity value, mat refers to the name of the new material used for perturbation, erg is the energy bin lower and upper bounds for a particular energy bin 
                    for k in range(len(erg_bins)-1):
                        pert_string='kpert'+str(k+1)+' cell='+cell_materials[x]+' rho='+str(round(float(list_of_materials[x][3])+((float(list_of_materials[x][2][i])*(1+float(density_change))-float(list_of_materials[x][2][i]))),9))+' mat='+str(list_of_materials[x][0])+str(i+1)+' '+ reactions[j]+' erg= '+str(erg_bins[k])+' '+str(erg_bins[k+1])+'\n'
                        temp.append(pert_string)
                    temp_list_pert.append(temp)
                temp_list_pert_mat.append(temp_list_pert)    
            list_of_pert.append(temp_list_pert_mat)
    #This section is for creating the pert cards need for source driven sensitivities
    #very similar to kpert above for source driven problems
    elif run_type=='pert' :
        list_of_pert=[]
        for x in range(len(list_of_materials)):
            temp_list_pert_mat=[]
            
            for i in range(len(list_of_materials[x][1])):
                temp_list_pert=[]
                for j in range(len(reactions)):
                    temp=[]
                    # the method=2 ensures the first order sensitivities are calculated
                    for k in range(len(erg_bins)-1):
                         pert_string='pert'+str(k+1)+':n cell='+cell_materials[x]+' rho='+str(round(float(list_of_materials[x][3])+((float(list_of_materials[x][2][i])*(1+float(density_change))-float(list_of_materials[x][2][i]))),9))+' mat='+str(list_of_materials[x][0])+str(i+1)+' '+ reactions[j]+' method=2'+' erg= '+str(erg_bins[k])+' '+str(erg_bins[k+1])+'\n'
                         temp.append(pert_string)
                    temp_list_pert.append(temp)
                temp_list_pert_mat.append(temp_list_pert)    
            list_of_pert.append(temp_list_pert_mat)
    #This section creates cards for ksen keff sensitivities
    elif run_type=='ksen':
       list_of_pert=[]
       for x in range(len(list_of_materials)):
            temp_list_pert_mat=[]
            for i in range(len(list_of_materials[x][1])):
                temp_list_pert=[]
                for j in range(len(reactions)):
                    temp=[]
                    #XS denotes cross section sensitiviites, iso is for the isotope name (1001 for H-1)
                    for k in range(len(erg_bins)-1):
                         pert_string='ksen'+str(k+1)+' XS'+' ISO= ' +str(list_of_materials[x][1][i])+' '+'cell='+cell_materials[x]+' '+str(reactions[j]) +' erg= '+str(erg_bins[k])+' '+str(erg_bins[k+1])+'\n'
                         temp.append(pert_string)
                    temp_list_pert.append(temp)
                temp_list_pert_mat.append(temp_list_pert)    
            list_of_pert.append(temp_list_pert_mat)
       
    else:
        #This is an error report for an improper run name.
        print('No run type selected. Please use ksen for keff, kpert for keff, or pert for sdef.')
                        
    return list_of_pert,list_of_new_materials,elements_out
###output
#list_of_pert- a list containing the perturbation cards to be added to an mcnp input file
#list_of_new_materials- A list containing text strings for each new material that is used in perturbation cards
#elements_out- a list of lists where the outer list represents the base materials in the system and the inner list is the elements of that material in MCNP format. 



########## This function creates submission scripts and submits MCNP jobs to the UTK NE cluster.
###inputs
#job_name- the naming string for the given run. Example:if you wanted to run job_in.txt the input would be 'job' 
def submit_to_cluster(job_name):
    # The folder this is called in must contain a mcnp input titled job_name_in.txt
    # This line writes the script text to submit the job
    cluster_input_string=' #!/bin/bash \n #PBS -V \n #PBS -q corei7 \n #PBS -l nodes=1:ppn=8 \n hostname \n module load MCNP6/2.0 \n RTP="/tmp/runtp--".`date "+%R%N"` \n cd $PBS_O_WORKDIR \n mcnp6 TASKS 8 I=%%%INPUT%%%_in.txt O=%%%INPUT%%%_out.txt runtpe=$RTP \n grep -a "final result" %%%INPUT%%%_out.txt > %%%INPUT%%%_done.dat \n rm $RTP'
    # This replace the %%%input%%% text with the given job_name 
    script_string = cluster_input_string.replace("%%%INPUT%%%", job_name)
    # This line creates the name for the submission script for the cluster
    script_file_string = job_name + "_script.txt"
    #These lines write the scriptfile given with script_string
    script_file = open(script_file_string, 'w')
    script_file.write(script_string)
    script_file.close()
    #Waiting for 5 seconds due to issues where the cluster can take a few seconds to recognize the script file
    time.sleep(5)
    
    ##This line looks to see if the output file already exists and removes the file if it does
    # This is because MCNP will name the job with the next alphabet letter ie job_ouu.txt instead of job_out.txt
    if os.path.exists(job_name+'_out.txt'):
        os.system('rm '+job_name+'_out.txt')
        print('file removed out')
        #waiting to verify the cluster knows the file is no longer present
        time.sleep(25)
        
    # This line submits the job via qsub command
    # the part of the string '>> cluster_submissions.txt' will create an output file with the text printed to the command window
    current_Dir = os.getcwd()
    os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir + ' && qsub ' + script_file_string+ '>> cluster_submissions.txt'+'"')
    output_string=job_name+'_out.txt'
    #this loop looks for a file created at the end of a MCNP run
    while not os.path.exists(job_name+'_done.dat'):
        template=open('cluster_submissions.txt','r')
        #This reads through the text file created above to see if an error has occured
        for line in template:
            if 'socket_connect_unix failed' in line:
                #if this line is found the job is resubmitted
                os.remove('cluster_submissions.txt')
                time.sleep(5)
                os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir + ' && qsub ' + script_file_string+ '>> cluster_submissions.txt'+'"')
            else:
                #if not then the code waits for the output file
                print('Job not complete waiting 15 seconds')
                time.sleep(15)
        template.close()
    #This shows the job is finished
    print('file found')
    #waits some time to ensure files are finished updating
    time.sleep(15)
    #removes job_done.dat and cluster_submissions files
    os.remove(job_name+'_done.dat')
    os.remove('cluster_submissions.txt')
    print('file removed')
    #deletes temp MCNP file
    if os.path.exists('srctp'):
        os.remove('srctp')
    return output_string
###output
#The output of this function is the name of the output file created with the cluster submission





#######This function pulls the data for a covariance matrix from a csv file and creates a matrix that can be used by python
###inputs
#element- a text string that represents the MCNP style string of an element with the library specified. ie 1001.70c is valid, 1001 is not currently
#reaction- a texte string that meets the format needed for the perturbation card style of formating. ie rxn=1 is valid, 1 is not
def pull_covariance_matrix(element,reaction):
    #removes the .xxc portion of the isotope name    
    element_call=element[:-4]
    #removes the rxn= from the reaction input
    rxn=reaction[4:]
    #looks for a file named element call_reaction.csv in the folder covariance_library
    #using the example inputs: covariance_library/1001_1.csv
    matrix_file=open('covariance_library/' +str(element_call)+'_'+str(rxn)+'.csv',encoding='utf-8-sig')
    cov_matrix_temp=[]
    #reads the csv file using the csv library
    file_mat=csv.reader(matrix_file)
    
    #this line reads the lines of the csv and creates a list for every row
    for row in file_mat:
        temp=[]
        for i in row:
            temp.append(float(i))
            
        cov_matrix_temp.append(temp)
    
    #converts a list of lists to an array for matrix multiplication when the covariance matrix is used                    
    cov_matrix=np.array(cov_matrix_temp)
    
    return cov_matrix
###output
#cov_matrix- An numpy matrix that contains the elements of the given covariance matrix






##### This function reads the output file of a kpert run of a keff eignevalue run for a sensitivity analysis
##### then does a sandwich rule with the covariance to report the variance of the system
###inputs 
#out_file- a string for the name of the output file from a MCNP run with kpert cards
#energy_b- list containing the energy bin structure used in the analysis
#element- string containing the isotope name for this kpert run of MCNP perturbation,ie 1001.70c
#reaction- string containg the reaction used for this perturbation in KPERT language, ie rxn=1
#perturb- a float that represents the value of perturbation for the system, if 1.01 is perturbed then perturb= 0.01
def read_kpert_code(out_file,energy_b,element,reaction,perturb):
    #prints the name of the output file name and sets it to a variable
    print(out_file)
    file_output=open(out_file,'r')
    in_tally = False
    in_data = False

    i=0
    vector=[]
    in_data=False
    in_tally=False
    # reading through the output file while we look for particular lines in the text.
    for line in file_output:
        # this signifies we have entered the kpert results in the MCNP file
        if 'kpert          delta-rho' in line:
            in_tally= True
            continue
            
            
        if in_tally==True:
            #there is a space after the first check so we do this to make sure we are in a data line
            if'1' in line:
                in_data=True

            if in_data==True:
                #splits the text into individual values
                line_split=line.split()
                #the output line_split[1] represents delta k so we divide by the perturbation value
                #This builds the sensitivity vector as a list
                vector.append([float(line_split[1])/float(perturb)])
                i+=1            
                #this check will determine if we exceed the number of data points, then ends the data scanning
                if i<len(energy_b)-1:
                    continue
                else:
                    in_data=False
                    in_tally=False
                        
    file_output.close()
    # This converts the sensitivity vector to an numpy array for matrix operations
    sens_vector=np.array(vector)
    # Transposing the sensitivity vector
    sens_vect_t=np.transpose(sens_vector)
    # Pulls the covariance matrix using the function pull_covariance_matrix
    cov_mat=pull_covariance_matrix(element,reaction)
    #Doing the first multiplication of the S^t*C*S sandwich rule
    temp_vect=np.matmul(sens_vect_t,cov_mat)
    #Second multiplication of the sandwich rule
    uncert=np.matmul(temp_vect,sens_vector)
    # creating a list style output for later use
    uncertainty=([uncert])            
            
    return uncertainty
###output
#uncertainty- a list that contains the variance (using sandwich rule) of the system based on the isotope and reaction specified





##### This function reads the output of a keff calculation to pull sensitivities and
##### then uses the sandwich rule to report the variance of the system
###input
#out_file- the name of the output file from MCNP
#energy_b- list containing the energy bin structure used in the analysis
#element- string containing the isotope name for this kpert run of MCNP perturbation,ie 1001.70c
#reaction- string containg the reaction used for this perturbation in KPERT language, ie rxn=1
def read_ksen_code(out_file,energy_b,element,reaction):
    #This outputs the name of the file being read
    print(out_file)
    #creating a variable the represents the text from the output file
    file_output=open(out_file,'r')
    in_tally = False
    in_data = False


    vector=[]
    in_data=False
    in_tally=False
    # this loop reads through the output file 
    for j in range(len(energy_b)-1):
        file_output.close()
        file_output=open(out_file,'r')
        for line in file_output:
            # This looks for the string the represents tha start of sensitivity data for the current energy bin
            if 'sensitivity profile    ' in line and str(j+1) in line:
                in_tally= True
                continue
            
            
            if in_tally==True:
                # once the sensitivity profile of interest is found, we find where the data starts with this line
                if'isotope       reaction       sensitivity   rel. unc.' in line:
                    in_data=True
                    continue
                
                if in_data==True:
                    # this verifies we are in the numbers of the data
                    if '.' in line:  
                        line_split=line.split()
                        # the third element of the line contains the sensitivity, dk/dsig
                        vector.append(float(line_split[2]))          
                        in_data=False
                        in_tally=False
                        break
                        
    file_output.close()
    # creating the sensitivity vector in energy bins, row vector
    sens_vector=np.array(vector)
    # creating a transpose version of the sensitivity vector, column vector
    sens_vect_t=np.transpose(sens_vector)
    # Pulling the covariance matrix from data
    cov_mat=pull_covariance_matrix(element,reaction)
    # First half of the matrix multiplicatons for the SCS sandwich rule
    temp_vect=np.matmul(sens_vect_t,cov_mat)
    # second multiplication
    uncert=np.matmul(temp_vect,sens_vector)
    # report uncertainty as a matrix value
    uncertainty=([uncert])            
            
    return uncertainty
###output
#uncertainty- a list that contains the variance (using sandwich rule) of the system based on the isotope and reaction specified





##### This function reads the output of a keff calculation to pull sensitivities and
##### then uses the sandwich rule to report the variance of the system
###input
#out_file- the name of the output file from MCNP
#energy_b- list containing the energy bin structure used in the analysis
#element- string containing the isotope name for this kpert run of MCNP perturbation,ie 1001.70c
#reaction- string containg the reaction used for this perturbation in KPERT language, ie rxn=1
#perturb- The value used for perturbation, ie if 1.001 if perturbed, perturb=0.001
def read_sdef_code(out_file,element,reaction,energy_b,perturb):
    #Initializing temp files
    tally_base=[0,0]
    # prints the file being used for data analysis
    print(out_file)
    #creates a variable for the output file text
    file_output=open(out_file,'r')
    in_tally = False
    in_data = False
    flux=0
    std=0
    # reads the text fir the file
    for line in file_output:
        #looks for the first f14 tally of the MCNP output file
        if '1tally       14'  in line:
            in_tally = True
            
        if in_tally:
            # cell 1 is the cell the f14 tally is loctaed in
            #this line appears directly above the data line
            if 'cell  1' in line:
                in_data=True
                continue
            if in_data:
                in_tally=False
                line_split=line.split()
                # This pulls the unperturbed f tally value from the MCNP output file
                # This value is needed to calculatate the sensitivity
                flux=float(line_split[0])
                std=float(line_split[1])
                break
                
    file_output.close()
    in_tally=False
    in_data=False

    tally_base[0]=flux
    tally_base[1]=std
    tallies=[tally_base]
    
    # This loop reads the the perturbed f tally values and calculates the sensitivities
    tallies.append([])
    vector=[]
    # This sweeps through the energy bins
    for j in range(len(energy_b)-1):
        file_output.close()
        file_output=open(out_file,'r')

        for line in file_output:
            # This is the call that denotes the current perturbation related to the current energy bin
            if 'the following output gives the predicted change in a tally for perturbation' in line and ' '+str(j+1)+'.' in line:
                in_tally= True
                continue     
            # This denotes the cell we are taking the tally in. This could change based on the type of tally chosen to inspect
            if in_tally==True:
                if'cell  1' in line:
                    in_data=True
                    continue
                #This line does that math required to calculate the sensitivity value
                if in_data==True:
                    line_split=line.split()
                    # change in f tally
                    c_1=(float(line_split[0]))
                    # original f tally
                    c_0=float(tallies[0][0])
                    # creating/ appending sensitivity vector
                    vector.append([c_1/c_0/perturb])
                    in_data=False
                    in_tally=False
    file_output.close()
    # Converts the sensitivity vector from a list to an array in order to matrix multiply
    sens_vector=np.array(vector)
    # Transposing the sensitivity vector
    sens_vect_t=np.transpose(sens_vector)
    # pulling the covariance matrix from a library
    cov_mat=pull_covariance_matrix(element,reaction)
    # First multiplication of the sandwich rule SCS
    temp_vect=np.matmul(sens_vect_t,cov_mat)
    #second multiplication of the sandwich rule
    uncert=np.matmul(temp_vect,sens_vector)
    #Making the output a list
    uncertainty=([uncert])
    return uncertainty
###output
#uncertainty- a list that contains the variance (using sandwich rule) of the system based on the isotope and reaction specified





##### This function takes the perturbation cards created with the make_pert_cards function and executes them, then reads
##### the data with the particular read function. This function also creates a csv file for easy data analysis, as well as
##### a dictonary containing the same data
### Inputs
#reaction_list- This is the list of reactions that the run will analyze(see reaction_calls in initializations)
#pert_calls- This is the output from the make_pert_cards function called list_of_pert
#new_materials- This is the output from the make_pert_cards function called list_of_new_materials
#names- This is a list of typical names given to the MT card reactions. This should mirror reaction_list such as:
    #reaction_list=['rxn=1']    names=['Total']. The naming does not matter and is used for making tables make more sense
#element_list- This is a list of list where the elements of the outer list represent the material. Each material is a list 
    #containg the isotopes present in said material
#template_file- This is the file name the sensitivity analysis will be done on.
#energies- A list containing the energy bin structure for the system
#density_change- The value of perturbation used
#sub_job- Given as True, where True means to submit MCNP jobs and False means to read data from files that have already ran.
    #This is useful for debugging
def submit_ksen(reaction_list,pert_calls,new_materials,names,element_list,template_file,energies,density_change,sub_job=True):
    # These lines are used for making the CSV file have columns that line up properly
    top_string=','
    top_string_two=','
    top_string_three='Reactions,'
    table=[]
    # This loop appends the csv file with the materials being used within the system
    # This loop also initializes the dictionary of uncertainty values
    for i in range(len(pert_calls)):
        # notes the material we the isotopes are in
        top_string+='Material'+str(i+1)
        # this loop tells cyles throught the isotopes used
        for j in range(len(pert_calls[i])):
            #pulls the isotope name from the fisrt pert call
            line_split=pert_calls[i][j][0][0].split()
            word=line_split[3]
            #this removes the library name of the material is one is present
            if '.' in word:
                element_string=word[:-4]
############# create a list of names to call from and do a if x=y then add isotope name here
            top_string_two+=element_string+',,'
###############
            # This section formats the text in the table
            top_string+=',,'
            top_string_three+='keff,,'
            # Initializing the dictionary
            temp_dict={'Material':str(i+1),'Element':element_string}
            table.append(temp_dict)
            
    #This section opens a csv file, writes the header information, then closes the csv file
    output_file=open("data_output_ksen.csv", 'w')
    output_file.write(top_string+'\n'+top_string_two+'\n'+top_string_three+'\n')
    output_file.close()
    
    
    #This loop is the major section of this function
    # This loop is where the MCNP files are created, then executes the files on the utk ne cluster through submit_to_cluster
    # After the results are read with the read_ksen_code function
    #### This loop is under construction where I may add parallelization to this section
    data=[]
    vect=[]
    m=1
    for i in range(len(pert_calls)):
        data.append([])
        for j in range(len(pert_calls[i])):
            #####multiprocess this section first
            for k in range(len(pert_calls[i][j])):
                job_num='job_'+str(i)+'_'+str(j)+'_'+str(k)
                job=open(job_num+'_in.txt','w')
                write_string=''
                for l in range(len(pert_calls[i][j][k])):
                    write_string+=pert_calls[i][j][k][l]

                template=open(template_file,'r')
                for line in template:
                    job.write(line)
                    if 'c End Materials' in line:
                        job.write('c PERT Materials'+'\n')
                        for material in new_materials:
                            job.write(material)
                        job.write('c End PERT Materials'+'\n')    
                template.close()
                job.write(write_string)
                job.close()
                if sub_job==True:
                    out=submit_to_cluster(job_num)
                else:
                    out=job_num+'_out.txt'
                
                react=reaction_list[k]
                element_cov=element_list[i][j]
                vect.append(read_ksen_code(out,energies,element_cov,react))
                if m==5:
                    data[i].append(vect)
                    vect=[]
                    m=0
                m=m+1
                
    list_of_data=[]
    for i in range(len(data[0][0])-1):
        list_of_data.append([])
        
    list_of_data.append([])
    list_of_data[0]=','
    for i in range(len(data)):
        #each material
        for j in range(len(data[i])):
            #each element for that material
            for k in range(len(reaction_calls)):

                if i==0 and j==0:
                    list_of_data[k]=(str(names[k])+','+str(float(data[i][j][k][0]))+',,')
                else:
                    list_of_data[k]+=str(float(data[i][j][k][0]))+',,'
    
    for i in range(len(list_of_data)):
        output_file=open("data_output_k.csv", 'a')
        output_file.write(list_of_data[i]+'\n')
        output_file.close()
            
    number_of_elements=0
    for i in range(len(data)):
        
        for j in range(len(data[i])):
            table[number_of_elements]['Rxn']=[]
            for k in range(len(data[i][j])):
                table[number_of_elements]['Rxn'].append({'Name':names[k],'Uncertainty':(float(data[i][j][k][0][0]))})
                    
            number_of_elements+=1
    
    full_table=copy.deepcopy(table)
    for out in table:
        removal_list=[]
        for react in out['Rxn']:
            if react['Uncertainty']==0:
                removal_list.append(react)
                
        for elements in removal_list:
            out['Rxn'].remove(elements)
    
    
            
    return table,full_table

# def test_multi(pert_calls,template_file,new_materials,sub_job):
#     data=[]
#     vect=[]
#     m=1
#     for i in range(len(pert_calls)):
#         data.append([])
#         for j in range(len(pert_calls[i])):
#             for k in range(len(pert_calls[i][j])):
#                 job_num='job_'+str(i)+'_'+str(j)+'_'+str(k)
#                 job=open(job_num+'_in.txt','w')
#                 write_string=''
#                 for l in range(len(pert_calls[i][j][k])):
#                     write_string+=pert_calls[i][j][k][l]

#                 template=open(template_file,'r')
#                 for line in template:
#                     job.write(line)
#                     if 'c End Materials' in line:
#                         job.write('c PERT Materials'+'\n')
#                         for material in new_materials:
#                             job.write(material)
#                         job.write('c End PERT Materials'+'\n')    
#                 template.close()
#                 job.write(write_string)
#                 job.close()
#                 if sub_job==True:
#                     out=submit_to_cluster(job_num)
#                 else:
#                     out=job_num+'_out.txt'
                
#                 react=reaction_list[k]
#                 element_cov=element_list[i][j]
#                 vect.append(read_ksen_code(out,reaction_list,energies,element_cov,react,density_change))
#                 if m==5:
#                     data[i].append(vect)
#                     vect=[]
#                     m=0
#                 m=m+1


##### This function takes the perturbation cards created with the make_pert_cards function and executes them, then reads
##### the data with the particular read function. This function also creates a csv file for easy data analysis, as well as
##### a dictonary containing the same data
### Inputs
#reaction_list- This is the list of reactions that the run will analyze(see reaction_calls in initializations)
#pert_calls- This is the output from the make_pert_cards function called list_of_pert
#new_materials- This is the output from the make_pert_cards function called list_of_new_materials
#names- This is a list of typical names given to the MT card reactions. This should mirror reaction_list such as:
    #reaction_list=['rxn=1']    names=['Total']. The naming does not matter and is used for making tables make more sense
#element_list- This is a list of list where the elements of the outer list represent the material. Each material is a list 
    #containg the isotopes present in said material
#template_file- This is the file name the sensitivity analysis will be done on.
#energies- A list containing the energy bin structure for the system
#density_change- The value of perturbation used
#sub_job- Given as True, where True means to submit MCNP jobs and False means to read data from files that have already ran.
    #This is useful for debugging
def submit_kpert(reaction_list,pert_calls,new_materials,names,element_list,template_file,energies,density_change,sub_job=True):
    # These lines are used for making the CSV file have columns that line up properly
    top_string=','
    top_string_two=','
    top_string_three='Reactions,'
    table=[]
    #
    for i in range(len(pert_calls)):
        top_string+='Material'+str(i+1)
        for j in range(len(pert_calls[i])):
            line_split=pert_calls[i][j][0][0].split()
            word=line_split[3]
            element_string=''
            for x in range(len(word)):
                if x>=4:
                    element_string+=word[x]
            top_string_two+=element_string+',,'
            top_string+=',,'
            top_string_three+='keff,,'
            temp_dict={'Material':str(i+1),'Element':element_string}
            table.append(temp_dict)
        
    output_file=open("data_output_k.csv", 'w')
    output_file.write(top_string+'\n'+top_string_two+'\n'+top_string_three+'\n')
    output_file.close()
    
    
    
    data=[]
    vect=[]
    m=1
    for i in range(len(pert_calls)):
        data.append([])
        for j in range(len(pert_calls[i])):
            for k in range(len(pert_calls[i][j])):
                job_num='job_'+str(i)+'_'+str(j)+'_'+str(k)
                job=open(job_num+'_in.txt','w')
                write_string=''
                for l in range(len(pert_calls[i][j][k])):
                    write_string+=pert_calls[i][j][k][l]

                template=open(template_file,'r')
                for line in template:
                    job.write(line)
                    if 'c End Materials' in line:
                        job.write('c PERT Materials'+'\n')
                        for material in new_materials:
                            job.write(material)
                        job.write('c End PERT Materials'+'\n')    
                template.close()
                job.write(write_string)
                job.close()
                if sub_job==True:
                    out=submit_to_cluster(job_num)
                else:
                    out=job_num+'_out.txt'
                
                react=reaction_list[k]
                element_cov=element_list[i][j]
                vect.append(read_ksen_code(out,reaction_list,energies,element_cov,react,density_change))
                if m==5:
                    data[i].append(vect)
                    vect=[]
                    m=0
                m=m+1
                
    list_of_data=[]
    for i in range(len(data[0][0])-1):
        list_of_data.append([])
        
    list_of_data.append([])
    list_of_data[0]=','
    for i in range(len(data)):
        #each material
        for j in range(len(data[i])):
            #each element for that material
            for k in range(len(reaction_calls)):

                if i==0 and j==0:
                    list_of_data[k]=(str(names[k])+','+str(float(data[i][j][k][0]))+',,')
                else:
                    list_of_data[k]+=str(float(data[i][j][k][0]))+',,'
    
    for i in range(len(list_of_data)):
        output_file=open("data_output_k.csv", 'a')
        output_file.write(list_of_data[i]+'\n')
        output_file.close()
            
    number_of_elements=0
    for i in range(len(data)):
        
        for j in range(len(data[i])):
            table[number_of_elements]['Rxn']=[]
            for k in range(len(data[i][j])):
                table[number_of_elements]['Rxn'].append({'Name':names[k],'Uncertainty':(float(data[i][j][k][0][0]))})
                    
            number_of_elements+=1
    
    full_table=copy.deepcopy(table)
    for out in table:
        removal_list=[]
        for react in out['Rxn']:
            if react['Uncertainty']==0:
                removal_list.append(react)
                
        for elements in removal_list:
            out['Rxn'].remove(elements)
    
    
            
    return table,full_table




def submit_pert(reaction_list,pert_calls,new_materials,names,element_list,template_file,energies,density_change,sub_job=True):
    top_string=','
    top_string_two=','
    top_string_three='Reactions,'
    table=[]
    for i in range(len(pert_calls)):
        top_string+='Material'+str(i+1)
        for j in range(len(pert_calls[i])):
            line_split=pert_calls[i][j][0][0].split()
            word=line_split[3]
            element_string=''
            for x in range(len(word)):
                if x>=4:
                    element_string+=word[x]
            top_string_two+=element_string+',,'
            top_string+=',,'
            top_string_three+='uncertainty,,'
            temp_dict={'Material':str(i+1),'Element':element_string}
            table.append(temp_dict)
        
    output_file=open("data_output_sdef.csv", 'w')
    output_file.write(top_string+'\n'+top_string_two+'\n'+top_string_three+'\n')
    output_file.close()
    
    
    
    data=[]
    vect=[]
    m=1
    # This section builds the input decks 
    # This section will be where to make runs faster
    for i in range(len(pert_calls)):
        data.append([])
        for j in range(len(pert_calls[i])):
            for k in range(len(pert_calls[i][j])):
                job_num='job_'+str(i)+'_'+str(j)+'_'+str(k)
                job=open(job_num+'_in.txt','w')
                write_string=''
                for l in range(len(pert_calls[i][j][k])):
                    write_string+=pert_calls[i][j][k][l]

                template=open(template_file,'r')
                for line in template:
                    job.write(line)
                    if 'c End Materials' in line:
                        job.write('c PERT Materials'+'\n')
                        for material in new_materials:
                            job.write(material)
                        job.write('c End PERT Materials'+'\n')    
                template.close()
                job.write(write_string)
                job.close()
                if sub_job==True:
                    out=submit_to_cluster(job_num)
                else:
                    out=job_num+'_out.txt'
                    
                
                react=reaction_list[k]
                element_cov=element_list[i][j]
                vect.append(read_sdef_code(out,element_cov,react,energies,density_change))
                if m==5:
                    data[i].append(vect)
                    vect=[]
                    m=0
                m=m+1
            
    list_of_data=[]
    for i in range(len(data[0][0])):
        list_of_data.append([])
        
    list_of_data.append([])
    list_of_data[0]=','
    for i in range(len(data)):
        #each material
        for j in range(len(data[i])):
            #each element for that material
            list_of_data[0]+=str(float(data[i][j][0][0]))+',,'
            for k in range(len(reaction_calls)):

                if i==0 and j==0:
                    list_of_data[k+1]=(str(names[k])+','+str(float(data[i][j][k][0]))+',,')
                else:
                    list_of_data[k+1]+=str(float(data[i][j][k][0]))+',,'
    
    for i in range(len(list_of_data)):
        output_file=open("data_output_sdef.csv", 'a')
        output_file.write(list_of_data[i]+'\n')
        output_file.close()
            
    number_of_elements=0
    for i in range(len(data)):
        
        for j in range(len(data[i])):
            table[number_of_elements]['Rxn']=[]
            for k in range(len(data[i][j])):
                table[number_of_elements]['Rxn'].append({'Name':names[k],'Uncertainty':(float(data[i][j][k][0][0]))})
                    
            number_of_elements+=1
    
    full_table=copy.deepcopy(table)
    for out in table:
        removal_list=[]
        for react in out['Rxn']:
            if react['Uncertainty']==0:
                removal_list.append(react)
                
        for elements in removal_list:
            out['Rxn'].remove(elements)
    
            
    return table,full_table

    
    
[calls,material_update,element_list]=make_pert_cards('base_k_temp.txt',0.01,reaction_calls,energy_bins,'ksen')
table_output,full=submit_ksen(reaction_calls,calls,material_update,reaction_names,element_list,'base_k_temp.txt',energy_bins,0.01,False)

print('done')
