import os
import time
import numpy as np
import shutil
import pickle
import matplotlib.pyplot as plt
from ERRORR_reader import ERRORR_tools
plt.rcParams.update({'font.size': 16})

class TOFFEE_class:
    ##### Initialization for a simple test problem
    def __init__(self):
        ### These are things the user can change to their specific need
        # Using NJOY and provided ENDF6 formated covariance data
        self.covariance_method = 'NJOY'
        # runs the MCNP files based on the source file
        self.run_mcnp = True
        self.mcnp_file_name = 'test_case.inp'  
        # SCALE 56 group structure in MeV
        self.energy_bin_structure = [1.00E-11,4.00E-09,1.00E-08,2.53E-08,4.00E-08,5.00E-08,6.00E-08,8.00E-08,1.00E-07,1.50E-07,2.00E-07,2.50E-07,3.25E-07,3.50E-07,3.75E-07,4.50E-07,6.25E-07,1.01E-06,1.08E-06,1.13E-06,5.00E-06,6.25E-06,6.50E-06,6.88E-06,7.00E-06,2.05E-05,2.12E-05,2.18E-05,3.60E-05,3.71E-05,6.50E-05,6.75E-05,1.01E-04,1.05E-04,1.16E-04,1.18E-04,1.88E-04,1.92E-04,2.25E-03,3.74E-03,1.70E-02,2.00E-02,5.00E-02,2.00E-01,2.70E-01,3.30E-01,4.70E-01,6.00E-01,7.50E-01,8.61E-01,1.20E+00,1.50E+00,1.85E+00,3.00E+00,4.30E+00,6.43E+00,2.00E+01]
        # evaluates uncertatinty for all available covariance data
        self.reactions_to_evaluate = ['rxn=1', 'rxn=2', 'rxn=4', 'rxn=16', 'rxn=18', 'rxn=102', 'rxn=103', 'rxn=104', 'rxn=105', 'rxn=106', 'rxn=107']
        # names for the reaction used for formating
        self.reaction_names = ['total', 'elastic', 'inelastic','n,2n', 'fission','n,gamma','n,p', 'n,d', 'n,t', 'n,3he', 'n,alpha']
        # percent change in density of the nuclides being perturbed
        self.density_change = 0.01

    
        ### These are things we keep so the code doesn't need to execute them multiple times
        # default neutron library
        self.n_lib = ''
        
        self.path = os.getcwd()
        
        
        if not os.path.exists(self.path+'/working_dir'):
            os.mkdir(self.path+'/working_dir')
        if not os.path.exists(self.path+'working_dir/njoy_run_file.txt'):
            shutil.copyfile('njoy_run_file.txt','working_dir/njoy_run_file.txt')
        if not os.path.exists(self.path+'/sensitivity_plots'):
            os.mkdir(self.path+'/sensitivity_plots')
        if not os.path.exists(self.path+'/covariance_plots'):
            os.mkdir(self.path+'/covariance_plots')
            
            
        if not os.path.exists(self.path+'cov_info.txt'):
            open('cov_info.txt','w')
        else:
            os.remove('cov_info.txt')
            open('cov_info.txt','w')
         
        
    ##### This function creates MCNP strings to add to the MCNP input file
    def make_mcnp_cards(self):
        #opening the mcnp file to read the materials and elements of each material
        template = open(self.mcnp_file_name,'r')
        #creating variables to build later/ setting defaults
        list_of_materials = []
        elements_in_materials = []
        in_materials = False
        j = 0
        #####This loop will read the MCNP input file to create a list of materials and elements of that material
        #This loop steps through the MCNP input line by line
        for line in template:
            #This function creates a list of the elements of the current line
            line_split = line.split()
            #This line looks for an identifier within the MCNP code the tells where the materials list begins
            if "c Materials" in line:
                in_materials = True
                continue
            # this determines if we are in the list of materials
            # if we are, then we begin to create a list of the materials within the system.
            if in_materials == True:
                #This if is true if we are at the first element of a material. This initiallizes the temparary variables:
                #list_of_elements-the element names in MCNP form. EXP:1001.70c is Hydrogen-1 with the 70c cross section library values used
                #element_material- the atom fraction given in the MCNP input material card. The fraction must be in the atom fraction form.    
                if list(line_split[0])[0] == 'm' or list(line_split[0])[0] == 'M':
                    if line_split[0][1]=='t' or line_split[0][1]=='T' in line:
                        continue
                    if 'nlib' in line:
                        if j == 0:
                            self.n_lib = '.'+line.split('=')[1][:-1] 
                            list_of_elements = []
                            element_material = []
                            # the number of the m card given i.e.M21 in MCNP
                            mat_num = line_split[0]
                            j += 1
                        else: 
                            mat_num_list = list(mat_num)
                            temp_num = ''
                            for j in range(len(mat_num_list)):
                                if j>0:
                                    temp_num+= str(mat_num_list[j])
                            temp_list = [temp_num,list_of_elements,element_material]
                            list_of_materials.append(temp_list)
                            list_of_elements = []
                            element_material = []
                            mat_num = line_split[0]
                            j += 1
                        continue
                    #initialization of the first material card within the system
                    if j == 0:
                        list_of_elements = []
                        element_material = []
                        # the number of the m card given i.e.M21 in MCNP
                        mat_num = line_split[0]
                        list_of_elements.append(line_split[1])
                        element_material.append(line_split[2])
                        j += 1
                        #All other materials in the system start with this section
                    else:
                        mat_num_list = list(mat_num)
                        temp_num = ''
                        for j in range(len(mat_num_list)):
                            if j>0:
                                temp_num+= str(mat_num_list[j])
                        temp_list = [temp_num,list_of_elements,element_material]
                        list_of_materials.append(temp_list)
                        list_of_elements = []
                        element_material = []
                        mat_num = line_split[0]
                        list_of_elements.append(line_split[1])
                        element_material.append(line_split[2])
                        j+= 1
                #This line determines if we are at the end of the materials list to end the building of the material library            
                elif'c End Materials' in line:
                    mat_num_list = list(mat_num)
                    temp_num = ''
                    for i in range(len(mat_num_list)):
                        if i>0:
                            temp_num+= str(mat_num_list[i])
                    temp_list = [temp_num,list_of_elements,element_material]
                    list_of_materials.append(temp_list)
                    list_of_elements = []
                    in_materials = False
                #This line skips commented lines
                elif line_split[0] == 'c' or line_split[0] == 'C':
                    continue
                
                #This line adds all additional materials within the system. This condition is most common once in materials
                else:
                    list_of_elements.append(line_split[0])
                    element_material.append(line_split[1])
        template.close()
        #####This section makes a list of the cells where each of the materials are located. It also saves the normal density of the material
        #This section opens the MCNP input
        template = open(self.mcnp_file_name,'r')
        cell_materials = []
        cell_materials_check = []
        #initializes elements of a list for each material within the system
        for materials in list_of_materials:
            cell_materials.append('')
            cell_materials_check.append(False)
        #adds the cell numbers that contain each material
        for line in template:
            # This looks for a call in the MCNP template that tells when the cells section starts
            if line[0] == 'c' or line[0] == ' ':
                continue
            if line == '\n':
                break
            
            line_split=line.split()
            if line_split[0].isdigit():
                for i in range(len(list_of_materials)):
                    if line_split[1] == list_of_materials[i][0]:
                        if cell_materials_check[i] == False:
                            cell_materials[i]+= line_split[0]
                            cell_materials_check[i] = True
                            list_of_materials[i].append(line_split[2])
                        else:
                            cell_materials[i]+= ','+line_split[0]
        
        for i in reversed(range(len(cell_materials_check))):
            if not cell_materials_check[i]:
                del list_of_materials[i]
                del cell_materials[i]

        ##### This converts mass stuctured materials to atomic stuctures
        list_of_materials=self.convert_material_form(list_of_materials)
        
        #elements_present()
        #####This section creates the new added materials with a perturbation in each nuclide. Only used for PERT runs
        list_of_new_materials = []
        elements_used = []
        for x in range(len(list_of_materials)):
            list_temp = ''
            element_temp = []
            for i in range(len(list_of_materials[x][1])):
                for j in range(len(list_of_materials[x][1])):
                    #Creates the first line of the materials that do not have the perturbation on this element
                    if j == 0 and j != i:
                        list_temp+= 'm'+str(list_of_materials[x][0])+str(i+1)+'    '+list_of_materials[x][1][j] + self.n_lib +' '+list_of_materials[x][2][j]+'\n'
                    #Creates the first line of the materials that do have the perturbation on this element
                    elif j == 0 and j == i:
                        list_temp+= 'm'+str(list_of_materials[x][0])+str(i+1)+'    '+list_of_materials[x][1][j] + self.n_lib +' '+str(round(float(list_of_materials[x][2][j])*(1+float(self.density_change)),9))+'\n'
                        element_temp.append(list_of_materials[x][1][j]+ self.n_lib )
                        if not( list_of_materials[x][1][j]+ self.n_lib  in elements_used):
                            elements_used.append(list_of_materials[x][1][j]+ self.n_lib )
                    #Creates the other lines of the materials that do the perturbation on this element
                    elif j!=0 and j == i:
                        list_temp+= '       '+list_of_materials[x][1][j]+ self.n_lib +' '+str(round(float(list_of_materials[x][2][j])*(1+float(self.density_change)),9))+'\n'
                        element_temp.append(list_of_materials[x][1][j]+ self.n_lib )
                        if not( list_of_materials[x][1][j]+ self.n_lib  in elements_used):
                            elements_used.append(list_of_materials[x][1][j]+ self.n_lib )
                    #Creates the other lines of the materials that do not have the perturbation on this element
                    else:
                        list_temp+= '       '+list_of_materials[x][1][j]+ self.n_lib +' '+list_of_materials[x][2][j]+'\n'
            elements_in_materials.append(element_temp)
            list_of_new_materials.append(list_temp)
        
        
        #####This section creates the added MCNP calls that are used to find the sensititvities
        #This section is for creating the pert cards need for source driven sensitivities
        if self.sensitivity_method == 'PERT' :
            list_of_pert = []
            for x in range(len(list_of_materials)):
                temp_list_pert_mat = []
                
                for i in range(len(list_of_materials[x][1])):
                    temp_list_pert = []
                    for j in range(len(self.reactions_to_evaluate)):
                        temp = []
                        for k in range(len(self.energy_bin_structure)-1):
                            pert_string = 'pert'+str(k+1)+':n cell = '+cell_materials[x]+' rho = '+str(round(float(list_of_materials[x][3])+((float(list_of_materials[x][2][i])*(1+float(self.density_change))-float(list_of_materials[x][2][i]))),9))+' mat = '+str(list_of_materials[x][0])+str(i+1)+' '+ self.reactions_to_evaluate[j]+' method = 2'+' erg = '+str(self.energy_bin_structure[k])+' '+str(self.energy_bin_structure[k+1])+'\n'
                            temp.append(pert_string)
                        temp_list_pert.append(temp)
                    temp_list_pert_mat.append(temp_list_pert)    
                list_of_pert.append(temp_list_pert_mat)
        #This section creates cards for ksen keff sensitivities
        elif self.sensitivity_method == 'KSEN':
           list_of_pert = []
           pert_string='ksen1 xs rxn='
           for reaction in self.reactions_to_evaluate:
               pert_string+=reaction[4:]+' '
            
           pert_string+='\n     erg = '
           for e in self.energy_bin_structure:
               pert_string+=str(e)+'\n           '
            
           list_of_pert.append(pert_string[:-9])
           
        else:
            #This is an error report for an improper run name.
            print('No run type selected. Please use KSEN for keff or PERT for sdef.')
            return
        
        # list of cards that needed to be added to the end of the MCNP file
        self.sensitivity_cards = list_of_pert
        # list of materials that need to be added the the MCNP file
        self.materials_to_add = list_of_new_materials
        # list of nuclides present in each material
        self.materials_present = elements_in_materials
        
        self.elements_present=[]
        for material in self.materials_present:
            for element in material:
                if not element in self.elements_present:
                    self.elements_present.append(element)
        return
    
    
    ##### This function converts mass fractions into atomic fractions
    def convert_material_form(self,list_of_materials):
        new_list=[]

        for material in list_of_materials:
            if float(material[2][0])>0:
                new_list.append(material)
                continue
            
            temp_list=[]
            for i in range(len(material[1])):
                materialid=material[1][i]
                
                if '.' in material[1][i]:
                    materialid=material[1][i][:-4]
                    
                temp=[materialid]
                
                A=float(material[1][i][-3:])
                
                awt= -1*float(material[2][i])/A
                
                temp.append(awt)
                temp_list.append(temp)
       
            sum_awt=0
            for l in temp_list:
                sum_awt+=l[1]
        
            for i in range(len(temp_list)):
                material[2][i]=temp_list[i][1]/sum_awt
            
            eff_A=0
            for i in range(len(material[1])):
                materialid=material[1][i]
                
                if '.' in material[1][i]:
                    materialid=material[1][i][:-4]
                    
                A=float(material[1][i][-3:])
                
                eff_A+=A*material[2][i]
                
            material[3]=str(-1*float(material[3])/eff_A*6.02*10**23/10**24)
            
            for i in range(len(material[2])):
                material[2][i]=str(material[2][i]*float(material[3]))
                
            new_list.append(material)
            
        return new_list
    
    
    ##### This function evaluates the NJOY needed to build the covariance library
    def make_NJOY_covariance(self):
        x = ERRORR_tools()
        x.energy_structure = self.energy_bin_structure
        x.energy_structure = [round(a*10**6,6) for a in x.energy_structure]
        nuclides_needed=[]
        for mat in self.materials_present:
            for element in mat:
                if not element in nuclides_needed:
                    nuclides_needed.append(element)
        for element in nuclides_needed:
            element_a=element
            if '.' in element:
                element_a = element[:-4]
            element_a=str(int(element_a))
            x.nuclide_a = element_a[-3:]
            x.nuclide_z = element_a[:-3]
            while len(x.nuclide_z)<3:
                x.nuclide_z = '0'+x.nuclide_z
            i = 0
            with open('element_letters.csv',encoding = 'utf-8-sig') as file:
                for row in file:
                    i+= 1
                    if i == int(x.nuclide_z):
                        x.nuclide_letter = str(row)[:-1]
            
            reaction_list=[]
            for rxn in self.reactions_to_evaluate:
                rxn_num = rxn[4:]
                reaction_list.append(rxn_num)
            x.reaction_list = reaction_list
            print(x.nuclide_letter,x.nuclide_a)
            x.build_library()


    ##### This function takes the cards written by make_mcnp_cards and adds them to the input file
    def write_mcnp_input_file(self):
        
        if self.sensitivity_method == 'KSEN':
            shutil.copyfile(self.mcnp_file_name,'perturbed_mcnp.inp')
            lines_to_append=self.sensitivity_cards[0]
            
            self.job_name='perturbed_mcnp'
            file = open(self.job_name+'.inp','a')
            file.write(lines_to_append)
            file.close()
            self.mcnp_pertubed_file = 'perturbed_mcnp.inp'

        else:
            lines_to_append=''
            for material in self.materials_to_add:
                lines_to_append += material
            
            material_lines = lines_to_append
            old_nucl = 0
            name_string_old=''
            self.job_names = []
            self.out_names = []
            self.list_of_batches = []
            for card_set in self.cards_to_add:
                string=card_set[1].split('_')
                curr_nucl=int(string[1])
                name_string = ''.join(string[:2])
                string=''.join(string)
                                
                if old_nucl != curr_nucl:
                    self.job_name='perturbed_mcnp_'+name_string_old
                    shutil.copyfile(self.mcnp_file_name,'working_dir/'+self.job_name+'.inp')
                    file = open('working_dir/'+self.job_name+'.inp','a')
                    file.write(lines_to_append)
                    file.close()
                    self.job_names.append(self.job_name+'.inp')
                    self.out_names.append(self.job_name+'.out')
                    self.list_of_batches.append('mcnp_batch_script_pert_'+name_string_old)
                    lines_to_append = material_lines
                
                for card in card_set[0]:
                    temp_text = card
                    temp_text = temp_text.replace(':',string+':')
                    lines_to_append+=(temp_text)
                    
                old_nucl = curr_nucl
                name_string_old = name_string
                    
            
            
            self.job_name='perturbed_mcnp_'+name_string_old
            shutil.copyfile(self.mcnp_file_name,'working_dir/'+self.job_name+'.inp')
            file = open('working_dir/'+self.job_name+'.inp','a')
            file.write(lines_to_append)
            file.close()
            self.job_names.append(self.job_name+'.inp')
            self.out_names.append(self.job_name+'.out')
            self.list_of_batches.append('mcnp_batch_script_pert_'+name_string_old)
            self.mcnp_pertubed_file = 'perturbed_mcnp.inp'


    ##### This function writes the batch script for a given string in pert cards
    def write_batch_script(self, string):
        
        with open('mcnp_batch_script_pert','r') as file:
            filedata = file.read()
            filedata = filedata.replace('%%num%%',string)
             
        with open('working_dir/mcnp_batch_script_pert_'+string,'w') as batch:
            batch.write(filedata)
            batch.close()


    ##### This function creates submission scripts and submits MCNP jobs to the UTK NE cluster.
    def submit_to_cluster(self):
        ##This line looks to see if the output file already exists and removes the file if it does
        # This is because MCNP will name the job with the next alphabet letter ie job_ouu.txt instead of job_out.txt
        if self.sensitivity_method == 'KSEN':
            if os.path.exists(self.output_file_name):
                os.system('rm '+self.output_file_name)
                print('file removed out')
                #waiting to verify the cluster knows the file is no longer present
                time.sleep(5)
                
            # This line submits the job via qsub command
            # the part of the string '>> cluster_submissions.txt' will create an output file with the text printed to the command window
            os.system('qsub mcnp_batch_script >> cluster_submissions.txt')
            time.sleep(5)
            template = open('cluster_submissions.txt','r')
            for line in template:
                 if 'socket_connect_unix failed' in line:
                     #if this line is found the job is resubmitted
                     os.remove('cluster_submissions.txt')
                     time.sleep(5)
                     os.system('qsub mcnp_batch_script >> cluster_submissions.txt')
                 else:
                    #if not then the code waits for the output file
                    print('Job successfully submitted')
            template.close()
            os.remove('cluster_submissions.txt')
            os.system('cd ..')
            
        elif self.sensitivity_method == 'PERT':
            for file in self.out_names:
                if os.path.exists('working_dir/' + file):
                    os.system('rm '+'working_dir/' + file)
                    print('file removed out')
                    #waiting to verify the cluster knows the file is no longer present
                    time.sleep(5)
                
            # This line submits the job via qsub command
            # the part of the string '>> cluster_submissions.txt' will create an output file with the text printed to the command window
            for curr_batch in self.list_of_batches:
                os.system('qsub '+'working_dir/'+curr_batch+' >> cluster_submissions.txt')
                time.sleep(5)
                template = open('cluster_submissions.txt','r')
                for line in template:
                     if 'socket_connect_unix failed' in line:
                         #if this line is found the job is resubmitted
                         os.remove('cluster_submissions.txt')
                         time.sleep(5)
                         os.system('qsub mcnp_batch_script >> cluster_submissions.txt')
                     else:
                        #if not then the code waits for the output file
                        print(curr_batch+' successfully submitted')
                template.close()
                os.remove('cluster_submissions.txt')
                os.system('cd ..')


    ##### This function verifies that a given job has completed on the cluster and returns a boolean
    def check_cluster_job(self):
        print('Checking',self.job_name)
        #this loop looks for a file created at the end of a MCNP run
        if not os.path.exists(self.job_name+'_done.dat'):
            job_finish = False   
        else:
            #This shows the job is finished
            print('file found')
            time.sleep(1)
            #removes job_done.dat and cluster_submissions files
            # this will remove the _done.dat file if that is desired. Recommend that you comment out when debugging and uncomment once
            #finished with your coding changes
            # alls you to run the read mechinism on the same output files without having to remake the _done.dat files
            #os.remove(job_name+'_done.dat')
            print('file removed')
            #deletes temp MCNP file
            if os.path.exists('srctp'+self.job_name):
                os.remove('srctp_'+self.job_name)
            job_finish = True
        return job_finish


    ##### This function pulls the data for a covariance matrix from a csv file and creates a matrix that can be used by python
    def pull_covariance_matrix(self):
        x = ERRORR_tools()
        element=self.element
        #removes the .xxc portion of the nuclide name    
        x.element_call = self.element[:-4]
        
        element_a=element
        if '.' in element:
            element_a = element[:-4]
            
        x.nuclide_a = element_a[-3:]
        x.nuclide_z = element_a[:-3]
        while len(x.nuclide_z)<3:
            x.nuclide_z = '0'+x.nuclide_z
        i = 0
        with open('element_letters.csv',encoding = 'utf-8-sig') as file:
            for row in file:
                i+= 1
                if i == int(x.nuclide_z):
                    x.nuclide_letter = str(row)[:-1]
                    break
        x.MT_1 = self.MT_1
        x.MT_2 = self.MT_2
        #converts a list of lists to an array for matrix multiplication when the covariance matrix is used                    
        self.covariance_matrix = x.get_covariance_matrix()


    ##### This function reads the output file of a ksen run of a keff eignevalue run for a sensitivity analysis
    def read_ksen_code(self):
        # if type(self.covariance) ==  int:
        #     return
        #creating a variable the represents the text from the output file
        vector = []
        in_data = False
        in_tally = False
        
        with open(self.output_file_name,'r') as file:
            lines=enumerate(file)
            
            # this loop reads through the output file
            for index,line in lines:        
                # This looks for the string the represents tha start of sensitivity data for the current energy bin
                if self.element + ' ' + self.reaction_call in line:
                    in_tally =  True
                    continue
                
                
                if in_tally == True:
                    # once the sensitivity profile of interest is found, we find where the data starts with this line
                    if'energy range (MeV)         sensitivity   rel. unc.' in line:
                        in_data = True
                        continue
                    
                    if in_data == True:
                        # this verifies we are in the numbers of the data
                        if '.' in line:  
                            start_index=index       
                            in_data = False
                            in_tally = False
                            break
                        
        with open(self.output_file_name,'r') as file:
            counter=0
            lines=enumerate(file)
            
            for index,line in lines:
                if index == start_index:
                    line_split = line.split()
                    # the third element of the line contains the sensitivity, dk/dsig
                    vector.append(float(line_split[2]))
                    start_index+=1
                    counter+=1
                    if counter == (len(self.energy_bin_structure)-1):
                        break
                            
        # creating the sensitivity vector in energy bins, row vector
        self.sensitivity_vector = np.array(vector)
      

    ##### This function reads the output file of a pert run of a source driven run for a sensitivity analysis
    def read_pert_code(self,i,j,k):
        #Initializing temp files
        tally_base = [0,0]
        self.output_file_name='working_dir/perturbed_mcnp_'+str(i)+str(j)+'.out'
        #creates a variable for the output file text
        file_output = open(self.output_file_name,'r')
        in_tally  =  False
        in_data  =  False
        flux = 0
        std = 0
        # reads the text in the file
        ############ This section should be highlighted as the area to read the MCNP tally   
        for line in file_output:
            #looks for the first f14 tally of the MCNP output file
            if '1tally        1'  in line:
                in_tally  =  True
                
            if in_tally:
                # cell 1 is the cell the f14 tally is loctaed in
                #this line appears directly above the data line
                if 'surface  1' in line:
                    in_data = True
                    continue         
                if in_data:
                    in_tally = False
                    line_split = line.split()
                    # This pulls the unperturbed f tally value from the MCNP output file
                    # This value is needed to calculatate the sensitivity
                    flux = float(line_split[0])
                    std = float(line_split[1])
                    break
                    
        file_output.close()
        in_tally = False
        in_data = False
        tally_base[0] = flux
        tally_base[1] = std
        tallies = [tally_base]
        
        # This loop reads the the perturbed f tally values and calculates the sensitivities
        tallies.append([])
        vector = []
        card_string= str(i) + str(j) + str(k)
        
        file_output.close()
        file_output = open(self.output_file_name,'r')
            
        # This sweeps through the energy bins
        for j in range(len(self.energy_bin_structure)-1):
            
            for line in file_output:
                # This is the call that denotes the current perturbation related to the current energy bin
                if 'the following output gives the predicted change in a tally for perturbation' in line and ' '+str(j+1)+card_string+'.' in line:
                    in_tally =  True
                    continue     
                # This denotes the cell we are taking the tally in. This could change based on the type of tally chosen to inspect
                if in_tally == True:
                    if 'surface  1' in line:
                        in_data = True
                        continue
                    #This line does that math required to calculate the sensitivity value
                    if in_data == True:
                        line_split = line.split()
                        # change in f tally
                        c_1 = (float(line_split[0]))/self.density_change
                        # original f tally
                        c_0 = float(tallies[0][0])
                        # creating/ appending sensitivity vector
                        vector.append(c_1/c_0)
                        in_data = False
                        in_tally = False
                        break
        file_output.close()
        # Converts the sensitivity vector from a list to an array in order to matrix multiply
        self.sensitivity_vector = np.array(vector)
        

    ##### This function conducts the sandwich rule using the sensitivities and covariance matrix
    def sandwich_rule(self):
        # Transposing the sensitivity vector
        sens_vect_t = np.transpose(self.sensitivity_vector_1)
        # First multiplication of the sandwich rule SCS
        temp_vect = np.matmul(sens_vect_t,self.covariance_matrix)
        #second multiplication of the sandwich rule
        uncert = np.matmul(temp_vect,self.sensitivity_vector_2)
        #Making the output a list
        self.variance = float(uncert)
     
    ##### This function deterimines the mode of the MCNP file
    def determine_mode(self):
        with open(self.mcnp_file_name,'r') as file:
            text=file.read()
            val=text.find('kcode')
            val1=text.find('KCODE')
            if val1 != -1 or val != -1:
                self.sensitivity_method = 'KSEN'
            else:
                self.sensitivity_method = 'PERT'
        
        
    ##### This function builds the mcnp sensitivity file and runs the file
    def run_mcnp_sensitivity(self):
        
        if os.path.exists('temp'):
            os.remove('temp')
            time.sleep(3)
        
            
        # This determines the run type 
        self.determine_mode()
        
        # This creates the cards needed to calculate the sensitvities
        self.make_mcnp_cards()
        
        # This determines the method for getting covariance data
        if self.covariance_method == "NJOY":
            # This uses NJOY to create a NJOY covariance library
            self.make_NJOY_covariance()
        
        # This removes needless cards if a covariance matrix is not available
        self.cards_to_add=[]
        
        if self.sensitivity_method == 'PERT':
            for i in range(len(self.materials_present)):
                for j in range(len(self.materials_present[i])):
                    for k in range(len(self.reactions_to_evaluate)):
                                            
                        # This gets the covariance data for the nuclide
                        self.element = self.materials_present[i][j]
                        self.reaction = self.reactions_to_evaluate[k]
                        
                        pert_add_on = str(i)+'_'+str(j)+'_'+str(k)
                        self.cards_to_add.append([self.sensitivity_cards[i][j][k],pert_add_on])
        
        elif self.sensitivity_method == 'KSEN':
            self.cards_to_add.append([self.sensitivity_cards[0],''])
            
                        
        # This writes the sensitivity cards into the input file
        self.write_mcnp_input_file()
        
        
        if self.sensitivity_method == 'PERT':
            for run in self.list_of_batches:
                string = run.split('_')[4]
                self.write_batch_script(string)
            
        # This saves the name of the output file if you use the batch script provided
        self.output_file_name = self.job_name+'.out'
        
        # This submits the job
        if self.run_mcnp:
            self.submit_to_cluster()
            
        with open('temp','wb') as file:
            pickle.dump(self,file,pickle.HIGHEST_PROTOCOL)    
        
    
    ##### This function reads the sensitivity vectors and evaluates the uncertainty values
    def evaluate_variance(self):
        
        self.variance_data_table = []
        if self.sensitivity_method == 'KSEN':
            for element in self.elements_present:
                self.element = element
                
                for i in range(len(self.reactions_to_evaluate)):
                    self.reaction = self.reactions_to_evaluate[i]
                    self.reaction_call=self.reaction_names[i]
                    
                    
                    self.read_ksen_code()
                    
                    rxn_num = self.reaction[4:]
                    self.MT_1 = rxn_num
                    self.MT_2 = rxn_num
                    self.pull_covariance_matrix()
                    
                    self.sensitivity_vector_1=self.sensitivity_vector
                    self.sensitivity_vector_2=self.sensitivity_vector
                    self.sandwich_rule()
                    print(self.element +'_' + self.MT_1 + '_' + self.MT_2)
                
                    # This builds a dictionary of uncertainty data 
                    temp_dictionary = {'nuclide':self.element,'reaction':self.reaction_call,'sensitivity':self.sensitivity_vector,'covariance matrix':self.covariance_matrix,'variance':float(self.variance)}
                    self.variance_data_table.append(temp_dictionary)
            
        elif self.sensitivity_method == 'PERT':    
            for card_set in self.cards_to_add:
                pert_add_on = card_set[1]
                
                ints=pert_add_on.split('_')
                i=int(ints[0])
                j=int(ints[1])
                k=int(ints[2])
                
                
                self.read_pert_code(i,j,k)
        
                
                # This gets the covariance matrix and calculates uncertainty
                self.element = self.materials_present[i][j]
                self.reaction = self.reactions_to_evaluate[k]
                self.reaction_call = self.reaction_names[k]
                
                rxn_num = self.reaction[4:]
                self.MT_1 = rxn_num
                self.MT_2 = rxn_num
                self.sensitivity_vector_1=self.sensitivity_vector
                self.sensitivity_vector_2=self.sensitivity_vector
                
                self.pull_covariance_matrix()
                self.sandwich_rule()
                print(self.element +'_' + self.MT_1 + '_' + self.MT_2)
                
                
                added=False
                # This builds a dictionary of uncertainty data
                for i in range(len(self.variance_data_table)):
                    if self.variance_data_table[i]['nuclide'] == self.element:
                        if self.variance_data_table[i]['reaction'] == self.reaction_call:
                            self.variance_data_table[i]['sensitivity'] +=  self.sensitivity_vector  
                            added=True
                            break
                 
                if not added:
                    temp_dictionary = {'nuclide':self.element,'reaction':self.reaction_call,'sensitivity':self.sensitivity_vector.transpose(),'covariance matrix':self.covariance_matrix,'variance':self.variance}
                    self.variance_data_table.append(temp_dictionary)
        
                
            if added:    
                for i in range(len(self.variance_data_table)):
                    sens_1=self.variance_data_table[i]['sensitivity']
                    sens_2=np.transpose(sens_1)
                    temp=np.matmul(sens_2,self.variance_data_table[i]['covariance matrix'])
                    uncert=np.matmul(temp,sens_1)
                    self.variance_data_table[i]['variance']=float(uncert)
            
        #### This section does the non-autocovariance matrices
        for element in self.elements_present:
            for i in range(len(self.reactions_to_evaluate)):
                for j in range(len(self.reactions_to_evaluate)):
                    mt1=self.reactions_to_evaluate[i]
                    mt2=self.reactions_to_evaluate[j]
                    
                    self.MT_1=mt1[4:]
                    self.MT_2=mt2[4:]
                    
                    if float(self.MT_1) >= float(self.MT_2):
                        continue
                    self.element=element
                    self.pull_covariance_matrix()
                    
                    for k in range(len(self.variance_data_table)):
                        if self.variance_data_table[k]['nuclide'] == self.element:
                            if not 'reaction' in self.variance_data_table[k]:
                                continue
                            if self.variance_data_table[k]['reaction'] == self.reaction_names[i]:
                                self.sensitivity_vector_1 = self.variance_data_table[k]['sensitivity']
                                
                            elif self.variance_data_table[k]['reaction'] == self.reaction_names[j]:
                                self.sensitivity_vector_2 = self.variance_data_table[k]['sensitivity']
                                
                                
                    self.sandwich_rule()
                    print(self.element +'_' + self.MT_1 + '_' + self.MT_2)
                
                    temp_dictionary = {'nuclide':self.element,'reaction 1':self.reaction_names[i],'reaction 2':self.reaction_names[j],'covariance matrix':self.covariance_matrix,'variance':self.variance}
                    self.variance_data_table.append(temp_dictionary)
                    
        
    ##### This fucntion prints the nuclides, reactions, and total uncertainty
    def print_total_variance(self):
        tot=0
        if self.sensitivity_method == 'KSEN':
            self.read_ksen_tally()
        else:    
            self.read_pert_tally()
            
        with open('Variance_data.txt','w') as file:
            self.variance_data_table_sorted=sorted(self.variance_data_table, key=lambda x: x['variance'],reverse=True)
            for d in self.variance_data_table_sorted:
                if 'reaction' in d:
                    file.write(f"Variance for the nuclide {d['nuclide']} for reaction {d['reaction']} is: {d['variance']}\n")
                    if not d['reaction'] == 'total':
                        tot+= d['variance']
                else:
                    file.write(f"Covariance for the nuclide {d['nuclide']} for reactions {d['reaction 1']} x {d['reaction 2']} is: {2*d['variance']}\n")
                    if not d['reaction 1'] == 'total':
                        tot+= 2*d['variance']
            file.write(f"Result: {self.tally} \u00B1 {np.sqrt(tot)}")
            
        self.total_variance=tot
       
        
    ##### This fucntion prints the nuclides, reactions, and total uncertainty
    def read_ksen_tally(self):
        with open(self.output_file_name) as file:
            for line in file.readlines():
                if "final result" in line:
                    self.tally=line.split()[2]
                    
                    
    ##### This fucntion prints the nuclides, reactions, and total uncertainty
    def read_pert_tally(self):
        in_tally = False
        in_data = False
        with open(self.output_file_name) as file:
            for line in file.readlines():
                
                if '1tally        1'  in line:
                    in_tally = True
                
                if in_tally:
                    if 'surface  1' in line:
                        in_data = True
                        continue  
                    if in_data:
                        line_split = line.split()
                        self.tally = float(line_split[0])
                        return
                    
                    
    ##### This function is used to automate the process of doing the uncertainty qunatification process for a set of given reaction rates for all of the data that is available
    def run_and_evaluate(self):
        
        # This generates mcnp data needed for the sensititvity vectors
        self.run_mcnp_sensitivity()
        
        # This evaluates the uncertainty
        self.evaluate_variance()
    
        
    ##### This function plots all sensitivity vectors   
    def plot_sensitivity(self,typ):
        x = [e*10**6 for e in self.energy_bin_structure[:-1]]
        
        width = []
        leth = []
        erg = []
        for i in range(len(self.energy_bin_structure)-1):
            width.append((self.energy_bin_structure[i+1]-self.energy_bin_structure[i])*10**6)
            leth.append(1/np.log(self.energy_bin_structure[i+1]/self.energy_bin_structure[i]))
            erg.append(1/(self.energy_bin_structure[i+1]-self.energy_bin_structure[i]))
        
        
        for data in self.variance_data_table:
            if 'reaction 1' in data:
                continue
            if typ == 'dU':    
                y = np.multiply(data['sensitivity'],leth)
            elif typ == 'dE':
                y = np.multiply(data['sensitivity'],erg)
            elif typ == 'E':
                y = data['sensitivity']
            nuc=data['nuclide'].split('.')[0]
            rec=data['reaction']
            
            plt.figure(1)
            plt.bar(x, y, width=width, align='edge')
            plt.xscale('Log')
            plt.xlabel('Energy [eV]')
            if typ == 'dU':    
                plt.ylabel('Sensitivity per unit lethargy')
            elif typ == 'dE':
                plt.ylabel('Sensitivity per MeV')
            elif typ == 'E':
                plt.ylabel('Sensitivity')
            plt.title(nuc+' '+rec)
            plt.tight_layout()
            plt.savefig('sensitivity_plots/'+nuc+'_'+rec+'.png')
            plt.close(1)
        
        
    ##### This function plots all covariance matrices 
    def plot_matrix(self):
        for data in self.variance_data_table:
            if 'reaction' in data:
                nuc=data['nuclide'].split('.')[0]
                rec=data['reaction']
                cov=data['covariance matrix']
                
                
                figure = plt.figure()
                axes = figure.add_subplot(111)
                
                caxes = axes.matshow(cov)
                axes.xaxis.set_ticks_position('bottom')
                figure.colorbar(caxes)
                plt.title(nuc+' '+rec+' covariance')
                plt.gca().invert_yaxis()
                plt.tight_layout()
                plt.savefig('covariance_plots/'+nuc+'_'+rec+'.png')
                plt.close()
            else:
                nuc=data['nuclide'].split('.')[0]
                rec=data['reaction 1']
                rec_2=data['reaction 2']
                cov=data['covariance matrix']
                
                figure = plt.figure()
                axes = figure.add_subplot(111)
                
                caxes = axes.matshow(cov)
                axes.xaxis.set_ticks_position('bottom')
                figure.colorbar(caxes)
                plt.title(nuc+' '+rec+' '+rec_2+' covariance')
                plt.gca().invert_yaxis()
                plt.savefig('covariance_plots/'+nuc+'_'+rec+'_'+rec_2+'.png')
                plt.tight_layout()
                plt.close()
                
                
    ##### This function plots the variances from highest to lowest
    def plot_variance(self,typ):
        list_of_var = []
        list_of_names = []
        num = 0
        
        for data in self.variance_data_table_sorted:
            if data['variance'] < 1e-8:
                break
            
            
            if typ == 'var':
                list_of_var.append((data['variance']))          
            elif typ == 'std':
                list_of_var.append(np.sqrt(data['variance']))
                
                
            nuc=data['nuclide'].split('.')[0]
            if 'reaction' in data:
                rec=data['reaction']
                list_of_names.append(nuc+' '+rec)
            else:
                rec1=data['reaction 1']
                rec2=data['reaction 2']
                list_of_names.append(nuc+' '+rec1+' '+rec2)
            
            num+=1
            if num==5:
                break
        
        rng=np.linspace(1,num,num) 
        
        plt.figure(1)
        plt.bar(rng,list_of_var, width=0.5)
        plt.xticks(rng,list_of_names,rotation=90)
        if typ == 'var':
            plt.ylabel('Variance')        
        elif typ == 'std':
            plt.ylabel('Uncertainty')
        plt.tight_layout()
        plt.savefig('variance.png')
        plt.close()
    
if __name__  == '__main__':
    main = TOFFEE_class()
    time_1 = time.perf_counter()

    main.run_mcnp=False
    main.covariance_method=''
    main.sensitivity_method='PERT'
    main.mcnp_file_name = 'jezebel.inp'
    main.run_mcnp_sensitivity()

    # with open('temp','rb') as file:
    #     main = pickle.load(file) 
    # main.evaluate_variance()
    
    # with open('dict_jezebel','wb') as file:
    #     pickle.dump(main,file,pickle.HIGHEST_PROTOCOL)    
      
    # with open('dict_jezebel','rb') as file:
    #     main = pickle.load(file) 
    
    # main.print_total_variance()
    # main.plot_sensitivity()
    # main.plot_matrix()
    # main.plot_variance()

    
    print('done in ' + str(round(time.perf_counter()-time_1,3)) + ' seconds')
