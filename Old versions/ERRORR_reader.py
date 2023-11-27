#This program is a tool to create an input file and read the output of the ERRORR file from NJOY
import shutil
import linecache
import time
import os
import numpy as np

class ERRORR_tools:
    #isotope_z- This refers to the z value of the isotope
    #isotope_letter- This refers to the character abbrivation of the isotope
    #isotope_a- This refers to the a value of the isotope
    #MT_1- This refers to the MT value for the first reaction
    #MT_2- This refers to the MT value for the second reaction if you are not looking at self covariance
    #File_name- This refers to the name of the endf file used must be of the type "tapexx" ie "tape30"
    #energy_structure- This refers to the energy bin scheme used for the njoy files
    def __init__(self):
        self.isotope_z="092"
        self.isotope_letter='U'
        self.isotope_a="238"
        self.MT_1=1
        self.MT_2=1
        self.file_name='tape30'
        self.endf_folder='endf_neutron_libraries/'
        self.energy_structure=[1.00E-5,4.00E-03,1.00E-02,2.53E-02,4.00E-02,5.00E-02,6.00E-02,8.00E-02,1.00E-01,1.50E-01,2.00E-01,2.50E-01,3.25E-01,3.50E-01,3.75E-01,4.50E-01,6.25E-01,1.01E+00,1.08E+00,1.13E+00,5.00E+00,6.25E+00,6.50E+00,6.88E+00,7.00E+00,2.05E+01,2.12E+01,2.18E+01,3.60E+01,3.71E+01,6.50E+01,6.75E+01,1.01E+02,1.05E+02,1.16E+02,1.18E+02,1.88E+02,1.92E+02,2.25E+03,3.74E+03,1.70E+04,2.00E+04,5.00E+04,2.00E+05,2.70E+05,3.30E+05,4.70E+05,6.00E+05,7.50E+05,8.61E+05,1.20E+06,1.50E+06,1.85E+06,3.00E+06,4.30E+06,6.43E+06,2.00E+07]
        self.covariance_matrix=np.zeros((int(len(self.energy_structure))-1,int(len(self.energy_structure))-1))
        path=os.getcwd()
        if not os.path.exists(path+'/working_dir'):
            os.mkdir(path+'/working_dir')
        if not os.path.exists(path+'/njoy_covariance'):
            os.mkdir(path+'/njoy_covariance')
      
    #This function is used to create an input file for NJOY on the NEcluster
    #inputs
    #isotope_z- This refers to the z value of the isotope
    #isotope_letter- This refers to the character abbrivation of the isotope
    #isotope_a- This refers to the a value of the isotope
    #MT_1- This refers to the MT value for the first reaction
    def create_input_file(self): 
        ## This section copies the endf file from the library assumed to be located in the current directory
        
        endf_file_name=self.endf_folder+'n-'+self.isotope_z+'_'+self.isotope_letter+'_'+self.isotope_a+'.endf'
        shutil.copyfile(endf_file_name,'working_dir/'+ self.file_name)
        ## This section converts the energy bins into the format for NJOY
        energy_bin_length=len(self.energy_structure)-1
        energy_bin_data=''
        temp_string=''
        for erg_bin in self.energy_structure:
            current_string=str(erg_bin)+' '
            if len(temp_string+current_string)<72:
                temp_string+=current_string
            else:
                energy_bin_data+=temp_string+'\n'
                temp_string=current_string
        energy_bin_data+=temp_string
        ## This section gets the endf data name from the tape created 
        linecache.clearcache()
        endf_data_name=linecache.getline('working_dir/'+self.file_name, 2)[66:70]

        
        ## This section writes the njoy run file
        with open('njoy_template.txt','r') as file:
            data=file.read()
            
            data=data.replace('%%isotope_endf_value%%',endf_data_name)
            data=data.replace('%%energy_bins_length%%',str(energy_bin_length))
            data=data.replace('%%energy_bins%%',energy_bin_data)
            data=data.replace('%%mt_number%%',str(self.MT_1))
            
                
        with open('working_dir/njoy_run_file.txt','w') as file:
            file.write(data)
            
    #This function is used to create an input file for NJOY on the NEcluster
    #inputs
    #isotope_z- This refers to the z value of the isotope
    #isotope_letter- This refers to the character abbrivation of the isotope
    #isotope_a- This refers to the a value of the isotope
    #MT_1- This refers to the MT value for the first reaction
    def create_input_file_specific(self): 
        ## This section copies the endf file from the library assumed to be located in the current directory
        
        endf_file_name=self.endf_folder+'n-'+self.isotope_z+'_'+self.isotope_letter+'_'+self.isotope_a+'.endf'
        shutil.copyfile(endf_file_name,'working_dir/'+ self.file_name)

        ## This section converts the energy bins into the format for NJOY
        energy_bin_length=len(self.energy_structure)-1
        energy_bin_data=''
        temp_string=''
        for erg_bin in self.energy_structure:
            current_string=str(erg_bin)+' '
            if len(temp_string+current_string)<72:
                temp_string+=current_string
            else:
                energy_bin_data+=temp_string+'\n'
                temp_string=current_string
        energy_bin_data+=temp_string
        ## This section gets the endf data name from the tape created 
        endf_data_name=linecache.getline('working_dir/'+self.file_name, 2)[66:70]

        
        ## This section writes the njoy run file
        with open('njoy_template_specific.txt','r') as file:
            data=file.read()
            
            data=data.replace('%%isotope_endf_value%%',endf_data_name)
            data=data.replace('%%energy_bins_length%%',str(energy_bin_length))
            data=data.replace('%%energy_bins%%',energy_bin_data)
            data=data.replace('%%mt_number%%',str(self.MT_1))
            
               
        with open('working_dir/njoy_run_file.txt','w') as file:
            file.write(data)
        
    def submit_njoy_run(self):
        # This line submits the job via qsub command
        # the part of the string '>> cluster_submissions.txt' will create an output file with the text printed to the command window
        current_Dir = os.getcwd()
        os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir +'/working_dir' + ' && qsub njoy_run_file.txt >> cluster_submissions.txt'+'"')
        time.sleep(5)
        template=open('working_dir/cluster_submissions.txt','r')
        for line in template:
            if 'socket_connect_unix failed' in line:
                #if this line is found the job is resubmitted
                os.remove('working_dir/cluster_submissions.txt')
                time.sleep(5)
                os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir +'/working_dir' + ' && qsub njoy_run_file.txt >> cluster_submissions.txt')
            else:
                #if not then the code waits for the output file
                outfile='njoy_run_file.txt.e'+(line.split('.')[0])
                print('Job successfully submitted')
        template.close()
        os.remove('working_dir/cluster_submissions.txt')
        job_finish=True
        while job_finish:
            if os.path.exists('working_dir/'+outfile):
                #This shows the job is finished
                os.remove('working_dir/'+outfile)
                print('file removed')
                job_finish=False
        
        
    def read_endf_style_number(self,line):
        #The loop determines if a value is positive or negative and the magnitude of the exponential
        if '-' in line:
            num_list=line.split('-')
            if len(num_list)==3:
                #This is the case for a negative number with a negative exponent
                num=-1*float(num_list[1])*10**(-1*float(num_list[2]))
            else:
                if '+' in num_list[1]:
                    #This is the case for a negative number with a positive exponent
                    num_list=line.split('+')
                    num=float(num_list[0])*10**(float(num_list[1]))
                else:
                    #This is the case for a positive number with a negative exponent
                    num=float(num_list[0])*10**(-1*float(num_list[1]))
        elif '+' in line:
            #This is the case for a positive number with a positive exponent
            num_list=line.split('+')
            num=float(num_list[0])*10**(float(num_list[1]))
        else:
            num=float(line.strip() or 0)
            
        return num
        
    def read_output_file(self):
        with open('working_dir/tape33','r') as file:
            matrix_size=len(self.energy_structure)-1
            self.covariance_matrix=np.zeros((matrix_size,matrix_size))
            in_covariance=False
            in_data=False
            get_row_number=True
            for line in file:
                if not in_covariance:
                    if line[70:72]=='33':
                        in_covariance=True
                    continue
                if in_covariance and not in_data:
                    if int(line[72:75])==int(self.MT_1) and int(line[35:46])==int(self.MT_2):
                        in_data=True
                        
                    continue
                if in_data:
                    if get_row_number:
                        get_row_number=False
                        i_index=int(line[55:66])
                        j_index=int(line[35:46])
                        j_values=int(line[45:56])
                        continue
                    matrix_elements=[line[:11],line[11:22],line[22:33],line[33:44],line[44:55],line[55:66]]
                    matrix_elements=[self.read_endf_style_number(a) for a in matrix_elements]
                    while j_values>0 and matrix_elements!=[]:
                        self.covariance_matrix[i_index-1,j_index-1]=round(matrix_elements[0],9)
                        j_values-=1
                        j_index+=1
                        matrix_elements.pop(0)
                    if j_values==0:
                        get_row_number=True
                        if i_index==matrix_size:
                            break
                        
    def get_covariance_matrix(self):
        covariance_file_name='njoy_covariance/'+str(self.isotope_letter)+'_'+str(self.isotope_a)+'_'+str(self.MT_1)+'_'+str(self.MT_2)+'.csv'
        if os.path.exists(covariance_file_name):
            self.covariance_matrix=np.loadtxt(covariance_file_name, delimiter=',')
        else:
            if self.MT_1==1 or self.MT_1==102:
                self.create_input_file()
            else:
                #self.create_input_file_specific()
                self.create_input_file()
            self.submit_njoy_run()
            self.read_output_file()
            file=open(covariance_file_name,'w')
            file.close()
            np.savetxt(covariance_file_name, self.covariance_matrix, delimiter=',')
        
if __name__== '__main__':
    x=ERRORR_tools()
    x.isotope_a='016'
    x.isotope_z='008'
    x.isotope_letter='O'
    x.MT_1=1
    x.MT_2=1
    x.get_covariance_matrix()
    print(x.covariance_matrix[0,0])
    
    


