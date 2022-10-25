# This code will read a endf file and convert it to a np.array
import numpy as np
from os.path import exists
import pathlib
import time


##### This function converts an ENDF written number to a float in python
# ENDF does not typically use the 10^ when writing a number so we must convert them into a real number in python
# ENDF number take the form '1.65378-6' which means 1.65378e-6
###inputs
#line - a number written in the format of endf numbers
def str_2_float(line):
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
        num=line
        
    return num
### output
#num - a float converted from endf data

##### This function creates a string to represent the name of an endf file
# ENDF files take the form: 'n-092_U_235.endf' where the first number is Z and the second number is A
###inputs
#mcnp_number - a element written like it would be used in a MCNP material
def get_file_name(mcnp_number):
    atomic_number=int(mcnp_number[:-3])
    i=0
    # this opens a csv contains pairs of element letters with for assosiated Z's
    with open('element_letters.csv',encoding='utf-8-sig') as file:
        for row in file:
            i+=1
            if i==atomic_number:
                letter=str(row)[:-1]
    # for Z's <100 the number has leading zeros
    # A's do not have to account for this because they are written as such in MCNP
    if atomic_number<10:
        atomic_number_string='00'+str(atomic_number)
    elif atomic_number<100:
        atomic_number_string='0'+str(atomic_number)
    else:
        atomic_number_string=str(atomic_number)
        
    #This line puts all the pieces together to create the ENDF format file name
    endf_file_name='n-'+atomic_number_string+'_'+letter+'_'+mcnp_number[-3:]+'.endf'
    
    
    return endf_file_name
### output
#endf_file_name -  The name of a given endf file written as a string

##### This function creates a covariance matrix for ENDF data where the matrix is written directly
###inputs
#mcnp_number - a element written like it would be used in a MCNP material
#reaction_number- an integer representing the MT number of the reaction desired
def get_cov_mat_norm(mcnp_number,reaction_number):
    # This determines the file name for the isotope of interest
    filename=get_file_name(mcnp_number)
    energy_bins=[]
    number_of_sets=0 
    
    # This gets the element written in the form that it can show up in the ENDF data
    # This should include all forms, but different isotopes could deviate from this form and cause a the data to not be found
    number = int(mcnp_number)
    expon=len(str(number))-1
    if expon>3:
        element_string=str(number/10**expon)+'00'+'+'+str(expon)
        element_string_1=str(number/10**expon)+'01'+'+'+str(expon)
    else:
        element_string=str(number)+'.00000'
        element_string_1=str(number/10**expon)+'000'+'+'+str(expon)
    
    # This determines if the ENDF file exists for an isotope and returns empty data if it does not
    line_num=0
    current_col=0
    in_tally_test=False
    in_tally=False
    path_name=str(pathlib.Path(__file__).parent.resolve())
    if exists(path_name+'/endf_neutron_libraries/'+filename)==False:
        output_matrix=0
        energy_bins=0
        return [output_matrix,energy_bins]
    
    
    # This begins the search for the data in the ENDF file for the isotope 
    file=open('endf_neutron_libraries/'+filename)
    for line in file:
        
        #This checks if the line represents the location of the reaction desired for the given isotope
        if int(line[-4:])==reaction_number:
            # This is to verify if the reaction number is at the start of a data set
            if line[:len(element_string)+1][1:]==element_string:
                in_tally_test=True
                continue
            elif line[:len(element_string)+1][1:]==element_string_1:
                in_tally_test=True
                continue
            
        # This makes sure the data is covariance data
        if in_tally_test==True and in_tally==False:
            str_rxn=str(reaction_number)
            if len(str_rxn)==1:
                str_add='  '+str_rxn
            elif len(str_rxn)==2:
                str_add=' '+str_rxn
            else:
                str_add=str_rxn
                
            if  ' 0.000000+0 0.000000+0          0        '+str_add in line:
                in_tally=True
                line_num=2
                continue
            else:
                in_tally_test=False
                continue
        
        # This helps the code keep up with how much data is left to read
        if in_tally==True:
            line_num+=1
            list_line=list(line)
            #line_num=int(''.join(list_line[-6:-1]))
        
        # The third line of the data gives information about the size of the matrix and energy bins 
        if line_num==3:
            
            temp_str=''.join(list_line[55:66])
            number_of_energy_bins=int(str_2_float(temp_str))
            number_of_matrix_points=int(number_of_energy_bins*(number_of_energy_bins-1)/2)
            current_col=int(number_of_energy_bins-1)
            
            
            
            
        # This section converts the data into the desired format. 
        if line_num>3:
            # 6 data points are given per line so we read the data first
            
            #first number in the line
            temp_str=''.join(list_line[0:11])
            first_num=str_2_float(temp_str)
            
            # second number in the line
            temp_str=''.join(list_line[11:22])
            second_num=str_2_float(temp_str)
            
            # third number in the line
            temp_str=''.join(list_line[22:33])
            third_num=str_2_float(temp_str)
            
            # fourth number in the line
            temp_str=''.join(list_line[33:44])
            fourth_num=str_2_float(temp_str)
            
            # fifth number in the line
            temp_str=''.join(list_line[44:55])
            fifth_num=str_2_float(temp_str)
            
            # sixth number in the line
            temp_str=''.join(list_line[55:66])
            sixth_num=str_2_float(temp_str)
            
            # This section of the code verifies where each data point should go within the energy bins and covariance matrix
            # After his section the covariance data will be a single list of size (N)(N+1)/2
            if number_of_energy_bins>0:
                
                # This section is for reading energy bin data. The data is the N+1 points for an N x N matrix
                if (number_of_energy_bins)/6>=1:
                    
                    temp_list=[first_num,second_num,third_num,fourth_num,fifth_num,sixth_num]
                    energy_bins+=temp_list
                    number_of_energy_bins-=6
                    
                    if number_of_energy_bins==0:
                        number_of_sets+=1
                        data=[]
                        
                elif (number_of_energy_bins)==5:
                    
                    temp_list=[first_num,second_num,third_num,fourth_num,fifth_num]
                    energy_bins+=temp_list
                    
                    number_of_sets+=1
                    number_of_matrix_points-= 1
                    data=[sixth_num]
                    number_of_energy_bins=0
                    
                elif (number_of_energy_bins)==4:
                    temp_list=[first_num,second_num,third_num,fourth_num]
                    energy_bins+=temp_list
                    
                    number_of_sets+=1
                    number_of_matrix_points-= 2
                    data=[fifth_num,sixth_num]
                    number_of_energy_bins=0
                    
                elif (number_of_energy_bins)==3:
                    
                    temp_list=[first_num,second_num,third_num]
                    energy_bins+=temp_list

                    number_of_sets+=1
                    number_of_matrix_points-= 3
                    data=[fourth_num,fifth_num,sixth_num]
                    number_of_energy_bins=0
                    
                elif (number_of_energy_bins)==2:
                    
                    temp_list=[first_num,second_num]
                    energy_bins+=temp_list

                    number_of_sets+=1
                    number_of_matrix_points-= 4
                    data=[third_num,fourth_num,fifth_num,sixth_num]
                    number_of_energy_bins=0
                    
                elif (number_of_energy_bins)==1:
                    
                    temp_list=[first_num]
                    energy_bins+=temp_list
                    
                    number_of_sets+=1
                    number_of_matrix_points-= 5
                    data=[second_num,third_num,fourth_num,fifth_num,sixth_num]
                    number_of_energy_bins=0
            
            # This section reads the data one by one until all of the data has been added to the covariance matrix
            elif number_of_matrix_points>0:
                data.append(first_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(second_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(third_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(fourth_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(fifth_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(sixth_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break      
    file.close()
    
    len_vecs=current_col
    # This checks if the data exists
    if len_vecs<1:
        output_matrix=0
        energy_bins=0
    else:
        i=1
        matrix_temp=[]
        # This converts the data to a lower trianglur matrix
        while current_col>0:
            data_temp=[]
            for i in range(current_col):
                data_temp.append(data[i])
            data=data[current_col:]
            matrix_temp.append(data_temp)
            current_col-=1
            
        #This is an additional verification of the data
        if len(energy_bins)>1:
            for i in range(len(energy_bins)):
                if i>0:
                    if energy_bins[i]<=energy_bins[i-1]:
                        output_matrix=0
                        energy_bins=0
                        return [output_matrix,energy_bins]
                    
        # This mirrors the data across the diagonal to get the proper covariance data as an numpy array
        output_matrix=np.zeros([len_vecs,len_vecs])
        for i in range(len(matrix_temp)):
             for j in range(len(matrix_temp[i])):
                 output_matrix[i+j,j]=float(matrix_temp[i][j])
                 output_matrix[j,i+j]=float(matrix_temp[i][j])    
    return [output_matrix,energy_bins]
### output
#output_matrix - The relative covariance matrix for the given element and reaction
#energy_bins - The energy bin structure for the covariance matrix

##### This function creates a covariance matrix for ENDF data where the matrix is written as a sum of other reactions
###inputs
#mcnp_number - a element written like it would be used in a MCNP material
#reaction_number- an integer representing the MT number of the reaction desired
def get_cov_mat_sum(mcnp_number,reaction_number):
    ###See get_cov_mat_norm for comments that are the same across both functions
    filename=get_file_name(mcnp_number)
    energy_bins=[]
    number = int(mcnp_number)
    expon=len(str(number))-1
    if expon>3:
        element_string=str(number/10**expon)+'00'+'+'+str(expon)
        element_string_1=str(number/10**expon)+'01'+'+'+str(expon)
    else:
        element_string=str(number)+'.00000'
        element_string_1=str(number/10**expon)+'000'+'+'+str(expon)
    
    
    line_num=0
    in_tally_test=False
    in_tally=False
    path_name=str(pathlib.Path(__file__).parent.resolve())
    if exists(path_name+'/endf_neutron_libraries/'+filename)==False:
        output_matrix=0
        energy_bins=0
        return [output_matrix,energy_bins]
    
    file=open('endf_neutron_libraries/'+filename)
    data=[]
    for line in file:
        if int(line[-4:])==reaction_number:
            if line[:len(element_string)+1][1:]==element_string:
                in_tally_test=True
                continue
            elif line[:len(element_string)+1][1:]==element_string_1:
                in_tally_test=True
                continue
        
        if in_tally_test==True and in_tally==False:
            str_rxn=str(reaction_number)
            if len(str_rxn)==1:
                str_add='  '+str_rxn
            elif len(str_rxn)==2:
                str_add=' '+str_rxn
            else:
                str_add=str_rxn
                
            if  ' 0.000000+0 0.000000+0          0        '+str_add in line:
                in_tally=True
                line_num=2
                continue
            else:
                in_tally_test=False
                continue
        

        if in_tally==True:
            line_num+=1
            list_line=list(line)

        # All above comments see get_cov_mat_norm
        # Line 4 gives the number of data points in this set of data
        if line_num==4:
            temp_str=''.join(list_line[44:55])
            number_of_matrix_points=int(temp_str)
            
            
            
            
        # this reads the data and adds it to a list
        if line_num>4:
            
            #first number in the line
            temp_str=''.join(list_line[0:11])
            first_num=str_2_float(temp_str)
            
            # second number in the line
            temp_str=''.join(list_line[11:22])
            second_num=str_2_float(temp_str)
            
            # third number in the line
            temp_str=''.join(list_line[22:33])
            third_num=str_2_float(temp_str)
            
            # fourth number in the line
            temp_str=''.join(list_line[33:44])
            fourth_num=str_2_float(temp_str)
            
            # fifth number in the line
            temp_str=''.join(list_line[44:55])
            fifth_num=str_2_float(temp_str)
            
            # sixth number in the line
            temp_str=''.join(list_line[55:66])
            sixth_num=str_2_float(temp_str)
            
            if number_of_matrix_points>0:
                data.append(first_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(second_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(third_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(fourth_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(fifth_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
                data.append(sixth_num)
                number_of_matrix_points-=1
                if number_of_matrix_points==0:
                    break
    file.close()
    
    
    list_of_cov_matrices=[]
    # The data is listed as pairs
    # The first number is the multiplier of the matrix
    # The second number is the MT number
    for i in range(len(data)):
        # The odd numbers represent MT numbers and i-1 is the multiplier
        if i%2==1:
            tmp_matrix=get_cov_mat_norm(mcnp_number,int(data[i]))
            tmp_matrix[0]*=data[i-1]
            list_of_cov_matrices.append(tmp_matrix)
    
    # This section verifies if the covariance matrices have the same structure as on another
    # If the do not, they can not be simple added together and a 0 matrix must be returned
    # In the future, NJOY can help make sure the data is uniform
    default_erg_bins=list_of_cov_matrices[0][1]
    default_matrix=list_of_cov_matrices[0][0]*0
    for x in list_of_cov_matrices:
        if x[1]==default_erg_bins:
            default_matrix+=x[0]
        else:
            return[0,0]
    return [default_matrix,default_erg_bins]
### output
#output_matrix - The relative covariance matrix for the given element and reaction
#energy_bins - The energy bin structure for the covariance matrix

##### This function creates a covariance matrix for ENDF data where the matrix is written as a sum of other reactions
###inputs
#mcnp_number - a element written like it would be used in a MCNP material
#reaction_number- an integer representing the MT number of the reaction desired
def determine_type(mcnp_number,reaction_number):
    filename=get_file_name(mcnp_number)
    out='No Data'
    number = int(mcnp_number)
    expon=len(str(number))-1
    if expon>3:
        element_string=str(number/10**expon)+'00'+'+'+str(expon)
        element_string_1=str(number/10**expon)+'01'+'+'+str(expon)
    else:
        element_string=str(number)+'.00000'
        element_string_1=str(number/10**expon)+'000'+'+'+str(expon)
    
    in_tally_test=False
    in_tally=False
    path_name=str(pathlib.Path(__file__).parent.resolve())
    if exists(path_name+'/endf_neutron_libraries/'+filename)==False:
        return out
    
    file=open('endf_neutron_libraries/'+filename)
    for line in file:
        if int(line[-4:])==reaction_number:
            if line[:len(element_string)+1][1:]==element_string:
                in_tally_test=True
                continue
            elif line[:len(element_string)+1][1:]==element_string_1:
                in_tally_test=True
                continue
        
        if in_tally_test==True and in_tally==False:
            str_rxn=str(reaction_number)
            if len(str_rxn)==1:
                str_add='  '+str_rxn
            elif len(str_rxn)==2:
                str_add=' '+str_rxn
            else:
                str_add=str_rxn
                
            if  ' 0.000000+0 0.000000+0          0        '+str_add in line:
                out= str(float(line[44:55]))
                break
            else:
                in_tally_test=False
                continue
    
    
    return out
### output
#output_matrix - The relative covariance matrix for the given element and reaction
#energy_bins - The energy bin structure for the covariance matrix

def get_cov_mat(mcnp_number,reaction_number):
    cov_type=determine_type(mcnp_number,reaction_number)
    if cov_type=='1.0':
        [output_matrix,energy_bins]=get_cov_mat_sum(mcnp_number,reaction_number)
    elif cov_type=='No Data':
        return [0,0]
    else:
        [output_matrix,energy_bins]=get_cov_mat_norm(mcnp_number,reaction_number)
    
    return [output_matrix,energy_bins]

time_1=time.time()
num='40090'
test=get_cov_mat(num,16)

print(time.time()-time_1)