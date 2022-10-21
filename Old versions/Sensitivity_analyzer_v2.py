import os
import time
import copy
import numpy as np
import csv


output_file=open("data_output_k.csv", 'w')
output_file.write("Isotopes,U234,,,U235,,,U236,,,U238,,,H1,,,O16\nReactions,keff,std,,keff,std,,keff,std,,keff,std,,keff,std,,keff,std\n")
output_file.close()

list_of_rxn=['rxn 1', 'rxn 2', 'rxn 4', 'rxn 18', 'rxn 102' ]
reaction_names=['Total', 'Elastic', 'Total inelastic', 'Total Fission','(n-gamma)']
reaction_calls=['rxn=1', 'rxn=2', 'rxn=4', 'rxn=18', 'rxn=102' ]


material_names=[['1001.70c','H1'],['8016.70c','O16'],['92234.70c','U234'],['92235.70c','U235'],['92236.70c','U236'],['92238.70c','U238']]

#energy_bins=[1.00E-10,5.00E-10,7.50E-10,0.000000001,1.20E-09,1.50E-09,0.000000002,2.50E-09,0.000000003,0.000000004,0.000000005,7.50E-09,0.00000001,2.53E-08,0.00000003,0.00000004,0.00000005,0.00000006,0.00000007,0.00000008,0.00000009,0.0000001,0.000000125,0.00000015,0.000000175,0.0000002,0.000000225,0.00000025,0.000000275,0.0000003,0.000000325,0.00000035,0.000000375,0.0000004,0.00000045,0.0000005,0.00000055,0.0000006,0.000000625,0.00000065,0.0000007,0.00000075,0.0000008,0.00000085,0.0000009,0.000000925,0.00000095,0.000000975,0.000001,0.00000101,0.00000102,0.00000103,0.00000104,0.00000105,0.00000106,0.00000107,0.00000108,0.00000109,0.0000011,0.00000111,0.00000112,0.00000113,0.00000114,0.00000115,0.000001175,0.0000012,0.000001225,0.00000125,0.0000013,0.00000135,0.0000014,0.00000145,0.0000015,0.00000159,0.00000168,0.00000177,0.00000186,0.00000194,0.000002,0.00000212,0.00000221,0.0000023,0.00000238,0.00000247,0.00000257,0.00000267,0.00000277,0.00000287,0.00000297,0.000003,0.0000031,0.0000032,0.0000035,0.00000373,0.0000041,0.0000047,0.000005,0.0000054,0.000006,0.00000625,0.0000065,0.00000675,0.000006875,0.000007,0.00000715,0.0000081,0.0000091,0.00001,0.0000115,0.0000119,0.0000129,0.0000144,0.000016,0.000017,0.0000185,0.0000194,0.00002,0.0000205,0.0000212,0.00002175,0.0000225,0.000025,0.0000275,0.00003,0.00003125,0.00003175,0.00003325,0.00003375,0.000035,0.0000355,0.000036,0.000037,0.00003713,0.00003727,0.00003763,0.000038,0.0000391,0.0000396,0.000041,0.0000424,0.000044,0.0000452,0.0000483,0.0000506,0.0000534,0.000058,0.000061,0.000063,0.000065,0.0000675,0.000072,0.000076,0.00008,0.0000817,0.00009,0.000097,0.0001012,0.000105,0.000108,0.000113,0.000116,0.0001175,0.000119,0.000122,0.000143,0.00017,0.00018,0.0001877,0.0001885,0.0001915,0.000193,0.000202,0.0002074,0.0002095,0.00022,0.00024,0.000285,0.000305,0.00055,0.00067,0.000683,0.00095,0.00115,0.0015,0.00155,0.0018,0.0022,0.00225,0.0025,0.003,0.00374,0.0039,0.0057,0.00803,0.0095,0.013,0.017,0.02,0.03,0.045,0.05,0.052,0.06,0.073,0.075,0.082,0.085,0.1,0.1283,0.149,0.2,0.27,0.33,0.4,0.42,0.44,0.47,0.492,0.55,0.573,0.6,0.67,0.679,0.75,0.82,0.8611,0.875,0.9,0.92,1.01,1.1,1.2,1.25,1.317,1.356,1.4,1.5,1.85,2.354,2.479,3,4.304,4.8,6.434,8.187,10,12.84,13.84,14.55,15.68,17.33,20]
energy_bins_temp=[1.00E-5,4.00E-03,1.00E-02,2.53E-02,4.00E-02,5.00E-02,6.00E-02,8.00E-02,1.00E-01,1.50E-01,2.00E-01,2.50E-01,3.25E-01,3.50E-01,3.75E-01,4.50E-01,6.25E-01,1.01E+00,1.08E+00,1.13E+00,5.00E+00,6.25E+00,6.50E+00,6.88E+00,7.00E+00,2.05E+01,2.12E+01,2.18E+01,3.60E+01,3.71E+01,6.50E+01,6.75E+01,1.01E+02,1.05E+02,1.16E+02,1.18E+02,1.88E+02,1.92E+02,2.25E+03,3.74E+03,1.70E+04,2.00E+04,5.00E+04,2.00E+05,2.70E+05,3.30E+05,4.70E+05,6.00E+05,7.50E+05,8.61E+05,1.20E+06,1.50E+06,1.85E+06,3.00E+06,4.30E+06,6.43E+06,2.00E+07]
energy_bins=[round(a/10**6,12) for a in energy_bins_temp]

def make_pert_cards(template_name,density_change,reactions,erg_bins,run_type):
    template=open(template_name,'r')
    list_of_materials=[]
    list_of_elements=[]
    elements_out=[]
    in_materials=False
    j=0
    for line in template:
        line_split=line.split()
        if "c Materials" in line:
            in_materials=True
            continue
            
        if in_materials==True:
            if list(line_split[0])[0]=='m':
                if j==0:
                    list_of_elements=[]
                    element_material=[]
                    mat_num=line_split[0]
                    list_of_elements.append(line_split[1])
                    element_material.append(line_split[2])
                    j+=1
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
            elif line_split[0]=='c':
                continue
            else:
                list_of_elements.append(line_split[0])
                element_material.append(line_split[1])
    template.close()
    
    template=open(template_name,'r')
    in_cells=False
    cell_materials=[]
    cell_materials_check=[]
    for materials in list_of_materials:
        cell_materials.append('')
        cell_materials_check.append(False)
    
    for line in template:
        if 'c Cells' in line:
            in_cells= True
            continue
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
        if 'c End Cells'in line:
            in_cells=False
    
    list_of_new_materials=[]
    for x in range(len(list_of_materials)):
        list_temp=''
        element_temp=[]
        for i in range(len(list_of_materials[x][1])):
            for j in range(len(list_of_materials[x][1])):
                if j==0 and j!=i:
                    list_temp+='m'+str(list_of_materials[x][0])+str(i+1)+'    '+list_of_materials[x][1][j]+' '+list_of_materials[x][2][j]+'\n'
                elif j==0 and j==i:
                    list_temp+='m'+str(list_of_materials[x][0])+str(i+1)+'    '+list_of_materials[x][1][j]+' '+str(round(float(list_of_materials[x][2][j])*(1+float(density_change)),9))+'\n'
                    element_temp.append(list_of_materials[x][1][j])
                elif j!=0 and j==i:
                    list_temp+='       '+list_of_materials[x][1][j]+' '+str(round(float(list_of_materials[x][2][j])*(1+float(density_change)),9))+'\n'
                    element_temp.append(list_of_materials[x][1][j])
                else:
                    list_temp+='       '+list_of_materials[x][1][j]+' '+list_of_materials[x][2][j]+'\n'
        elements_out.append(element_temp)    
        list_of_new_materials.append(list_temp)
    
    if run_type=='kpert':
        list_of_pert=[]
        for x in range(len(list_of_materials)):
            temp_list_pert_mat=[]
            for i in range(len(list_of_materials[x][1])):
                temp_list_pert=[]
                for j in range(len(reactions)):
                    temp=[]
                    for k in range(len(erg_bins)-1):
                        pert_string='kpert'+str(k+1)+' cell='+cell_materials[x]+' rho='+str(round(float(list_of_materials[x][3])+((float(list_of_materials[x][2][i])*(1+float(density_change))-float(list_of_materials[x][2][i]))),9))+' mat='+str(list_of_materials[x][0])+str(i+1)+' '+ reactions[j]+' erg= '+str(erg_bins[k])+' '+str(erg_bins[k+1])+'\n'
                        temp.append(pert_string)
                    temp_list_pert.append(temp)
                temp_list_pert_mat.append(temp_list_pert)    
            list_of_pert.append(temp_list_pert_mat)
    elif run_type=='pert' :
        list_of_pert=[]
        for x in range(len(list_of_materials)):
            temp_list_pert_mat=[]
            for i in range(len(list_of_materials[x][1])):
                temp_list_pert=[]
                for j in range(len(reactions)):
                    temp=[]
                    for k in range(len(erg_bins)-1):
                         pert_string='pert'+str(k+1)+':n cell='+cell_materials[x]+' rho='+str(round(float(list_of_materials[x][3])+((float(list_of_materials[x][2][i])*(1+float(density_change))-float(list_of_materials[x][2][i]))),9))+' mat='+str(list_of_materials[x][0])+str(i+1)+' '+ reactions[j]+' method=2'+' erg= '+str(erg_bins[k])+' '+str(erg_bins[k+1])+'\n'
                         temp.append(pert_string)
                    temp_list_pert.append(temp)
                temp_list_pert_mat.append(temp_list_pert)    
            list_of_pert.append(temp_list_pert_mat)
    else:
        print('No run type selected. Please use kpert for keff or pert for sdef.')
                        
                
    return list_of_pert,list_of_new_materials,elements_out

def submit_to_cluster(job_name):
    list_of_jobs_submitted = []
    # The folder this is called in must contain a mcnp input titled job_name_in.txt
    cluster_input_string=' #!/bin/bash \n #PBS -V \n #PBS -q corei7 \n #PBS -l nodes=1:ppn=8 \n hostname \n module load MCNP6/2.0 \n RTP="/tmp/runtp--".`date "+%R%N"` \n cd $PBS_O_WORKDIR \n mcnp6 TASKS 8 I=%%%INPUT%%%_in.txt O=%%%INPUT%%%_out.txt runtpe=$RTP \n grep -a "final result" %%%INPUT%%%_out.txt > %%%INPUT%%%_done.dat \n rm $RTP'
    #print("Creating cluster script for  " + job_name + " job: " + job_name)
    script_string = cluster_input_string.replace("%%%INPUT%%%", job_name)
    script_file_string = job_name + "_script.txt"
    script_file = open(script_file_string, 'w')
    script_file.write(script_string)
    script_file.close()
    time.sleep(5)
    
    #print("Submitting MCNP job: " + job_name)
    list_of_jobs_submitted.append(job_name)
    if os.path.exists(job_name+'_out.txt'):
        os.system('rm '+job_name+'_out.txt')
        print('file removed out')
        time.sleep(25)
        
    # This line actually submits the job
    #while os.path.exists(job_name+'_out.txt')
    #os.remove(job_name+'_out.txt')
    #os.system('qsub ' + script_file_string)
    current_Dir = os.getcwd()
    os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir + ' && qsub ' + script_file_string+ '>> cluster_submissions.txt'+'"')
    output_string=job_name+'_out.txt'
    while not os.path.exists(job_name+'_done.dat'):
        template=open('cluster_submissions.txt','r')
        for line in template:
            if 'socket_connect_unix failed' in line:
                os.remove('cluster_submissions.txt')
                time.sleep(5)
                os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir + ' && qsub ' + script_file_string+ '>> cluster_submissions.txt'+'"')
            else:
                print('Job not complete waiting 15 seconds')
                time.sleep(15)
        template.close()
    
    print('file found')
    time.sleep(15)
    os.remove(job_name+'_done.dat')
    os.remove('cluster_submissions.txt')
    print('file removed')
    if os.path.exists('srctp'):
        os.remove('srctp')
    return output_string

def pull_covariance_matrix(element,reaction):
    
    element_call=element[:-4]
    rxn=reaction[4:]
    matrix_file=open('covariance_library/' +str(element_call)+'_'+str(rxn)+'.csv',encoding='utf-8-sig')
    cov_matrix_temp=[]
    
    file_mat=csv.reader(matrix_file)
    
    for row in file_mat:
        temp=[]
        for i in row:
            temp.append(float(i))
            
        cov_matrix_temp.append(temp)
                        
    cov_matrix=np.array(cov_matrix_temp)
    
    return cov_matrix

def read_k_code(string,reactions,energy_b,element,reaction,perturb):
    #string is the name of a txt file that has the output of a kcode run
    print(string)
    file_output=open(string,'r')
    in_tally = False
    in_data = False

    i=0
    vector=[]
    in_data=False
    in_tally=False
    for line in file_output:
        if 'kpert          delta-rho' in line:
            in_tally= True
            continue
            
            
        if in_tally==True:
            if'1' in line:
                in_data=True

            if in_data==True:
                line_split=line.split()
                vector.append([float(line_split[1])/float(perturb)])
                i+=1            
                if i<len(energy_b)-1:
                    continue
                else:
                    in_data=False
                    in_tally=False
                        
    file_output.close()
      
    sens_vector=np.array(vector)
    sens_vect_t=np.transpose(sens_vector)
    cov_mat=pull_covariance_matrix(element,reaction)
    temp_vect=np.matmul(sens_vect_t,cov_mat)
    uncert=np.matmul(temp_vect,sens_vector)
    file_output.close()
    uncertainty=([uncert])            
            
    return uncertainty

def read_sdef_code(string,element,reaction,energy_b,perturb):
    #string is the name of a txt file that has the output of a kcode run
    tally_base=[0,0]
    print(string)
    file_output=open(string,'r')
    in_tally = False
    in_data = False
    flux=0
    std=0
    for line in file_output:
        if '1tally       14'  in line:
            in_tally = True
            
        if in_tally:
            if 'cell  1' in line:
                in_data=True
                continue
            if in_data:
                in_tally=False
                line_split=line.split()
                flux=float(line_split[0])
                std=float(line_split[1])
                break
                
    file_output.close()
    in_tally=False
    in_data=False

    tally_base[0]=flux
    tally_base[1]=std
    tallies=[tally_base]
    
    file_output=open(string,'r')
    tallies.append([])
    vector=[]
    for j in range(len(energy_b)-1):
        file_output.close()
        file_output=open(string,'r')

        for line in file_output:
            if 'the following output gives the predicted change in a tally for perturbation' in line and ' '+str(j+1)+'.' in line:
                in_tally= True
                continue     
            
            if in_tally==True:
                if'cell  1' in line:
                    in_data=True
                    continue

                if in_data==True:
                    line_split=line.split()
                    c_1=(float(line_split[0]))
                    c_0=float(tallies[0][0])
                        
                    vector.append([c_1/c_0/perturb])
                    in_data=False
                    in_tally=False
        
    sens_vector=np.array(vector)
    sens_vect_t=np.transpose(sens_vector)
    cov_mat=pull_covariance_matrix(element,reaction)
    temp_vect=np.matmul(sens_vect_t,cov_mat)
    uncert=np.matmul(temp_vect,sens_vector)
    file_output.close()
    uncertainty=([uncert])
    return uncertainty

def submit_kpert(reaction_list,pert_calls,new_materials,names,element_list,template_file,energies,density_change,sub_job=True):
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
                vect.append(read_k_code(out,reaction_list,energies,element_cov,react,density_change))
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

    
    
[calls,material_update,element_list]=make_pert_cards('base_k_temp.txt',0.01,reaction_calls,energy_bins,'kpert')
table_output,full=submit_kpert(reaction_calls,calls,material_update,reaction_names,element_list,'base_k_temp.txt',energy_bins,0.01,False)

print('done')
        



