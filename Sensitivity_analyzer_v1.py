import os
import time
import copy


output_file=open("data_output_k.csv", 'w')
output_file.write("Isotopes,U234,,,U235,,,U236,,,U238,,,H1,,,O16\nReactions,keff,std,,keff,std,,keff,std,,keff,std,,keff,std,,keff,std\n")
output_file.close()

list_of_rxn=['all','rxn 1', 'rxn 2', 'rxn -2', 'rxn 4', 'rxn 18', 'rxn 19','rxn 20','rxn 21','rxn 38','rxn 452','rxn 456','rxn 455','rxn 16','rxn 17','rxn 22','rxn 23','rxn 24','rxn 28','rxn 29','rxn 30','rxn 32','rxn 33','rxn 34','rxn 35','rxn 36','rxn 37','rxn 41','rxn 42','rxn 102','rxn 103','rxn 104','rxn 105','rxn 106','rxn 107','rxn 91' ]
reaction_names=['all','Total', 'Elastic', 'Capture', 'Total inelastic', 'Total Fission', 'First Chance Fission','Second Chance Fission','Third Chance Fission','Fourth Chance Fission','Total Fission v','Prompt Fission v','Delayed Fission v','(n-2n)','(n-3n)','(n-na)','(n-n3a)','(n-2na)','(n-np)','(n-n2a)','(n-2n2a)','(n-nd)','(n-nt)','(n-nhe-3)','(n-nd2a)','(n-nt2a)','(n-4n)','(n-2np)','(n-3np)','(n-gamma)','(n-p)','(n-d)','(n-t)','(n-he-3)','(n-a)','inelastic continum' ]
reaction_calls=['','rxn=1', 'rxn=2', 'rxn=-2', 'rxn=4', 'rxn=18', 'rxn=19','rxn=20','rxn=21','rxn=38','rxn=452','rxn=456','rxn=455','rxn=16','rxn=17','rxn=22','rxn=23','rxn=24','rxn=28','rxn=29','rxn=30','rxn=32','rxn=33','rxn=34','rxn=35','rxn=36','rxn=37','rxn=41','rxn=42','rxn=102','rxn=103','rxn=104','rxn=105','rxn=106','rxn=107','rxn=91' ]

template=open('base_k_temp.txt','r')

material_names=[['1001.70c','H1'],['8016.70c','O16'],['92234.70c','U234'],['92235.70c','U235'],['92236.70c','U236'],['92238.70c','U238']]

energy_bins=[1.00E-10,5.00E-10,7.50E-10,0.000000001,1.20E-09,1.50E-09,0.000000002,2.50E-09,0.000000003,0.000000004,0.000000005,7.50E-09,0.00000001,2.53E-08,0.00000003,0.00000004,0.00000005,0.00000006,0.00000007,0.00000008,0.00000009,0.0000001,0.000000125,0.00000015,0.000000175,0.0000002,0.000000225,0.00000025,0.000000275,0.0000003,0.000000325,0.00000035,0.000000375,0.0000004,0.00000045,0.0000005,0.00000055,0.0000006,0.000000625,0.00000065,0.0000007,0.00000075,0.0000008,0.00000085,0.0000009,0.000000925,0.00000095,0.000000975,0.000001,0.00000101,0.00000102,0.00000103,0.00000104,0.00000105,0.00000106,0.00000107,0.00000108,0.00000109,0.0000011,0.00000111,0.00000112,0.00000113,0.00000114,0.00000115,0.000001175,0.0000012,0.000001225,0.00000125,0.0000013,0.00000135,0.0000014,0.00000145,0.0000015,0.00000159,0.00000168,0.00000177,0.00000186,0.00000194,0.000002,0.00000212,0.00000221,0.0000023,0.00000238,0.00000247,0.00000257,0.00000267,0.00000277,0.00000287,0.00000297,0.000003,0.0000031,0.0000032,0.0000035,0.00000373,0.0000041,0.0000047,0.000005,0.0000054,0.000006,0.00000625,0.0000065,0.00000675,0.000006875,0.000007,0.00000715,0.0000081,0.0000091,0.00001,0.0000115,0.0000119,0.0000129,0.0000144,0.000016,0.000017,0.0000185,0.0000194,0.00002,0.0000205,0.0000212,0.00002175,0.0000225,0.000025,0.0000275,0.00003,0.00003125,0.00003175,0.00003325,0.00003375,0.000035,0.0000355,0.000036,0.000037,0.00003713,0.00003727,0.00003763,0.000038,0.0000391,0.0000396,0.000041,0.0000424,0.000044,0.0000452,0.0000483,0.0000506,0.0000534,0.000058,0.000061,0.000063,0.000065,0.0000675,0.000072,0.000076,0.00008,0.0000817,0.00009,0.000097,0.0001012,0.000105,0.000108,0.000113,0.000116,0.0001175,0.000119,0.000122,0.000143,0.00017,0.00018,0.0001877,0.0001885,0.0001915,0.000193,0.000202,0.0002074,0.0002095,0.00022,0.00024,0.000285,0.000305,0.00055,0.00067,0.000683,0.00095,0.00115,0.0015,0.00155,0.0018,0.0022,0.00225,0.0025,0.003,0.00374,0.0039,0.0057,0.00803,0.0095,0.013,0.017,0.02,0.03,0.045,0.05,0.052,0.06,0.073,0.075,0.082,0.085,0.1,0.1283,0.149,0.2,0.27,0.33,0.4,0.42,0.44,0.47,0.492,0.55,0.573,0.6,0.67,0.679,0.75,0.82,0.8611,0.875,0.9,0.92,1.01,1.1,1.2,1.25,1.317,1.356,1.4,1.5,1.85,2.354,2.479,3,4.304,4.8,6.434,8.187,10,12.84,13.84,14.55,15.68,17.33,20]

def make_pert_cards(template_name,density_change,reactions,erg_bins):
    template=open(template_name,'r')
    list_of_materials=[]
    list_of_elements=[]
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
                    mat_num=line_split[0]
                    list_of_elements.append(line_split[1])
                    j+=1
                else:
                    mat_num_list=list(mat_num)
                    temp_num=''
                    for j in range(len(mat_num_list)):
                        if j>0:
                            temp_num+=str(mat_num_list[j])
                    temp_list=[temp_num,list_of_elements]
                    list_of_materials.append(temp_list)
                    list_of_elements=[]
                    mat_num=line_split[0]
                    list_of_elements.append(line_split[1])
                    j+=1
                        
            elif'c End Materials' in line:
                mat_num_list=list(mat_num)
                temp_num=''
                for i in range(len(mat_num_list)):
                    if i>0:
                        temp_num+=str(mat_num_list[i])
                temp_list=[temp_num,list_of_elements]
                list_of_materials.append(temp_list)
                list_of_elements=[]
                in_materials=False
            elif line_split[0]=='c':
                continue
            else:
                list_of_elements.append(line_split[0])
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
    
                        
            
            
    
    list_of_pert=[]
    for x in range(len(list_of_materials)):
        temp_list_pert_mat=[]
        for i in range(len(list_of_materials[x][1])):
            temp_list_pert=[]
            for j in range(len(reaction_calls)):
                
                pert_string='kpert'+str(j+1)+' cell='+cell_materials[x]+' rho='+str(float(list_of_materials[x][2])-density_change)+' iso='+list_of_materials[x][1][i]+' '+ reactions[j]+ ' erg='
                for erg in erg_bins:
                    pert_string+='      '+str(erg)+'\n'
                temp_list_pert.append(pert_string)
                
                pert_string='kpert'+str(j+101)+' cell='+cell_materials[x]+' rho='+str(float(list_of_materials[x][2])+density_change)+' iso='+list_of_materials[x][1][i]+' '+ reactions[j]+ ' erg='
                for erg in erg_bins:
                    pert_string+='      '+str(erg)+'\n'
                temp_list_pert.append(pert_string)
            temp_list_pert_mat.append(temp_list_pert)
        list_of_pert.append(temp_list_pert_mat)
                        
                
    return list_of_pert

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
    os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir + ' && qsub ' + script_file_string+ '"')
    output_string=job_name+'_out.txt'
    i=0
    while not os.path.exists(job_name+'_done.dat'):
        i+=1
        if i<=100:
            print('Job not complete waiting 15 seconds')
            time.sleep(15)
        elif i==100:
            os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_Dir + ' && qsub ' + script_file_string+ '"')
            i=0
    print('file found')
    os.remove(job_name+'_done.dat')
    print('file removed')
    if os.path.exists('srctp'):
        os.remove('srctp')
    return output_string

def read_k_code(string,reactions,energy_b):
    #string is the name of a txt file that has the output of a kcode run
    tally_base=[0,0]
    print(string)
    file_output=open(string,'r')
    in_tally = False
    in_data = False
    keff=0
    std=0
    for line in file_output:
        if '            first half'  in line:
            in_tally = True
            
        if in_tally:
            if '     second half' in line:
                in_data=True
                continue
            if in_data:
                in_tally=False
                line_split=line.split()
                keff=float(line_split[2])
                
    file_output.close()
    file_output=open(string,'r')
    in_tally=False
    in_data=False

    tally_base[0]=keff
    tally_base[1]=std
    tallies=[tally_base]
    
    i=0
    
    while i in range(len(reactions)):
        in_data=False
        in_tally=False
        num_of_data=1
        tallies.append([])
        for line in file_output:
            if 'kpert   '  in line and str(i+1) in line:
                in_tally= True
                continue
            
            
            if in_tally==True:
                if 'group' in line:
                    in_data=True
                    continue
                if in_data==True and num_of_data<len(energy_b):
                    
                    line_split=line.split()
                    tallies[i+1].append([float(line_split[5]),float(line_split[4])])
                    num_of_data+=1
                    
                if in_data==True and num_of_data==len(energy_b):
                    in_data=False
                    num_of_data=1
                    in_tally=False
                    i=i+1
                    if i<len(reactions):
                        tallies.append([])
                        
    while i in range(len(reactions)):
        in_data=False
        in_tally=False
        num_of_data=1
        tallies.append([])
        for line in file_output:
            if 'kpert   '  in line and str(i+101) in line:
                in_tally= True
                continue
            
            
            if in_tally==True:
                if 'group' in line:
                    in_data=True
                    continue
                if in_data==True and num_of_data<len(energy_b):
                    
                    line_split=line.split()
                    tallies[i+1].append([float(line_split[5]),float(line_split[4])])
                    num_of_data+=1
                    
                if in_data==True and num_of_data==len(energy_b):
                    in_data=False
                    num_of_data=1
                    in_tally=False
                    i=i+1
                    if i<len(reactions):
                        tallies.append([])
            
                    
            
                    
                    
    file_output.close()
    return tallies

def submit_kpert(reaction_list,pert_calls,names,template_file,energies,sub_job=True):
    top_string=','
    top_string_two=','
    top_string_three='Reactions,'
    table=[]
    for i in range(len(pert_calls)):
        top_string+='Material'+str(i+1)
        for j in range(len(pert_calls[i])):
            line_split=pert_calls[i][j][0].split()
            word=line_split[3]
            element_string=''
            for x in range(len(word)):
                if x>=4:
                    element_string+=word[x]
            top_string_two+=element_string+',,,,'
            top_string+=',,,,'
            top_string_three+='keff,std,diff,,'
            temp_dict={'Material':str(i+1),'Element':element_string}
            table.append(temp_dict)
        
    output_file=open("data_output_k.csv", 'w')
    output_file.write(top_string+'\n'+top_string_two+'\n'+top_string_three+'\n')
    output_file.close()
    
    
    
    data=[]
    for i in range(len(pert_calls)):
        data.append([])
        for j in range(len(pert_calls[i])):
            job_num='job'+str(i)+'_'+str(j)
            job=open(job_num+'_in.txt','w')
            write_string=''
            for k in range(len(pert_calls[i][j])):
                write_string+=pert_calls[i][j][k]

            template=open(template_file,'r')
            for line in template:
                job.write(line)
            template.close()
            job.write(write_string)
            job.close()
            if sub_job==True:
                out=submit_to_cluster(job_num)
            else:
                out=job_num+'_out.txt'
            
            data[i].append(read_k_code(out,reaction_list,energies))
            
    list_of_data=[]
    for i in range(len(data[0][0])):
        list_of_data.append([])
    
    list_of_data[0]=','
    for i in range(len(data)):
        #each material
        for j in range(len(data[i])):
            #each element for that material
            list_of_data[0]+=str(data[0][0][0][0][0])+','+str(data[0][0][0][0][1])+',,,'
            for k in range(len(reaction_calls)):

                if i==0 and j==0:
                    list_of_data[k+1]=(str(names[k])+','+str(data[i][j][k+1][0])+','+str(data[i][j][k+1][1])+','+str(abs(data[0][0][0][0]-data[i][j][k+1][0]))+',,')
                else:
                    list_of_data[k+1]+=str(data[i][j][k+1][0])+','+str(data[i][j][k+1][1])+','+str(((data[0][0][0][0]-data[i][j][k+1][0])*10**5))+',,'
    
    for i in range(len(list_of_data)):
        if i==0:
            output_file=open("data_output_k.csv", 'a')
            output_file.write('Base'+list_of_data[i]+'\n')
            output_file.close()
        else:
            output_file=open("data_output_k.csv", 'a')
            output_file.write(list_of_data[i]+'\n')
            output_file.close()
            
    number_of_elements=0
    for i in range(len(data)):
        
        for j in range(len(data[i])):
            table[number_of_elements]['Rxn']=[]
            for k in range(len(data[i][j])):
                if k==0:
                    table[number_of_elements]['Rxn'].append({'Name':'Base','Diff':(data[0][0][0][0]-data[i][j][k][0])*10**5})
                else:
                    table[number_of_elements]['Rxn'].append({'Name':names[k-1],'Diff':(data[i][j][k][0]-data[0][0][0][0])*10**5})
                    
            number_of_elements+=1
    
    full_table=copy.deepcopy(table)
    for out in table:
        removal_list=[]
        for react in out['Rxn']:
            if react['Diff']==0:
                removal_list.append(react)
                
        for elements in removal_list:
            out['Rxn'].remove(elements)
    
            
    return table,full_table

def submit_pert(reaction_list,pert_calls,names,template_file,energies,sub_job=True):
    top_string=','
    top_string_two=','
    top_string_three='Reactions,'
    for i in range(len(pert_calls)):
        top_string+='Material'+str(i+1)
        for j in range(len(pert_calls[i])):
            line_split=pert_calls[i][j][0].split()
            word=line_split[3]
            element_string=''
            for x in range(len(word)):
                if x>=4:
                    element_string+=word[x]
            top_string_two+=element_string+',,,,'
            top_string+=',,,,'
            top_string_three+='keff,std,diff,,'
    
    
        
    output_file=open("data_output_sdef.csv", 'w')
    output_file.write(top_string+'\n'+top_string_two+'\n'+top_string_three+'\n')
    output_file.close()
    
    data=[]
    for i in range(len(pert_calls)):
        data.append([])
        for j in range(len(pert_calls[i])):
            job_num='job'+str(i)+'_'+str(j)
            job=open(job_num+'_in.txt','w')
            write_string=''
            for k in range(len(pert_calls[i][j])):
                write_string+=pert_calls[i][j][k]
            template=open(template_file,'r')
            for line in template:
                job.write(line)
            template.close()
            job.write(write_string)
            job.close()
            if sub_job==True:
                out=submit_to_cluster(job_num)
            else:
                out=job_num+'_out.txt'
            
            data[i].append(read_k_code(out,reaction_list,energies))
            
    
   
    return

    
    
calls=make_pert_cards('base_k_temp.txt',0.1,reaction_calls,energy_bins)
table_output,full=submit_kpert(reaction_calls,calls,reaction_names,'base_k_temp.txt',energy_bins,False)

print('done')
       
