# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 16:08:48 2021

@author: jarod
"""


from Bio import SeqIO
import os
import pandas as pd





os.chdir("C:/Users/jarod/Desktop/Basespace Analyses/GenBank")


#Loading Manifests


file_list=os.listdir("C:/Users/jarod/Desktop/Basespace Analyses/GenBank")
manifest = 'manifest'
GSS = 'GSS'
fasta = 'fasta'
lims = 'lims'

full_mani = pd.DataFrame([])

for x in range(len(file_list)):
    
    if manifest in file_list[x]:
        
        manifest_csv = pd.read_csv(file_list[x])
        mani_frames = [full_mani,manifest_csv]
        full_mani = pd.concat(mani_frames)
        
    elif GSS in file_list[x]:
        
        GSS_csv = pd.read_csv(file_list[x])
        
    elif lims in file_list[x]:
        
        manifest_csv = pd.read_csv(file_list[x])
        mani_frames = [full_mani,manifest_csv]
        full_mani = pd.concat(mani_frames)
        

full_mani = full_mani.reset_index()
        
        
GSS_csv = GSS_csv[GSS_csv['Status'] == 'passed_qc']
GSS_csv = GSS_csv[GSS_csv['Seq_Quality'] == 'Pass']
GSS_csv = GSS_csv[GSS_csv['Dup_Seq'] != 'Remove']

#GISAID_sub = pd.DataFrame(index = range(len(GSS_csv)), columns = ['submitter','fn','covv_virus_name','covv_type','covv_passage','covv_collection_date','covv_location','covv_add_location','covv_host','covv_add_host_info','covv_sampling_strategy','covv_gender','covv_patient_age','covv_patient_status','covv_specimen','covv_outbreak','covv_last_vaccinated','covv_treatment','covv_seq_technology','covv_assembly_method','covv_coverage','covv_orig_lab','covv_orig_lab_addr','covv_provider_sample_id','covv_subm_lab','covv_subm_lab_addr','covv_subm_sample_id','covv_authors'])
GSS_csv = GSS_csv.reset_index()
GSS_csv = GSS_csv.drop(columns = ['index'])
n = 0



new_sequences = []
records = SeqIO.parse("genbank_seq.fasta", 'fasta')
n = 0


# =============================================================================
# for record in records:
#     for t in range(len(GSS_csv)):
#         if record.id == GSS_csv.LN[t]:
#             record.id = GSS_csv.AltCoV[t]
#             record.description = ''
#             new_sequences.append(record)
#             n = n + 1
#             print(n)
#             
#     SeqIO.write(new_sequences, "GenBank_submission.fasta", 'fasta')
# =============================================================================
    

    
tracker = pd.read_csv('source_modifiers.csv')

tracker.country = 'USA'
tracker.host = 'Human'

for x in range(len(tracker)):
    print('x' + str(x))
    for seq in range(len(GSS_csv)):
        
        if tracker.isolate[x] == GSS_csv.AltCoV[seq]:
            
            tracker.collection_date[x] = GSS_csv.Date_Collected[seq]
            
for y in range(len(tracker)):
    print('y' + str(y))
    for source in range(len(full_mani)):
        
        if tracker.isolate[y] == full_mani.investigator_sample_id[source]:
            
            tracker.isolation_source[y] = full_mani.specimen_id[source]
    
    