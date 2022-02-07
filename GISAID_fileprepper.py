# -*- coding: utf-8 -*-
"""
Created on Wed May 26 23:40:41 2021

@author: jarod
"""


from Bio import SeqIO
import os
import pandas as pd



#date = easygui.enterbox(msg = "Please provide the date samples were processed in the following format: Month.Day.Year.", title = "Don't worry, be happy - Bobby McFerrin")


date = '12.29.2021'

os.chdir("C:/Users/jarod/Desktop/Basespace Analyses/" + date)


#Loading Manifests


file_list=os.listdir("C:/Users/jarod/Desktop/Basespace Analyses/" + date)
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
        
    
        
        
GSS_csv = GSS_csv[GSS_csv['QC Status'] == 'passed_qc']
GISAID_sub = pd.DataFrame(index = range(len(GSS_csv)), columns = ['submitter','fn','covv_virus_name','covv_type','covv_passage','covv_collection_date','covv_location','covv_add_location','covv_host','covv_add_host_info','covv_sampling_strategy','covv_gender','covv_patient_age','covv_patient_status','covv_specimen','covv_outbreak','covv_last_vaccinated','covv_treatment','covv_seq_technology','covv_assembly_method','covv_coverage','covv_orig_lab','covv_orig_lab_addr','covv_provider_sample_id','covv_subm_lab','covv_subm_lab_addr','covv_subm_sample_id','covv_authors'])
GSS_csv = GSS_csv.reset_index()
GSS_csv = GSS_csv.drop(columns = ['index'])
n = 0
mani_test = pd.DataFrame([])
kingcounty = 'King County'

for y in range(len(GSS_csv)):
    
    mani_test = full_mani[full_mani['investigator_sample_id']==GSS_csv.AltCoV[y]]
    mani_test = mani_test.reset_index()
    
    GISAID_sub.submitter = str('jaroddy')
    GISAID_sub.fn = date + '_GISAID_Submission.fasta'
    GISAID_sub.covv_virus_name[y] = 'hCoV-19/USA/WA-' + GSS_csv.AltCoV[y] + '/2021'
    GISAID_sub.covv_type = 'betacoronavirus'
    GISAID_sub.covv_passage = 'Original'
    
    
    GISAID_sub.covv_collection_date[y] = mani_test.collection_date[0]
    
        
    
    GISAID_sub.covv_location = 'North America / USA / Washington'
    GISAID_sub.covv_host = 'Human'
    
    if mani_test.sex[0] == 'U':
        GISAID_sub.covv_gender[y] = 'unknown'
    else:
        GISAID_sub.covv_gender[y] = mani_test.sex[0]
        
    GISAID_sub.covv_patient_age[y] = mani_test.age[0]
    GISAID_sub.covv_patient_status = 'unknown'
    GISAID_sub.covv_seq_technology = 'Illumina'
    
    if kingcounty in mani_test.facility_name[0]:
        
        GISAID_sub.covv_orig_lab[y] = 'State Testing Facility'
        
    else:
        
        GISAID_sub.covv_orig_lab[y] = mani_test.facility_name[0]
        
    if isinstance(mani_test.facility_address[0], str) == True:
        
        GISAID_sub.covv_orig_lab_addr[y] = mani_test.facility_address[0]
        
    else:
        
        print('found')
        GISAID_sub.covv_orig_lab_addr[y] = 'unknown'
        
    GISAID_sub.covv_subm_lab = 'Altius Institute for Biomedical Research'
    GISAID_sub.covv_subm_lab_addr = '2211 Elliott Avenue, Suite 600 Seattle, WA 98121'
    GISAID_sub.covv_authors = 'Daniel Bates, Rebecca Bruders, Michael Buckley, Mark Frerker, Clem Green, Kneshay Harper, Matt Hartman, Audra Johnson, Jessica Kunder, Lauren Mitchell, Jemma Nelson, Alex Nguyen, Tobias Ragoczy, Joshua Richards, Jacob Rodriguez, John Stamatoyannopoulos, Eric Thorland, Julia Wald'


original_file = r"C:/Users/jarod/Desktop/Basespace Analyses/" + date + "/" + date + "_seq.fasta"
corrected_file = r"C:/Users/jarod/Desktop/Basespace Analyses/" + date + "/" + date + "_GISAID_Submission.fasta"

                 
with open(original_file) as original, open(corrected_file, 'w') as corrected:
        records = SeqIO.parse(original_file, 'fasta')
        for record in records:
            for t in range(len(GSS_csv)):
                if record.id == str(GSS_csv['#'][t]):
                    record.id = GISAID_sub.covv_virus_name[t]
                    record.description = ''
            SeqIO.write(record, corrected, 'fasta')

GISAID_sub.to_csv(date + '_GISAID_Submission.csv', index=False)
    

os.chdir("C:/Users/jarod/Desktop/Basespace Analyses/")
    
#import doh_file_prepper