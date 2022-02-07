# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 10:29:43 2021

@author: jarod
"""

#The purpose of this script is to combine information from a fasta file and 3 types of excel sheets into a format specified by the Department of Health

from Bio import SeqIO
import os
import pandas as pd

date = '12.29.2021'



os.chdir("C:/Users/jarod/Desktop/Basespace Analyses/" + date)

#Loading Manifests


file_list=os.listdir("C:/Users/jarod/Desktop/Basespace Analyses/" + date)
manifest = 'manifest'
GISAID = 'GISAID_Submission.csv'
fasta = 'fasta'
lims = 'lims'
GSS = 'GSS'

full_mani = pd.DataFrame([])
full_GISAID = pd.DataFrame([])

for x in range(len(file_list)):
    
    if manifest in file_list[x]:
        
        manifest_csv = pd.read_csv(file_list[x])
        mani_frames = [full_mani,manifest_csv]
        full_mani = pd.concat(mani_frames)
        
    elif GISAID in file_list[x]:
        
        GISAID_csv = pd.read_csv(file_list[x])
        GISAID_frames = [full_GISAID,GISAID_csv]
        full_GISAID = pd.concat(GISAID_frames)
        
    elif lims in file_list[x]:
        
        manifest_csv = pd.read_csv(file_list[x])
        mani_frames = [full_mani,manifest_csv]
        full_mani = pd.concat(mani_frames)
        
    elif GSS in file_list[x]:
        
        GSS_csv = pd.read_csv(file_list[x])
        
        
#Cleaning compiled GSS object, compiled manifest object and creating DOH request object

        
#GSS_csv = GSS_csv[GSS_csv['Dup_Seq'] != 'Remove']
doh_request = pd.read_csv('DOH_new_template_DOH_Submission_' + date + '.csv')
doh_request = pd.DataFrame(index = range(len(GSS_csv)), columns = [
    'LAB_ACCESSION_ID','GISAID_ID','SPECIMEN_COLLECTION_DATE',
    'SUBMITTING_LAB','SEQUENCE_REASON','SEQUENCING_STATUS',
    'PANGO_LINEAGE','FIRST_NAME','LAST_NAME',
    'MIDDLE_NAME','DOB','ALTERNATIVE_ID'])

doh_request.ALTERNATIVE_ID = GSS_csv.AltCoV
full_mani = full_mani.reset_index()
GSS_csv = GSS_csv.reset_index()

#Filling out the spreadsheet with information from various spreadsheets

for g in range(len(full_mani)):
    
    for h in range(len(doh_request.ALTERNATIVE_ID)):
        
        if str(full_mani.investigator_sample_id[g]) == str(doh_request.ALTERNATIVE_ID[h]):
            print(g) 
            doh_request.LAB_ACCESSION_ID[h] = full_mani.specimen_id[g]
            doh_request.GISAID_ID[h] = 'hCoV-19/USA/WA-Altius-' + str(full_mani.investigator_sample_id[g]) + '/2021'
            doh_request.SPECIMEN_COLLECTION_DATE[h] = full_mani.collection_date[g]
            doh_request.SUBMITTING_LAB[h] = 'Altius'
            

for g in range(len(GSS_csv)):
    
    for h in range(len(doh_request.ALTERNATIVE_ID)):
        
        if str(GSS_csv.AltCoV[g]) == str(doh_request.ALTERNATIVE_ID[h]):
            print(g)
            
            doh_request.SEQUENCING_STATUS[h] = GSS_csv['QC Status'][g]
            doh_request.PANGO_LINEAGE[h] = GSS_csv.Lineage[g]
            doh_request.SEQUENCE_REASON[h] = 'sentinel surveillance'
            
for h in range(len(doh_request)):
    
    if doh_request.SEQUENCING_STATUS[h] == 'passed_qc':
        
        doh_request.SEQUENCING_STATUS[h] = 'Complete'


#Modifying FASTA files to have IDs corresponding to the doh_request spreadsheet

records = SeqIO.parse(date + '_GISAID_Submission.fasta', 'fasta')
n = -1
record_set = pd.DataFrame(index = range(len(GSS_csv)), columns = ['sample'])

for record in records:
    n = n + 1
    temp_record = str(record.id)
    record_set.iloc[n] = str(temp_record)
    
record_set = record_set.dropna()

for t in range(len(doh_request)):
    temp_id = str(doh_request.ALTERNATIVE_ID[t])
    
    for record in range(len(record_set)):
        
        if temp_id in str(record_set.iloc[record]):
            
            print(temp_id)
            print(str(record_set.iloc[record]))
            
        elif temp_id == 'nan':
            
            doh_request.GISAID_ID[t] = 'nan'
            
    
for x in range(len(doh_request)):
    
    if doh_request.SEQUENCING_STATUS[x] == 'Failed':
        
        doh_request.GISAID_ID[x] = 'nan'



#Writing the final file
            
doh_request = doh_request.dropna(how = 'all')


doh_request.to_csv(date + '_DOH_Submission.csv', index=False)

