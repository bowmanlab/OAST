# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 07:31:16 2025

@author: jeff
"""

data_path = 'salcedo_merged_binset'
#data_path = 'testing'

import pandas as pd
import os
from Bio import SeqIO

## just count ORFs

i = 0

for f in os.listdir(data_path + '/annotations/prodigal_fna'):
    if f.endswith('fna'):
        for record in SeqIO.parse(data_path + '/annotations/prodigal_fna/' + f, 'fasta'):
            i += 1
            
print(i, 'ORFs to index')

j = 0
            
for f in os.listdir(data_path + '/annotations/kofamscan_processed'):
    if f.endswith('csv'):
        name = f.split('_processed')[0]
        df_temp = pd.read_csv(data_path + '/annotations/kofamscan_processed/' + f)
        df_temp['bin'] = name
        df_temp.index = name + '|' + df_temp['gene']
        
        df_temp.drop('Unnamed: 0', axis = 1, inplace = True)
        
        ## flag secondary annotations and eliminate.  This only impacts a few sequences.
        
        df_temp['secondary_annotation'] = 'False'
        df_temp.loc[df_temp.index.duplicated(keep = False), 'secondary_annotation'] = 'True'
        df_temp.drop_duplicates(subset = 'gene', inplace = True, keep = 'first')
        
        ## open fna file and add sequence info
        
        for record in SeqIO.parse(data_path + '/annotations/prodigal_fna/' + name + '_ORFs.fna', 'fasta'):
            j += 1
            
            if name + '|' + record.id not in df_temp.index:
                print('found', j, 'of', i)
                df_temp.loc[name + '|' + record.id, 'seq'] = str(record.seq) 
                df_temp.loc[name + '|' + record.id, 'length'] = len(record.seq)
            else:
                print(name + '|' + record.id, j, 'of', i)
                df_temp.loc[name + '|' + record.id, 'seq'] = str(record.seq)
                df_temp.loc[name + '|' + record.id, 'bin'] = name
                df_temp.loc[name + '|' + record.id, 'length'] = len(record.seq)
        
        try:
            df_main = pd.concat([df_main, df_temp], axis = 0, ignore_index= False)
        except NameError:
            df_main = df_temp
            
df_main.drop_duplicates(subset = 'seq', keep = False, inplace = True) 
            
df_main.to_csv('orca_salcedo_bins_MT.csv')

with open('orca_salcedo_bins_MT.fasta', 'w') as fasta_out:
    for index, row in df_main.iterrows():
        print('>' + index + '\n' + row.seq, file = fasta_out)
        

