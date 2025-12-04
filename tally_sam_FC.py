# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 21:16:51 2025

@author: jeff
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import gzip
import sys
import time
            
def add_covered_positions(line, contig, covered):
    seqlen = len(line[9])
    pos = int(line[3])
    pos_range = list(range(pos, pos + seqlen))
    covered.extend(list(map(lambda x: contig + '.' + str(x), pos_range)))
    return(covered)

try:
    os.mkdir('map_dets')
except FileExistsError:
    pass

n_reads_df = pd.read_csv('sample_read_counts.csv', index_col = 0)
bin_info = pd.read_csv('hiqh_qual_draft_info.csv', index_col = 0)

col_names = bin_info.index
        
data_rpm = pd.DataFrame(columns = col_names)
data_breadth = pd.DataFrame(columns = col_names)


mag = ''
contig = ''
i = 0
mag_i = 0

try:
    f = sys.argv[1]
    f = f.split('/')[-1]
except IndexError:
    f = '220809_MW_F8_DNA_Crust_A_S75_L002_map.q20.sam.gz'

name = f.split('_map.q20.sam.gz')[0]

with gzip.open('mag_map/' + f, 'rb') as sam_in:
    for line in sam_in.readlines():
        i += 1
        mag_i += 1
        line = line.decode('utf-8')
        line = line.split()
        lmag = line[2].split('|')[0]
        lcontig = line[2].split('|')[1]
        
        ## If first line just initiate values and continue to obtaining positions.
        
        if i == 1:
            mag = lmag
            contig = lcontig
            print(i, name, mag, contig, mag_i)
            covered = []
        
        ## If new mag, collect all data and reinitiate variables.
        
        if lmag != mag:
            time_i = time.time()
            depth = pd.Series(covered).value_counts()
            depth = depth.sort_index()
            frac_cov = float(len(depth)/bin_info.loc[mag, 'size'])
            data_breadth.loc[name, mag] = frac_cov
            
            ## plots really slow things down, turned off for now
            
            #plt.scatter(range(0, len(depth)), depth)
            #plt.title(name + '|' + mag + '|' + str(bin_info.loc[mag, 'size']))
            #plt.savefig('map_dets/' + name + '-' + mag + '.pdf')
            
            depth.to_csv('map_dets/' + name + '-' + mag + '.csv')
            
            time_f = time.time()
            time_delta = time_f - time_i
            print(i, time_delta)
            print(i, mag, contig, mag_i, len(covered))
            
            ## diagnostic tool to identify slow mags
            
            if time_delta > 10:
                #[][1]
                pass
            
            ## clear variables
            
            covered = []
            mag_i = 0
            mag = lmag
            contig = lcontig
            covered = add_covered_positions(line, contig, covered)
       
        ## if len(covered) has grown unreasonably stop adding to it and
        ## continue silent loop until new mag reached
        
        elif len(covered) > 10**8:
            continue            
            
        ## If new contig, obtain new contig name.
        
        elif lcontig != contig:
            contig = lcontig
            print(i, mag, contig, mag_i, len(covered))
            covered = add_covered_positions(line, contig, covered)
            
        else:
            covered = add_covered_positions(line, contig, covered)

data_breadth.to_csv(name + '_FC.csv')
        
        
            