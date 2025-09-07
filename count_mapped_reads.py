# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:12:25 2025

@author: jeff
"""

import sys
import pandas as pd

try:
    name = sys.argv[1]
    ref = sys.argv[2]
except IndexError:
    name = '20230701_CTD10_BTL2_RNA_C'
    ref = 'combined_bins_MT.jgicounts.csv'
    print('input error')
    
print(name)

columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL']    
#sam_in = pd.read_csv(name + '_map.q20.sam.gz', sep = '\t', usecols = list(range(0,11)), names = columns)
sam_in = pd.read_csv(name + '_combined_map.sam.gz', sep = '\t', usecols = list(range(0,11)), names = columns,  on_bad_lines = 'skip')
ref_in = pd.read_csv(ref, index_col = 0)

## Drop duplicates so that pairs aren't counted twice

sam_in.drop_duplicates(subset = 'QNAME', inplace = True)

counts = sam_in.RNAME.value_counts()
ref_in[name] = counts
ref_in[name].fillna(0, inplace = True)

## TPM calculatde according to https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

rpk = ref_in[name]/(ref_in['length']/1000)
per_million_sf = sum(rpk)/1000000
ref_in[name + '_TPM'] = rpk/per_million_sf

ref_in.to_csv(ref)
