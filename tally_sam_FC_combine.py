# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 19:53:48 2025

@author: jeff
"""

import pandas as pd
import glob

files = glob.glob('*FC.csv')

combined = pd.DataFrame()

for f in files:
    temp = pd.read_csv(f, index_col = 0)
    combined = pd.concat([combined, temp])

combined.to_csv('WATL_mg_combined_FC_all.csv')
