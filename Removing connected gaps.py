# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 13:54:43 2019

@author: katie
"""

import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt

#%%
# Ask user for directory where data are located
directory = input('Enter path for SPICEcore data: ')
os.chdir(directory)

# Load CFA file
cfa = pd.read_csv('SPICE_final_sync_unfiltered_24JULY2019.csv')
# Make sure data are in correct format
cfa = cfa.replace('#NAME?', np.nan)
cfa = cfa.astype('float')

#%%

# Fixing big gaps in depth between measurements
# To keep Python from connecting large gaps with a single line

# Subtract each row from the one before
diff = cfa.diff(periods = 1, axis = 0)
# Get rows where difference in depth is >= 6 cm - encoder precision?
bad_diff = diff.drop(diff[diff['Depth (m)'] < 0.06].index)
# Only keep the depth column.
bad_diff = bad_diff.loc[:, 'Depth (m)']
# Drop the NaNs in the depth column
bad_diff = bad_diff.dropna()
# Get indices of the rows to NaN
# NaN the first row after a gap in depths >= 6 cm
bad_rows = list(bad_diff.index.values)
# Change values in these rows to NaN
cfa.loc[bad_rows, :] = np.nan