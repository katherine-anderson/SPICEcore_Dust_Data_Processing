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
cfa = pd.read_csv('Cleaned_CFA_Phase1_2019-08-22.csv')
# Make sure data are in correct format
del cfa['Unnamed: 0']

#%%

# Fixing big gaps in depth between measurements
# To keep Python from connecting large gaps with a single line

# Make a copy of the depth and flow data. Don't drop the NaN'd rows.
depth_diff = cfa.loc[:, 'Depth (m)':'Flow Rate'].copy()
# Subtract each depth value from the one before
depth_diff = depth_diff.diff(periods = 1)
# Get rows where difference in depth is >= 6 cm - encoder precision?
bad_diff = depth_diff[(depth_diff['Depth (m)'] >= 0.03)]
# Remove incides of rows with NaN'd flow rate & Abakus data from this list
bad_rows = bad_diff.loc[:, 'Flow Rate'].dropna()
# Get the indices of these rows
bad_rows = list(bad_rows.index.values)
# Change the non-boolean values in these rows to NaN
cfa.loc[bad_rows, 'Depth (m)':'AgeBP'] = np.nan
cfa.loc[bad_rows, 'Sum 1.1-12':'CPP']  = np.nan