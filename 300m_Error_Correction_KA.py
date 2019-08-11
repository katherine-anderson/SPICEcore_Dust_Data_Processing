# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 13:54:43 2019

@author: katie
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\SPICE\\Data')
cfa = pd.read_csv('SPICE_final_sync_unfiltered_24JULY2019.csv')
cfa = cfa.replace('#NAME?', np.nan)
cfa = cfa.astype('float')

#%%

# Correcting low Abakus values during melt day 7/19/2016
# Get median of the days before and after
# Apply correction to the 7/19 data to get its median values up to the medians from before/after

# Depths for each melt day
# 7/18: tubes 292 - 301; 290 - 302 m   ---- Day before
# 7/19: tubes 302 - 311; 302 - 312 m   ---- This is the bad day
# 7/20: tubes 312 - 321; 312 - 322 m   ---- Day after

# Select data from the day before
cfa_BandA = cfa.loc[cfa['Depth (m)'] >= 290].copy()
cfa_BandA = cfa_BandA.loc[cfa_BandA['Depth (m)'] <  302].copy()

# Select data from the day after
cfa_A = cfa.loc[cfa['Depth (m)'] >= 312].copy()
cfa_A = cfa_A.loc[cfa_A['Depth (m)'] <  322].copy()

# Combine the data from the day before and after
cfa_BandA = cfa_BandA.append(cfa_A)
# Get median values in each Abakus column for combined data
med_BandA = cfa_BandA.median(axis = 0)
med_BandA = med_BandA.drop(labels = ['Depth (m)', 'Flow Rate', 'ECM'])


# Select data from Jul 19, to correct
cfa_Jul19 = cfa.loc[cfa['Depth (m)'] >= 302].copy()
cfa_Jul19 = cfa_Jul19.loc[cfa_Jul19['Depth (m)'] <  312].copy()
# Get median values in each Abakus column for the day to correct
med_Jul19 = cfa_Jul19.median(axis = 0)
med_Jul19 = med_Jul19.drop(labels = ['Depth (m)', 'Flow Rate', 'ECM'])

# Divide the before/after medians by the Jul 19 medians
# Need to multiply Jul 19 values by the products to get them to the correct levels
# Only want positive, not-NaN correction values
correction = med_BandA.div(med_Jul19)
# Drop NaNs, where it divided by zero
correction = correction.dropna()
# Drop 0s-- don't want to multiply everything by 0
correction = correction[(correction > 0)]

# Get the list of the remaining columns to correct
columns = correction.index.values
# Select these columns in the Jul 19 data
cfa_Jul19 = cfa_Jul19.loc[:, columns]
# Multiply each value from Jul 19 by the correction for its column
new_Jul19 = cfa_Jul19.mul(correction, axis = 1)

# Update the corrected columns in the Jul 19 data
cfa_Jul19.loc[:, columns] = new_Jul19
# New medians now equal the before/after median

# Update original CFA data with the new values in the corrected dataframe
cfa.update(cfa_Jul19)
