# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 14:21:12 2019

@author: katie
"""

#%%
# Import modules
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# Load function definition script
os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\Git Projects\\SPICE_Data_Processing')
# Need to read this script for the 'select_cfa' function
exec(open('SPICE_Data_Processing_Functions.py').read())

# Change to SPICE data directory
os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\SPICE\\Data')
#%%
# Load cleaned CFA data, delete extra columns
cfa = pd.read_csv('Cleaned_CFA_Phase2_2019-08-26.csv')
del cfa['Unnamed: 0']
#%%
# Set up interpolation dataframe
interp = pd.DataFrame()

# Need to get lists of incrementing age ranges
start_age = []
end_age = []
# Get youngest age value to start
years = cfa.loc[0, 'AgeBP']

# Increment ages until the max age of the data
while years <= 54302:
    # Add the start age to the list
    start_age.append(years)
    # Increment by 500 years
    years = years + 500
    # Add the end age to the list
    end_age.append(years)

# Add lists to new interp dataframe columns
interp['Starting Age'] = start_age
interp['Ending Age']   = end_age
# Create empty column to indicate the # of years to interpolate over in that interval
interp['Years to Interp'] = np.nan

# Subset the data for different time intervals (arbitrary)
# Set different limits on the # of years to interpolate for each interval
holocene = interp[interp['Starting Age'] <= 12000].copy()
holocene['Years to Interp'] = 0.5

lgm = interp[interp['Starting Age'] > 12000].copy()
lgm = lgm[lgm['Starting Age'] <= 22000]
lgm['Years to Interp'] = 1

glacial = interp[interp['Starting Age'] > 22000].copy()
glacial['Years to Interp'] = 5

# Merge all 3 dataframes, now with the 'years to interp' data
dfs = [holocene, lgm, glacial]
interp = pd.concat(dfs)
#%%

# Create empty list to store the # of rows to interpolate over for each interval
interp_rows = []

# Loop through each time interval in the interpolation dataframe
for start, end, years in zip(interp['Starting Age'], interp['Ending Age'], interp['Years to Interp']):
    # Subset the CFA data for this age interval
    temp = select_cfa(cfa, start, end, 'AgeBP')
    # Only keep the age values
    temp = temp.loc[:, 'AgeBP'].copy()
    # Get the difference in age between each row
    age_diff = temp.diff()
    # Get the mean age difference between rows
    mean_diff = np.nanmean(age_diff)
    # Find how many rows are needed to reach the specified number of years
    number_rows = years / mean_diff
    # Round to nearest integer. This is the number of rows to interpolate over for this interval
    number_rows = np.round(number_rows, decimals = 0)
    # Change to integer
    number_rows = number_rows.astype(np.int64)
    # Append the number of rows
    interp_rows.append(number_rows)

# Add list to new column in the interpolation dataframe
interp['Rows to Interp'] = interp_rows    

#%%
# INTERPOLATION
# Create empty series for final concentration & CPP columns
final_conc = pd.Series()
final_cpp  = pd.Series()

# Loop through each interpolation interval
for start, end, rows in zip(interp['Starting Age'], interp['Ending Age'], interp['Rows to Interp']):
    # Subset the CFA data for this age interval
    temp = select_cfa(cfa, start, end, 'AgeBP')
    # Interpolate concentration up to x number of rows, calculated above
    # Limit_direction = backward, so it interpolates forward through time
    # This is where to play around with the interpolation method
    interp_conc = temp['Sum 1.1-12'].interpolate(method = 'spline', order = 3, limit = rows, limit_direction = 'backward')
    # Interpolate CPP
    interp_cpp = temp['CPP'].interpolate(method = 'spline', order = 3, limit = rows, limit_direction = 'backward')    
    # Append interpolated values to final lists before moving to next interval
    final_conc = final_conc.append(interp_conc)
    final_cpp  = final_cpp.append(interp_cpp)

# Create a new copy of the CFA data, for comparison
new = cfa.copy()
# Update concentration & CPP columns with interpolated values
new['Sum 1.1-12'] = final_conc
new['CPP'] = final_cpp    
#%%