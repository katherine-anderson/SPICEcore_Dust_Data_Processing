# ------------------------------------------------------------------------------------------------------
#                                    FINAL SPICE REMOVE ERRORS SCRIPT
# Removes melting errors in the raw CFA data, interpolates timescale, and adds descriptive columns
#
#    - Loads raw, unfiltered CFA with depth corrections
#    - Tracks the number of measurements NaN'ed in each step
#    1) NaNs bubbles using liquid conductivity data
#    2) NaNs liquid conductivity values < 0.6
#    3) NaNs measurements without positive flow rates
#    4) NaNs measurements with depth duplicates or decreases
#    5) NaNs any negative and infinite Abakus values
#    6) Applies correction to Abakus data from bad melt day (7/19/2016))
#    7) Creates timescale for CFA data (annual in Holocene, tie points for glacial)
#    8) Labels all measurements within specified depth ranges of core breaks
#    9) Labels all measurements within volcanic events and dust events
#   10) Calculates particle concentration and coarse particle percentage (CPP)
#   11) Exports cleaned dataset to CSV
#
# Katie Anderson and Aaron Chesler, 8/13/19
# ------------------------------------------------------------------------------------------------------

#%%
# ------------------------------------------------------------------------------------------------------
#                                           1: FILE PREPARATION
# ------------------------------------------------------------------------------------------------------

# Import needed modules & packages
import pandas as pd
import os 
import numpy as np
from datetime import date

# Ask user for directory where scripts are located
directory = input('Enter path for SPICEcore scripts: ')
os.chdir(directory)

# Run script with function definitions
exec(open('SPICE_Data_Processing_Functions.py').read())

# Ask user for directory where data are located
directory = input('Enter path for SPICEcore data: ')
os.chdir(directory)

# Load CFA file
cfa = pd.read_csv('SPICE_final_sync_unfiltered_24JULY2019.csv')
# Make sure data are in correct format
cfa = cfa.replace('#NAME?', np.nan)
cfa = cfa.astype('float')

# Load other needed files
volcanic_record = pd.read_excel('Full_final_volcanic_record_7August2019.xlsx')
breaks = pd.read_excel('core_breaks_full.xlsx')
annual_depths = pd.read_excel('SPICEcore_Timescale_4_24_2019.xlsx', sheet_name = 'Depth-Age Scale')
dust_events = pd.read_excel('Dust_Events.xlsx')

#%%
# ------------------------------------------------------------------------------------------------------
#                                           2: ERROR REMOVAL
# ------------------------------------------------------------------------------------------------------

# Get original length of the CFA dataset, so errors can be tracked
original_length = len(cfa)
print('\n\n---------------------------------------------------------------------------------')
print('Filtering errors from liquid conductivity, flow rate, depth, and Abakus data.')
print('Original CFA dataset length:', original_length)

# 1) Remove data reflecting bubbles with liquid conductivity values 

#    Note: Liquid conductivity is listed in the 'ECM' column of the CFA data
#    Do this before NaN'ing a bunch of rows
#    Loop through the data and NaN all rows where slopes indicate bubbles. Remove these rows.
threshold_bubbles = 25
bubbles = 0
for i in range(1, len(cfa['Depth (m)']) - 1):                      
    # Calculate the slope between the liquid conductivity at index i and the points before and after it
    if (cfa['Depth (m)'][i] - cfa['Depth (m)'][i - 1]) == 0 or (cfa['Depth (m)'][i + 1] - cfa['Depth (m)'][i]) == 0: 
        continue # Don't divide by 0
    else:
        slope1 = (cfa['ECM'][i] - cfa['ECM'][i - 1]) / (cfa['Depth (m)'][i] - cfa['Depth (m)'][i - 1])
        # No need to do the other calculations if the slope with the point before is above threshold
        if slope1 <= -threshold_bubbles:  
            slope2 = (cfa['ECM'][i + 1] - cfa['ECM'][i]) / (cfa['Depth (m)'][i + 1] - cfa['Depth (m)'][i]) 
            if slope2 >= threshold_bubbles: # If slope1 <= 0 and slope2 >= 0, it's a bubble                       
                bubbles = bubbles + 1
                cfa.loc[i] = np.nan # Change all values in row to NaN
                
print('\n\tBubble errors:               ', bubbles)
# Update dataset length
length = original_length - bubbles

# 2) Filter out data with liquid conductivity values < 0.6

# Drop all good rows, where liquid conductivity is >= 0.6, and count remaining bad rows
# Drop NaN values from bubble analysis first
bad_ecm = cfa['ECM'].dropna()
# Find rows with bad liquid conductivity values-- will NaN the good values
bad_ecm = bad_ecm.where(bad_ecm < 0.6)
# Get the indices of the non-NaN (a.k.a., the bad) rows
bad_rows = bad_ecm[bad_ecm.notnull()]
bad_rows = list(bad_rows.index.values)
# Change values in the bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tLiquid conductivity < 0.6:   ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# 3) Filter out data without positive flow rate values

# Drop all good rows, where flow rate isn't > 0, and count remaining bad rows
bad_flow = cfa['Flow Rate'].dropna()
# Find rows with bad flow rate values-- will NaN the good values
bad_flow = bad_flow.where(bad_flow <= 0)
# Get the indices of the non-NaN (a.k.a., the bad) rows
bad_rows = bad_flow[bad_flow.notnull()]
bad_rows = list(bad_rows.index.values)
# Change all values in the bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tNo/negative flow rate errors:', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# 4) Filter out rows where depth does not increase, rows with no depth value, 
#    and adjacent rows which have discontinous depth values

# Get a copy of the CFA depth values, dropping NaNs
depth_diff = cfa.loc[:, 'Depth (m)'].dropna()
# Subtract each depth value from the depth value before
depth_diff = depth_diff.diff(periods = 1)
# Drop the first row, which becomes NaN
depth_diff = depth_diff.dropna()
# Drop all good rows, where the diff is > 0
bad_depth = depth_diff.drop(depth_diff[depth_diff > 0].index)
# Get a list of all of the indices with bad depth measurements
bad_rows = list(bad_depth.index.values)
# Change values in the bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tDepth not increasing errors: ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# Make sure all NaN'd depths have NaN'd CFA data
depth_isnull = cfa['Depth (m)'].isnull()
# Select rows with NaN depth values
null_depth = depth_isnull[depth_isnull == True]
# Get the indices of these rows
null_depth = list(null_depth.index.values)
# Get the rows with NaN depth values, but non-NaN Abakus values
bad_rows = cfa.loc[null_depth, :].dropna(how = 'all')
# Convert to indices
bad_rows = list(bad_rows.index.values)
# Change values in the bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tRows without depth data:     ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# Remove the first row following a large gap in depth values
# Replace this row with NaN so the data are not plotted as continuous data across the interval

# Subtract each depth value from the one before
depth_diff = cfa.loc[:, 'Depth (m)'].diff(periods = 1)
# Get rows where difference in depth is >= 6 cm
# Don't drop the NaNs. Want to select only adjacent rows where depth jumps by 6 cm
bad_depth = depth_diff[(depth_diff >= 0.06)]
# Get the indices of these rows
bad_rows = list(bad_depth.index.values)
# Change the non-boolean values in these rows to NaN
cfa.loc[bad_rows, :]  = np.nan

print('\n\tRemoving ' + str(len(bad_rows)) + ' rows immediately following depth discontinuities.')
# Update dataset length
length = length - len(bad_rows)

# 5) Filter out any infinite or negative Abakus values

print('\tRemoving negative Abakus values.')
# Get rid of any individual inf. values (just in case)
cfa = cfa.replace([np.inf, -np.inf], np.nan)

# Get a copy of just the Abakus columns
abakus = cfa.loc[:,'1':'12'].copy()
# NaN all negative values (just the values, not the entire row)
abakus[abakus < 0] = np.nan
# Replace CFA abakus columns with corrected abakus columns
cfa.loc[:, '1':'12'] = abakus

# 6) Apply correction to Abakus data from 7/19/2016

print('\tCorrecting units for one melt day.')

cfa = correct_meltday(cfa)

# 7) Interpolate ages for the CFA rows

# Need to do this before adding in the volcanic events
print('\nInterpolating timescale.')

#Interpolate ages for SPICEcore timescale
years_interp = pd.Series(np.interp(cfa['Depth (m)'] ,annual_depths['Depth (m)'],annual_depths['Age (Years Before 1950)'] ))
cfa['AgeBP'] = years_interp

# 8) Label each CFA row near core breaks

print('Labelling core breaks.')

# Add Y/N 'Break?' column. Default to False.
cfa['Break?']     = False
# Add Y'N 'New Break?' column to record first row in each discrete core break range. Default False.
cfa['New Break?'] = False

# Get the row indices of all measurements near core breaks
# Inputs: CFA data, core break data, depth buffer around core breaks
# Buffer: +/- 3 cm of a core break. To change, edit number in the function call below.
break_rows, new_break_rows = label_core_breaks(cfa, breaks, 0.03)
# Change all 'Break?' values in those rows to True
cfa.loc[break_rows, 'Break?']         = True
cfa.loc[new_break_rows, 'New Break?'] = True

# 9) Label all measurements near volcanic events and dust events

print('Labelling volcanic events.')

# Create Y/N 'Volcanic Event?' column. Default to False
cfa['Volcanic Event?']     = False
# This column will indicate the first measurement for each event, as a way to count the events
cfa['New Volcanic Event?'] = False

# Get list of all indices occurring near volcanic events (by year, not depth)
# Function inputs: CFA data, volcanic record, + year buffer, - year buffer
# Buffers: -6/+2 years. To change, edit numbers in the function call below.
volc_rows, new_event_rows = label_volc_events(cfa, volcanic_record, 2, 6)
# Change all 'Volcanic Event?' values in those rows to True
cfa.loc[volc_rows, 'Volcanic Event?']          = True
cfa.loc[new_event_rows, 'New Volcanic Event?'] = True

print('Labelling dust events.')

# Add Y/N 'Dust Event?' column. Default to false.
cfa['Dust Event?'] = False
# Get the row indices of all measurements within dust events
dust_rows = label_dust_events(cfa, dust_events)
# Change all 'Dust Event?' values in those rows to True
cfa.loc[dust_rows, 'Dust Event?'] = True

# 10) Calculate particle concentration and CPP

print('Calculating particle concentration and CPP.')
# Select the rows that aren't completely NaN'd
valid_depth = depth_isnull[depth_isnull == False]
# Get indices of these rows
valid_depth = list(valid_depth.index.values)
# Select the non-NaN CFA rows for concentration (ignoring bin 1)
conc = cfa.loc[valid_depth, '1.1':'12'].copy()
# Get particle sums for these rows. Use skipna = True so NaN cells count as 0.
conc = conc.sum(axis = 1, skipna = True)
# Add concentration column to CFA dataframe
cfa['Sum 1.1-12'] = conc

# Add CPP column to CFA dataframe
cfa['CPP'] = find_cpp(cfa)

# 10) Export CFA file to CSV. Report final length.
print('\nFinished CFA error removal')
print('\tFinal CFA dataset length:', length)
cfa.to_csv('../Data/Cleaned_CFA_Phase1_' + str(date.today()) + '.csv')
      
print('\tData exported to CSV.')
print('---------------------------------------------------------------------------------')





