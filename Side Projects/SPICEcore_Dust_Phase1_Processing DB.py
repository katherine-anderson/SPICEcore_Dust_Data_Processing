# ------------------------------------------------------------------------------------------------------
#                     SPICEcore Dust Phase 1 Processing Script: Mechanical Error Removal
# Removes melting errors in the raw CFA data, interpolates a timescale, and adds descriptive columns
#
#    - Dust processing has 2 phases
#      1) Mechanical error removal 
#      2) Removal of outliers & anomalies
# 
# Phase 1 Dust Processing
#    - Loads raw, unfiltered continuous flow analysis (CFA) data with minor depth corrections
#    - Loads supporting datafiles
#    - Tracks the number of measurements NaN'ed in each step
#
#    1) NaNs data from air bubbles bubbles using liquid conductivity values
#    2) NaNs liquid conductivity values < 0.6 us
#    3) NaNs measurements without positive flow rates
#    4) NaNs measurements with depth duplicates or decreases
#    5) NaNs measurements with infinite or negative dust values
#    6) Applies correction to dust data from bad melt day (7/19/2016)
#    7) Applies timescale to the dust data (annual layers in the Holocene, volcanic tie points for the glacial)
#    8) Labels all measurements near core breaks
#    9) Labels all measurements within volcanic events and dust events
#   10) Calculates particle concentration and coarse particle percentage (CPP)
#   11) Exports cleaned dataset to CSV
#
# Katie Anderson and Aaron Chesler, 7/16/20
# ------------------------------------------------------------------------------------------------------
#%%
# ------------------------------------------------------------------------------------------------------
#                                           1: FILE PREPARATION
# ------------------------------------------------------------------------------------------------------
print('\n\n.......................................................')
print('  SPICEcore Dust Data Phase 1 Cleaning: Melter Errors')
print('.......................................................')

# Import needed modules & packages
import pandas as pd
import os 
import numpy  as np
from datetime import date

# Ask user for directory where scripts are located
directory = input('Enter path for SPICEcore dust scripts: ')
os.chdir(directory)

# Run script with function definitions
exec(open('SPICEcore_Dust_Processing_Functions.py').read())

# Ask user for directory where data are located
directory = input('Enter path for SPICEcore dust data: ')
os.chdir(directory)

# Load CSV CFA data as floats
cfa = pd.read_csv('CFA_Unfiltered_Synchronized_1_2_20.csv', dtype = 'float', index_col = 'Unnamed: 0')

# Load other needed files
volcanic_record = pd.read_excel('Full_final_volcanic_record_7August2019.xlsx')
breaks = pd.read_excel('core_breaks_full.xlsx')
annual_depths = pd.read_excel('SPICEcore_Timescale_4_24_2019.xlsx', sheet_name = 'Depth-Age Scale')
dust_events = pd.read_excel('Dust_Events.xlsx')

# Interpolate ages for glacial volcanic events
years_interp = pd.Series(np.interp(volcanic_record['Volcanic Depth (m)'], annual_depths['Depth (m)'], annual_depths['Age (Years Before 1950)']))
volcanic_record.loc[1209:, 'Start Year (b1950)'] = years_interp

#%%
# ------------------------------------------------------------------------------------------------------
#                                           2: ERROR REMOVAL
# ------------------------------------------------------------------------------------------------------

# Get original length of the CFA dataset, so errors can be tracked
original_length = cfa['1'].count()
print('\n\n---------------------------------------------------------------------------------')
print('Filtering errors from liquid conductivity, flow rate, depth, and Abakus data.')
print('Original CFA dataset length:', original_length)

bad_rows_master = pd.DataFrame(columns = ['Row', 'Error Type'])

# 1) Remove data reflecting bubbles with liquid conductivity values 

#    Note: Liquid conductivity is listed in the 'ECM' column of the CFA data
#    Do this before NaN'ing a bunch of rows
#    Loop through the data and NaN all rows where slopes indicate bubbles. Remove these rows.

error = 'Bubble Error'
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
            if slope2 >= threshold_bubbles: # If slope1 <= -25 and slope2 >= 25, it's a bubble                       
                bubbles = bubbles + 1
                cfa.loc[i, :] = np.nan # Change all values in row to NaN
                bad_rows_master = bad_rows_master.append({'Row': i, 'Error Type': 'Liquid Conductivity Error'}, ignore_index = True)
                # bad_rows_master['Error Type'].append('Bubble Error')
                
print('\n\tBubble errors:               ', bubbles)
# Update dataset length
length = original_length - bubbles

#%%
# 2) Filter out data with liquid conductivity values < 0.6

# Get bad rows
bad_rows = cfa[cfa['ECM'] < 0.6]
# Get indices of bad rows
bad_rows = bad_rows.index
# Change values in bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

for i in range(0, len(bad_rows)):
    bad_rows_master = bad_rows_master.append({'Row': bad_rows[i], 'Error Type': 'Liquid Conductivity Error'}, ignore_index = True)


print('\tLiquid conductivity < 0.6:   ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)
#%%
# 3) Filter out data without positive flow rate values

# Get bad rows
bad_rows = cfa[cfa['Flow Rate'] <= 0]
# Get indices of bad rows
bad_rows = bad_rows.index
# Change values in bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

for i in range(0, len(bad_rows)):
    bad_rows_master = bad_rows_master.append({'Row': bad_rows[i], 'Error Type': 'Flow Rate Error'}, ignore_index = True)


print('\tNo/negative flow rate errors:', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# 4) Filter out rows where depth does not increase and rows with no depth value

# Select data by depth, drop all rows with NaN depth values
depth_diff = cfa.loc[:, 'Depth (m)'].dropna()
# Subtract each depth value from the depth value in the row above
depth_diff = depth_diff.diff(periods = 1)
# Drop the first row, which becomes NaN
depth_diff = depth_diff.dropna()
# Drop all good rows, where the difference in depth is > 0
bad_depth = depth_diff.drop(depth_diff[depth_diff > 0].index)
# Get a list of all of the indices with bad depth measurements
bad_depth = list(bad_depth.index.values)
# Exclude the rows which have already been NaN'd from the depth error counts
# Select 'bad depth' rows with non-NaN flow rates
bad_rows = cfa.loc[bad_depth, 'Flow Rate'].dropna()
# Get a list of all of the indices with only bad depth measurements
bad_rows = list(bad_rows.index.values)
# Change ALL values in the bad depth rows to NaN
cfa.loc[bad_depth, :] = np.nan

for i in range(0, len(bad_rows)):
    bad_rows_master = bad_rows_master.append({'Row': bad_rows[i], 'Error Type': 'Depth Error'}, ignore_index = True)


print('\tDepth not increasing errors: ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# Make sure all NaN'd depths have NaN'd CFA data
depth_isnull = cfa['Depth (m)'].isnull()
# Select rows with NaN depth values
null_depth = depth_isnull[depth_isnull == True]
# Get the indices of these rows
null_depth = list(null_depth.index.values)
# Make sure to select NaN depth values and non-NaN Abakus values
bad_rows = cfa.loc[null_depth, :].dropna(how = 'all')
# Convert to indices
bad_rows = list(bad_rows.index.values)
# Change values in the bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

for i in range(0, len(bad_rows)):
    bad_rows_master = bad_rows_master.append({'Row': bad_rows[i], 'Error Type': 'Depth Error'}, ignore_index = True)

print('\tRows without depth data:     ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# 5) Filter out any infinite or negative Abakus values

# Get indices of rows with infs
inf_rows = cfa.index[np.isinf(cfa.loc[:,'1':'12']).any(1)]
# Change values in the bad rows to NaN
cfa.loc[inf_rows, :] = np.nan
# Update dataset length 
length = length - len(inf_rows)

for i in range(0, len(inf_rows)):
    bad_rows_master = bad_rows_master.append({'Row': inf_rows[i], 'Error Type': 'Abakus Error'}, ignore_index = True)


# Select rows where any of the Abakus values are negative
bad_rows = cfa[(cfa.loc[:, '1':'12'] < 0).any(1) == True]
# Get indices of these rows
bad_rows = bad_rows.index
# Change values in the bad rows to NaN
cfa.loc[bad_rows, :] = np.nan

for i in range(0, len(bad_rows)):
    bad_rows_master = bad_rows_master.append({'Row': bad_rows[i], 'Error Type': 'Abakus Error'}, ignore_index = True)

# Update dataset length
length = length - len(bad_rows)
print('\tRows with invalid dust data: ', len(bad_rows) + len(inf_rows))

# 6) Apply correction to Abakus data from 7/19/2016

print('\tCorrecting units for one melt day.')

cfa = correct_meltday(cfa)

# 7) Interpolate ages for the CFA rows

# Need to interpolate ages before adding in the volcanic events
print('Interpolating depth-age timescale.')

#Interpolate ages for SPICEcore timescale
years_interp = pd.Series(np.interp(cfa['Depth (m)'], annual_depths['Depth (m)'], annual_depths['Age (Years Before 1950)']))
cfa['AgeBP'] = years_interp
#%%
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
# Need at least 1 value to sum (skip NaN rows)
cfa['Sum 1.1-12'] = cfa.loc[:, '1.1':'12'].sum(axis = 1, min_count = 1)

# Add CPP column to CFA dataframe
cfa['CPP'] = find_cpp(cfa)
#%%
# 10) Export CFA file to CSV. Report final length.
print('\nFinished Phase 1 dust processing.')
print('\tFinal dataset length:', length)
cfa.to_csv('Cleaned_CFA_Phase1_' + str(date.today()) + '.csv')
      
print('\tData exported to CSV [Cleaned_CFA_Phase1_...].')
print('---------------------------------------------------------------------------------')
