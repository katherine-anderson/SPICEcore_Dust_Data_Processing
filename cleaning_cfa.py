# ------------------------------------------------------------------------------------------------------
#                                    FINAL SPICE REMOVE ERRORS SCRIPT
# Removes melting errors in the raw CFA data, interpolates timescale, and adds descriptive columns
#
#    - Loads raw, unfiltered CFA with depth corrections
#    - NaNs bubbles and liquid conductivity values < 0.6
#    - NaNs measurements without positive flow rates
#    - NaNs negative Abakus values
#    - NaNs measurements with depth duplicates or decreases
#    - Tracks the number of measurements NaN'ed in each step
#    - Labels all measurements within specified depth ranges of core breaks
#    - Creates timescale for CFA data (annual in Holocene, tie points for glacial)
#    - Labels all measurements within specified years of volcanic events
#    - Labels the starting row for each volcanic event
#    - Labels all measurements within known dust events
#    - Calculates particle concentration and CPP
#    - Exports cleaned dataset to CSV
#
# Katie Anderson and Aaron Chesler, 8/8/19
# ------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------
#                                           1: FILE PREPARATION
# ------------------------------------------------------------------------------------------------------

# Import needed modules & packages
import pandas as pd
import os 
import numpy as np
from datetime import date

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
annual_depths = pd.read_excel('SPICEcore_Timescale_4_24_2019.xlsx', sheetname = 'Depth-Age Scale')
dust_events = pd.read_excel('Dust_Events.xlsx')




# ------------------------------------------------------------------------------------------------------
#                                           2: ERROR REMOVAL
# ------------------------------------------------------------------------------------------------------

# Get original length of the CFA dataset, so errors can be tracked
original_length = len(cfa)
print('\n\n---------------------------------------------------------------------------------')
print('Filtering errors from liquid conductivity, flow rate, Abakus, and depth data.')
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

# 2) Filter out data with liquid conductivity values < 0.6

# Drop all good rows, where liquid conductivity is >= 0.6, and count remaining bad rows
bad_ecm = cfa.drop(cfa[cfa['ECM'] >= 0.6].index)
# Change all values in rows with liquid conductivity < 0.6 to NaN
bad_rows = list(bad_ecm.index.values)
# Change values in the bad depth rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tLiquid conductivity < 0.6:   ', len(bad_ecm) - bubbles)

# 3) Filter out data without positive flow rate values

# Drop all good rows, where flow rate isn't > 0, and count remaining bad rows
bad_flow = cfa.drop(cfa[cfa['Flow Rate'] > 0].index)
flow = len(bad_flow)
# Change all values in rows without positive flow rates to NaN
bad_rows = list(bad_flow.index.values)
# Change values in the bad depth rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tNo/negative flow rate errors:', len(bad_flow) - len(bad_ecm))

# 4) Filter out measurements where depth does not increase

# Subtract each row from the one before
# All NaNs stay NaN
diff = cfa.diff(periods = 1, axis = 0)
# Drop all good rows, where the diff is > 0
bad_depth = diff.drop(diff[diff['Depth (m)'] > 0].index)
# Delete all rows with NaN and count remaining bad rows
bad_depth = bad_depth.dropna()

# Get a list of all of the indices with bad depth measurements
bad_rows = list(bad_depth.index.values)
# Change values in the bad depth rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('\tDepth not increasing errors: ', len(bad_depth))

# 5) Filter out any infinite or negative Abakus values. Check that NaN'd rows are completely NaN'd

print('\tRemoving negative Abakus values.')
# Get rid of any individual inf. values (just in case)
cfa = cfa.replace([np.inf, -np.inf], np.nan)

# Get a copy of just the Abakus columns
abakus = cfa.loc[:,'1':'12'].copy()
# NaN all negative values (just the values, not the entire row)
abakus[abakus < 0] = np.nan
# Replace CFA abakus columns with corrected abakus columns
cfa.loc[:, '1':'12'] = abakus

# Make sure all NaN'd depths have NaN'd CFA data
depth_isnull = cfa['Depth (m)'].isnull()
bad_rows = depth_isnull[depth_isnull == True].index.values
cfa.loc[bad_rows, :] = np.nan

# THIS IS WHERE WE COULD PUT THE 300-M CORRECTION

# 6) Label each CFA row near core breaks

print('\nLabelling core breaks.')

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

#7) Interpolate ages for the CFA rows

# Need to do this before adding in the volcanic events
print('Interpolating timescale.')

#Interpolate ages for SPICEcore timescale
years_interp = pd.Series(np.interp(cfa['Depth (m)'] ,annual_depths['Depth (m)'],annual_depths['Age (Years Before 1950)'] ))
cfa['AgeBP'] = years_interp

#8) Label all measurements near volcanic events and dust events

print('Labelling volcanic events.')

# Create Y/N 'Volcanic Event?' column. Default to False
cfa['Volcanic Event?']     = False
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

# 9) Calculate particle concentration and CPP, add to dataframe

print('Calculating particle concentration and CPP.')
# Change this list if we decide to include smallest & largest bins
sum_columns = ['1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', 
               '1.9', '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.7', '2.9',
               '3.2', '3.6', '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', 
               '9', '10', '12']
# Set skipna to False, otherwise, rows with all NaNs will sum to 0
cfa['Sum 1.1-12'] = cfa[sum_columns].sum(axis = 1, skipna = False)

# Add CPP column to CFA dataframe. Function will ask to include/exclude bins 1 and 12
cfa['CPP'] = find_cpp(cfa)

# 10) Export CFA file to CSV. Report final length.
print('\nFinished CFA error removal')
print('\tFinal CFA dataset length:', (len(cfa) - len(bad_flow) - len(bad_depth)))
cfa.to_csv('../Data/Cleaned_CFA_Phase1_' + str(date.today()) + '.csv')
      
print('\tData exported to CSV.')
print('---------------------------------------------------------------------------------')





