#!/usr/bin/env python
# coding: utf-8

# In[1]:


# --------------------------------------------------------------------------------------
#                        SPICE REMOVE ERRORS SCRIPT 3.0
# Removes MEASUREMENT errors in the raw CFA data, labels CORE BREAKS & VOLCANIC EVENTS
#
#    - Loads raw, unfiltered CFA with depth corrections
#    - NaNs bubbles and liquid conductivity values < 0.6
#    - NaNs measurements without positive flow rates
#    - NaNs negative Abakus values
#    - NaNs measurements with depth duplicates or decreases
#    - Tracks the number of measurements NaN'ed in each step
#    - Labels all measurements within specified depth ranges of core breaks
#    - Creates timescale for CFA data (annual in Holocene, tie points for glacial)
#    - Labels all measurements within specific years of volcanic events
#    - Labels the starting row for each volcanic event
#    - Exports cleaned dataset to CSV
#
# Katie Anderson, 7/10/19
# ---------------------------------------------------------------------------------------


# In[2]:


from   scipy.io import loadmat
from   scipy    import interpolate
import numpy  as np
import pandas as pd
import csv
import os
from   datetime import date
import matplotlib.pyplot as plt
import statistics

# Go into the correct directory- this is a personal desktop folder that backs up to OneDrive
# Current setup: Master SPICE folder with scripts, data, and figure subfolders
os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\SPICE\\Scripts')

# Run script to establish function definitions
get_ipython().run_line_magic('run', '"SPICE Data Processing Functions.ipynb"')


# In[5]:


# --------------------CFA FILE PREP--------------------------#

# Load Matlab raw, unfiltered, depth-corrected CFA file
#mat = loadmat('../Data/CFA_Unfiltered_Synchronized.mat')
# Import main variable from Matlab file
#mdata = mat['FinalCFA']

# Original files I've been using. Give different results than the unfiltered_synchronized file...
#mat = loadmat('../Data/CFA_FullCore_Unfiltered.mat') 
mat = loadmat('../Data/CFA_Unfiltered_Synchronized_7_24_19.mat')
mdata = mat['FinalCFA']

# Create dataframe and add 1st column
cfa = pd.DataFrame({'Depth (m)':mdata[:,0]})
# Add remaining columns and data to the dataframe
cfa['Flow Rate'] = mdata[:,1]
cfa['ECM'] = mdata[:,2]
cfa['1']   = mdata[:,3]
cfa['1.1'] = mdata[:,4]
cfa['1.2'] = mdata[:,5]
cfa['1.3'] = mdata[:,6]
cfa['1.4'] = mdata[:,7]
cfa['1.5'] = mdata[:,8]
cfa['1.6'] = mdata[:,9]
cfa['1.7'] = mdata[:,10]
cfa['1.8'] = mdata[:,11]
cfa['1.9'] = mdata[:,12]
cfa['2']   = mdata[:,13]
cfa['2.1'] = mdata[:,14]
cfa['2.2'] = mdata[:,15]
cfa['2.3'] = mdata[:,16]
cfa['2.4'] = mdata[:,17]
cfa['2.5'] = mdata[:,18]
cfa['2.7'] = mdata[:,19]
cfa['2.9'] = mdata[:,20]
cfa['3.2'] = mdata[:,21]
cfa['3.6'] = mdata[:,22]
cfa['4']   = mdata[:,23]
cfa['4.5'] = mdata[:,24]
cfa['5.1'] = mdata[:,25]
cfa['5.7'] = mdata[:,26]
cfa['6.4'] = mdata[:,27]
cfa['7.2'] = mdata[:,28]
cfa['8.1'] = mdata[:,29]
cfa['9']   = mdata[:,30]
cfa['10']  = mdata[:,31]
cfa['12']  = mdata[:,32]

# Load corrected core breaks file
breaks = pd.read_csv('../Data/Core Breaks Full Core.csv')

# Load annual timescale for the Holocene
holocene_timescale = pd.read_excel('../Data/Holocene Timescale.xlsx')

# Load tie point timescale for the glacial
glacial_timescale = pd.read_excel('../Data/Glacial Tie Point Timescale.xlsx')

# Load combined Holocene and glacial volcanic records

# These are the 3*MAD from Dave, 101-pt running median
volcanic_record = pd.read_excel('../Data/Volcanic Record.xlsx')
# Select the glacial rows, which don't yet have ages
glacial_volc = volcanic_record[volcanic_record['Volcanic Depth (m)'].notnull()].copy()
# Interpolate ages for the glacial volcanic events
glacial_volc_age_interp = pd.Series(np.interp(glacial_volc['Volcanic Depth (m)'], 
                                              glacial_timescale['Bot D (m)'], glacial_timescale['Bot Year (b1950)']))
glacial_volc['Start Year (b1950)'] = glacial_volc_age_interp.values
# Add glacial volcanic ages to the complete volcanic record list
volcanic_record.loc[1209:5400, 'Start Year (b1950)'] = glacial_volc['Start Year (b1950)']


# In[6]:


# ----------------------------------------------------------------------------
#                                        PART 1:
#                   CFA DATA FILTERING AND MECHANICAL ERROR REMOVAL
# ----------------------------------------------------------------------------

original_length = len(cfa)
print('Original CFA dataset length:', original_length)
print('\nFiltering liquid conductivity, flow rate, Abakus, and depth data.\n')

# 1) Remove data reflecting bubbles with ECM 
#    Do this before NaN'ing a bunch of rows
#    DOM AND AARON DEFINE BUBBLES DIFFERENTLY. DOM: < 90% OF LAST/NEXT. AARON: STEEP +/- SLOPES (25)

#    Loop through the data and NaN all rows where slopes indicate bubbles. Remove these rows.
threshold_bubbles = 25
bubbles = 0
for i in range(1, len(cfa['Depth (m)']) - 1):                      
    # Calculate the slope between the ECM at index i and the points before and after it
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
print('Bubble errors:               ', bubbles)

# 2) Filter out data with ECM values < 0.6

# Drop all good rows, where ECM is >= 0.6, and count remaining bad rows
bad_ecm = cfa.drop(cfa[cfa['ECM'] >= 0.6].index)
ecm = len(bad_ecm)

# Change all values in rows with ECM < 0.6 to NaN
bad_rows = list(bad_ecm.index.values)
# Change values in the bad depth rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('Liquid conductivity < 0.6:   ', ecm - bubbles)

# 3) Filter out data without positive flow rate values

# Drop all good rows, where flow rate isn't > 0, and count remaining bad rows
bad_flow = cfa.drop(cfa[cfa['Flow Rate'] > 0].index)
flow = len(bad_flow)
# Change all values in rows without positive flow rates to NaN
bad_rows = list(bad_flow.index.values)
# Change values in the bad depth rows to NaN
cfa.loc[bad_rows, :] = np.nan

print('No/negative flow rate errors:', flow - ecm)

# 4) Filter out measurements where depth does not increase

# Subtract each row from the one before
# All NaNs stay NaN
diff = cfa.diff(periods = 1, axis = 0)
# Drop all good rows, where the diff is > 0
bad_depth = diff.drop(diff[diff['Depth (m)'] > 0].index)
# Delete all rows with NaN and count remaining bad rows
bad_depth = bad_depth.dropna()
depth = len(bad_depth)

print('Depth not increasing errors: ', depth)

# Get a list of all of the indices with bad depth measurements
bad_rows = list(bad_depth.index.values)
# Change values in the bad depth rows to NaN
cfa.loc[bad_rows, :] = np.nan

# 5) Filter out any inf. or negative Abakus values. Check that everything is NaN'd.
print('\nRemoving negative Abakus values.')
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

# 6) Label each CFA row near core breaks
print('\nLabelling core breaks.')

# Add Y/N 'Break?' column. Default to False.
cfa['Break?'] = False
# Add Y'N 'New Break?' column to record first row in each discrete core break range. Default False.
cfa['New Break?'] = False

# Get the row indices of all measurements near core breaks
# Inputs: CFA data, break list, depth range
break_rows, new_break_rows = label_core_breaks(cfa, breaks, 0.03)
# Change all 'Break?' values in those rows to True
cfa.loc[break_rows, 'Break?'] = True
cfa.loc[new_break_rows, 'New Break?'] = True

#7) Interpolate ages for the CFA rows
# Need to do this before adding in the volcanic events
print('\nInterpolating timescale.')

# Subset the Holocene CFA data
# This is the row where depths < 798 m (closest to end of annual timescale)
# Make a copy of this data to avoid chained assignment later
cfa_holocene = cfa.loc[0:214959, :].copy()

# Interpolate ages for Holocene annual timescale
age_interp = pd.Series(np.interp(cfa_holocene['Depth (m)'], 
                                 holocene_timescale['Depth (m)'], holocene_timescale['Age (yr b 1950)']))
# Save ages to copy of Holocene CFA. Keep until step 8) is complete for the whole core
cfa_holocene['Age b 1950'] = age_interp.values

# Subset the deep CFA data
# Start with the next row after the Holocene rows
cfa_deep = cfa.loc[214960:448667, :].copy()

# Interpolate ages for the deep core tie points
glacial_age_interp = pd.Series(np.interp(cfa_deep['Depth (m)'],
                                         glacial_timescale['Bot D (m)'], glacial_timescale['Bot Year (b1950)']))

# Create one series with the interpolated age for each row
ages = age_interp.append(glacial_age_interp)
# Add 'Age' column with age values
cfa['Age b 1950'] = ages.values

#8) Label all measurements near volcanic events (currently, only for the Holocene)

print('\nLabelling volcanic events.')

# Create Y/N 'Volcanic Event?' column. Default to False
cfa['Volcanic Event?'] = False
cfa['New Volcanic Event?'] = False

# Get list of all indices occurring near volcanic events (by year, not depth)
# Minus 2- and plus 2-year buffers
volc_rows, new_event_rows = label_volc_events(cfa, volcanic_record, 2, 6)
# Change all 'Volcanic Event?' values in those rows to True
cfa.loc[volc_rows, 'Volcanic Event?'] = True
cfa.loc[new_event_rows, 'New Volcanic Event?'] = True

# 9) Export CFA files to CSV. Report final length.
print('\nFinal CFA dataset length:', (len(cfa) - flow - depth))
cfa.to_csv('../Data/Cleaned_CFA_Phase1_' + str(date.today()) + '.csv')
      
print('\nData exported to CSV.')


# In[ ]:




