# --------------------------------------------------------------------------------------
#                        SPICE REMOVE CONTAMINATION SCRIPT
#               Removes anomalies and outliers from the error-free CFA data
#
#    - Reads CFA dataframe (with mechanical errors removed)
#    - Calculates CPP and particle concentration (bins 1.1-10) for the full core
#    - Counts and NaNs all rows with a 'hump-shaped' PSD anomaly
#    - Creates MAD outlier removal scenarios
#    - Removes integral outliers (after Aaron)
#    - Preserves known dust events during 'hump' and outlier removal
#
#    - Loads raw, unfiltered CFA with depth corrections. Gets CPP and particle conc.
#    - Creates a series of plots to evaluate different outlier removal scenarios
#    - Creates before & after plots for the data processing
#
# Aaron Chesler and Katie Anderson, 8/8/19
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
#                               1: FILE PREPARATION
# ---------------------------------------------------------------------------------------

# Import modules and packages
from   scipy.io import loadmat
import numpy  as np
from   numpy import trapz
import pandas as pd
import csv
import os
from   datetime import date
import matplotlib.pyplot as plt
import statistics

# Run function definitions script
exec(open('SPICE_Data_Processing_Functions.py').read())

# Ask user whether to run phase 1 cleaning (error removal) from this script
choice = input('Run SPICEcore error removal script from this file? Enter Y or N: ')
if choice == 'y' or choice == 'Y':
    exec(open('cleaning_cfa.py').read())
else: continue

# Ask user for directory where data are located
directory = input('Enter path for SPICEcore data: ')
os.chdir(directory)

# Load complete CFA file after mechanical error removal (-2/+6 yr volcanic buffer)
cfa = pd.read_csv('Cleaned_CFA_Phase1_2019-08-07.csv', header = 0)
del cfa['Unnamed: 0']

# Get the row indices of all measurements within dust events
dust_rows = cfa[(cfa['Dust Event?'] == True)].index.values.tolist()

# Get the row indices of all measurements within volcanic events
volc_rows = cfa[(cfa['Volcanic Event?'] == True)].index.values.tolist()



# ---------------------------------------------------------------------------------------
#                            PART 2: Outlier and Contamination Removal
# ---------------------------------------------------------------------------------------

original_length = 438212 # From error removal script. Update as needed.
print('CFA dataset length after error removal:', original_length)
print('Removing outliers and contamination signals.')

# 1) Identify and remove hump-shaped PSD anomalies

# Find humps for all CFA depths
# Inputs: CFA data, minimum depth, maximum depth. Currently using the full core.
humps = find_humps(cfa, 0, 1752)

# PRESERVE KNOWN DUST EVENTS
# Remove all rows in real dust events from the hump list
bad_rows = humps.index.difference(dust_rows)
# Remove all rows in real volcanic events from the hump list
bad_rows = humps.index.difference(volc_rows)

# NaN values in remaining rows, except boolean columns
contamination = pd.DataFrame()
cfa.loc[bad_rows, 'Depth (m)': '12'] = np.nan
cfa.loc[bad_rows, 'AgeBP'] = np.nan
cfa.loc[bad_rows, 'CPP']        = np.nan

# Print number of measurements removed
print('\tPSD anomalies removed: ', len(bad_rows))
# Update dataset length
original_length = original_length - len(bad_rows)

# 2) Identify and remove particle concentration & CPP outliers, using MAD

# Set # of measurements to use for background medians
window = 500
# Set threshold for accepted Median Absolute Deviations (MAD).
threshold = 2
# Make a new copy of the CFA data to ensure that the data cleaning doesn't change original data
new_cfa = cfa.copy()

# Remove overlapping concentration & CPP outliers
# Inputs: CFA data, dust event indices, volcanic event indices, background window size, and MAD threshold
new_cfa, num_outliers = remove_outliers_MAD(new_cfa, dust_rows, volc_rows, window, threshold)

print('\tMAD outliers removed: ', num_outliers)

# 3) Identify and remove particle concentration & CPP outliers, using 2-pt. integrals

outlier_rows, bad_rows = remove_outliers_integrals(new_cfa, 2, dust_rows, volc_rows)

# NaN values in remaining rows, except boolean columns
new_cfa.loc[bad_rows, 'Depth (m)':'12']   = np.nan
new_cfa.loc[bad_rows, 'AgeBP']       = np.nan
new_cfa.loc[bad_rows, 'Sum 1.1-12':'CPP'] = np.nan

print('\tIntegral outliers removed:', len(bad_rows))



# 4) Compute summary statistics before and after phase 2 cleaning

print('\n--Phase 1 Cleaning Results--')
# Input the before & after CFA data into the summary statistics function
summary_statistics(cfa)
print('\n--Phase 2 Cleaning Results--')
summary_statistics(new_cfa)

# 5) Export CFA file to CSV. Report final length.

print('\nFinal CFA dataset length:', original_length - num_outliers - len(bad_rows))

cfa.to_csv('../Data/Cleaned_CFA_Phase2_' + str(date.today()) + '.csv')
print('\nData exported to CSV.')





