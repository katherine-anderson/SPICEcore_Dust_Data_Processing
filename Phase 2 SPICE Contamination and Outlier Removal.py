# --------------------------------------------------------------------------------------
#                     SPICE REMOVE OUTLIERS AND CONTAMINATION SCRIPT
#               Removes anomalies and outliers from the error-free CFA data
#
#    - Reads CFA dataframe (with mechanical errors removed)
#    - Counts and NaNs all rows with a 'hump-shaped' PSD anomaly
#    - Removes outliers using MAD
#    - Removes outliers using integrals
#    - Preserves known dust events during 'hump' and outlier removal
#    - Prints summary statistics
#    - Saves cleaned and bad data to two separate files
#
# Aaron Chesler and Katie Anderson, 8/8/19
# ---------------------------------------------------------------------------------------
#%%
# ---------------------------------------------------------------------------------------
#                               1: FILE PREPARATION
# ---------------------------------------------------------------------------------------

# Import modules and packages
import numpy  as np
from   numpy import trapz
import pandas as pd
import os
from   datetime import date

# Ask user whether to run phase 1 cleaning (error removal) from this script
choice = input('Run SPICEcore error removal script from this file? Enter Y or N: ')
if choice == 'y' or choice == 'Y':
    exec(open('cleaning_cfa.py').read())

# Ask user for directory where data are located
else: 
    # Run function definitions script. It's included in the phase 1 script
    exec(open('SPICE_Data_Processing_Functions.py').read())
    # Get directory for data files
    directory = input('Enter path for SPICEcore data: ')
    os.chdir(directory)

# Load complete CFA file after mechanical error removal
# Ask user for CFA file to use
file = input('Enter name of the CFA file after error removal as .csv: ')
cfa_phase1 = pd.read_csv(file, header = 0)
del cfa_phase1['Unnamed: 0']
# Make separate copies of the CFA data before and after phase 2 cleaning
cfa = cfa_phase1.copy()

# Get the row indices of all measurements within dust events
dust_rows = cfa[(cfa['Dust Event?'] == True)].index.values.tolist()

# Get the row indices of all measurements within volcanic events
volc_rows = cfa[(cfa['Volcanic Event?'] == True)].index.values.tolist()

#%%
# ---------------------------------------------------------------------------------------
#                            PART 2: Outlier and Contamination Removal
# ---------------------------------------------------------------------------------------

print('\n\n-----------------------------------------------------------------------')
length = 438212 # From error removal script. Update as needed.
print('\n\nRemoving outliers and contamination signals.')
print('CFA dataset length after error removal:', length)

# 1) Identify and remove hump-shaped PSD anomalies

# Find humps for all CFA depths
# Inputs: CFA data, minimum depth, maximum depth. Currently using the full core.
humps = find_humps(cfa, 0, 1752)

# Remove all rows in real dust events from the hump list
bad_rows = humps.index.difference(dust_rows)
# Remove all rows in real volcanic events from the hump list
bad_rows = humps.index.difference(volc_rows)

# Save all bad data into a separate dataframe
bad_cfa = cfa.loc[bad_rows, :].copy()
bad_cfa['Error Type'] = 'PSD Hump Anomaly'

# NaN values in the bad rows in the original CFA data, except boolean columns
cfa.loc[bad_rows, ['Depth (m)', 'AgeBP', 'Flow Rate', 'ECM', '1', '1.1', '1.2', 
                   '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9', '2', '2.1', 
                   '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', '3.2', '3.6', '4', 
                   '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12',
                   'CPP', 'Sum 1.1-12']] = np.nan

# Print number of measurements removed
print('\tPSD anomalies removed: ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# 2) Identify and remove particle concentration & CPP outliers, using MAD

# Set # of measurements to use for background medians
window = 500
# Set threshold for accepted Median Absolute Deviations (MAD) (e.g., 2 * MAD)
threshold = 2

# Remove overlapping concentration & CPP outliers
# Inputs: CFA data, dust event indices, volcanic event indices, background window size, and MAD threshold
bad_rows = remove_outliers_MAD(cfa, dust_rows, volc_rows, window, threshold)

# Add bad data to the bad CFA dataframe
bad_cfa = bad_cfa.append(cfa.loc[bad_rows, :], sort = False)
# Label error type
bad_cfa['Error Type'].fillna('MAD Outlier', inplace = True)

# NaN values in the bad rows in the original CFA data, except boolean columns
cfa.loc[bad_rows, ['Depth (m)', 'AgeBP', 'Flow Rate', 'ECM', '1', '1.1', '1.2', 
                   '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9', '2', '2.1', 
                   '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', '3.2', '3.6', '4', 
                   '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12',
                   'CPP', 'Sum 1.1-12']] = np.nan

print('\tMAD outliers removed: ', len(bad_rows))

# Update dataset length
length = length - len(bad_rows)

# 3) Identify and remove particle concentration & CPP outliers, using 2-pt. integrals

# UPDATE THIS FUNCTION WITH AARON'S NEWEST VERSION

outlier_rows, bad_rows = remove_outliers_integrals(new_cfa, 2, dust_rows, volc_rows)

# Add bad data to the bad CFA dataframe
bad_cfa = bad_cfa.append(cfa.loc[bad_rows, :], sort = False)
# Label error type
bad_cfa['Error Type'].fillna('Integral Outlier', inplace = True)

# NaN values in remaining rows, except boolean columns
cfa.loc[bad_rows, ['Depth (m)', 'AgeBP', 'Flow Rate', 'ECM', '1', '1.1', '1.2', 
                   '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9', '2', '2.1', 
                   '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', '3.2', '3.6', '4', 
                   '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12',
                   'CPP', 'Sum 1.1-12']] = np.nan

print('\tIntegral outliers removed:', len(bad_rows)

# Update dataset length
length = length - len(bad_rows)

# ADD DATA INTERPOLATION, TAIL CORRECTIONS HERE?

# 4) Compute summary statistics before and after phase 2 cleaning, if requested

choice = input('Print summary statistics? Enter Y or N: ')
if choice == 'Y' or choice == 'y':

    print('\n--Phase 1 Cleaning Results--')
    # Input the before & after CFA data into the summary statistics function
    summary_statistics(cfa_phase1)
    print('\n--Phase 2 Cleaning Results--')
    summary_statistics(cfa)

# 5) Export CFA file to CSV. Report final length.
print('\n\nFinished CFA outlier & contamination removal')
print('\n\tFinal CFA dataset length:', length)

cfa.to_csv('Cleaned_CFA_Phase2_' + str(date.today()) + '.csv')
bad_cfa.to_csv('Bad_CFA_Phase2_' + str(date.today()) + '.csv')
print('\n\tData exported to CSV. Bad data saved in separate file.')
print('-----------------------------------------------------------------------')





