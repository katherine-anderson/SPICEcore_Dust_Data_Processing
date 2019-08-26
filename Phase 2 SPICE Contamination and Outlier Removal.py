# --------------------------------------------------------------------------------------
#                     SPICE REMOVE OUTLIERS AND CONTAMINATION SCRIPT
#               Removes anomalies and outliers from the error-free CFA data
#
#    - Gives user the option to run all of the data cleaning from this script
#    - Cleans anomalies and outliers from the CFA data after mechanical error removal
#      - Preserves data during known dust and volcanic events
#      - Saves 'bad' data into another file, labelled by error type
#      - NaNs 'bad' data in the CFA file and prints error counts
#      - Error types:
#        1) PSD hump anomaly
#        2) MAD outliers
#        3) Integral outliers
#        4) Manually-identified issues which remain
#    - Prints summary statistics
#    - Saves cleaned and 'bad' data to two separate files
#
# Aaron Chesler and Katie Anderson, 8/26/19
# ---------------------------------------------------------------------------------------
#%%
# ---------------------------------------------------------------------------------------
#                               1: FILE PREPARATION
# ---------------------------------------------------------------------------------------
print('\n\n...................................................................')
print('                 SPICECORE DUST DATA PROCESSING')
print('...................................................................')
# Import modules and packages
import numpy  as np
from   numpy import trapz
import pandas as pd
import os
from   datetime import date

# Ask user for directory where scripts are located
directory = input('Enter path for SPICEcore scripts: ')
os.chdir(directory)

# Ask user whether to run phase 1 cleaning (error removal) from this script
choice = input('Run SPICEcore error removal script from this file? Enter Y or N: ')
if choice == 'y' or choice == 'Y':
    
    # Run phase 1 cleaning script
    exec(open('cleaning_cfa.py').read())
    
    # Print header for next part
    print('\n\n...................................................................')
    print('  SPICEcore Dust Data Phase 2 Cleaning: Outliers and Contamination')
    print('...................................................................')
# Ask user for directory where data are located
else: 
    # Print header
    print('\n\n...................................................................')
    print('  SPICEcore Dust Data Phase 2 Cleaning: Outliers and Contamination')
    print('...................................................................')

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
# Make separate copies of the CFA data before and after phase 2 cleaning, for comparison
cfa = cfa_phase1.copy()

# Get the row indices of all measurements within dust events
# These rows will be preserved during subsequent data cleaning
dust_rows = cfa[(cfa['Dust Event?'] == True)].index.values.tolist()

# Get the row indices of all measurements within volcanic events
# These rows will be preserved during subsequent data cleaning
volc_rows = cfa[(cfa['Volcanic Event?'] == True)].index.values.tolist()

# Load file with depth intervals for manual data removal
manual = pd.read_excel('CFA_Manual_Cleaning.xlsx')

#%%
# ---------------------------------------------------------------------------------------
#                            PART 2: Outlier and Contamination Removal
# ---------------------------------------------------------------------------------------

print('\n\n-----------------------------------------------------------------------')
# Get length of dataset from phase 1 cleaning. Use this column to get an accurate count.
length = cfa['Sum 1.1-12'].count()
print('\n\nRemoving outliers and contamination signals.')
print('CFA dataset length after error removal:', length)

# 1) Identify and remove hump-shaped PSD anomalies

# Find humps for all CFA depths
# Inputs: CFA data, minimum depth, maximum depth. Currently using the full core.
humps = find_humps(cfa, 0, 1752)

# Remove all rows in real dust events from the hump list
bad_rows = humps.index.difference(dust_rows)
# Subsequently remove all rows in real volcanic events from the hump list
bad_rows = bad_rows.difference(volc_rows)

# Save all bad data into a separate dataframe
bad_cfa = cfa.loc[bad_rows, :].copy()
bad_cfa['Error Type'] = 'PSD Hump Anomaly'

# NaN values in the bad rows, except depth, age, & boolean columns
cfa.loc[bad_rows, ['Flow Rate', 'ECM', '1', '1.1', '1.2', '1.3', '1.4', '1.5', 
                   '1.6', '1.7', '1.8', '1.9', '2', '2.1', '2.2', '2.3', '2.4', 
                   '2.5', '2.7', '2.9', '3.2', '3.6', '4', '4.5', '5.1', '5.7', 
                   '6.4', '7.2', '8.1', '9', '10', '12', 'CPP', 'Sum 1.1-12']] = np.nan

# Print number of measurements removed
print('\tRows removed: ', len(bad_rows))
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

# NaN values in the bad rows, except depth, age, & boolean columns
cfa.loc[bad_rows, ['Flow Rate', 'ECM', '1', '1.1', '1.2', '1.3', '1.4', '1.5', 
                   '1.6', '1.7', '1.8', '1.9', '2', '2.1', '2.2', '2.3', '2.4', 
                   '2.5', '2.7', '2.9', '3.2', '3.6', '4', '4.5', '5.1', '5.7', 
                   '6.4', '7.2', '8.1', '9', '10', '12', 'CPP', 'Sum 1.1-12']] = np.nan

print('\tRows removed: ', len(bad_rows))

# Update dataset length
length = length - len(bad_rows)

# 3) Identify and remove particle concentration & CPP outliers, using 2-pt. integrals

# UPDATE THIS FUNCTION WITH AARON'S NEWEST VERSION

#outlier_rows, bad_rows = remove_outliers_integrals(cfa, 2, dust_rows, volc_rows)

# Add bad data to the bad CFA dataframe
#bad_cfa = bad_cfa.append(cfa.loc[bad_rows, :], sort = False)
# Label error type
#bad_cfa['Error Type'].fillna('Integral Outlier', inplace = True)

# NaN values in the bad rows, except depth, age, & boolean columns
#cfa.loc[bad_rows, ['Flow Rate', 'ECM', '1', '1.1', '1.2', '1.3', '1.4', '1.5', 
#                   '1.6', '1.7', '1.8', '1.9', '2', '2.1', '2.2', '2.3', '2.4', 
#                   '2.5', '2.7', '2.9', '3.2', '3.6', '4', '4.5', '5.1', '5.7', 
#                   '6.4', '7.2', '8.1', '9', '10', '12', 'CPP', 'Sum 1.1-12']] = np.nan

#print('\tRows removed:', len(bad_rows))

# Update dataset length
#length = length - len(bad_rows)

# 4) Remove remaining manually-identified issues
print('\nManually removing remaining issues.')
# Create empty dataframe to collect intervals to remove
remove_manually = pd.DataFrame()

# Loop through each depth interval in the manual removal file
for start, end in zip(manual['Depth Start (m)'], manual['Depth End (m)']):
    # Subset the CFA data for each depth interval
    selection = select_cfa(cfa, start, end, 'Depth (m)')
    # Append subsetted data to dataframe of data to remove manually
    remove_manually = remove_manually.append(selection, sort = False)

# Drop all rows where everything but depth has already been NaN'd
bad_rows = remove_manually.loc[:, 'Flow Rate'].dropna()
# Get indices of remaining rows
bad_rows = list(bad_rows.index.values)
# Add bad data to the bad CFA dataframe
bad_cfa = bad_cfa.append(cfa.loc[bad_rows, :], sort = False)
# Label error type
bad_cfa['Error Type'].fillna('Manual Removal', inplace = True)

# NaN values in the bad rows, except depth, age, & boolean columns
cfa.loc[bad_rows, ['Flow Rate', 'ECM', '1', '1.1', '1.2', '1.3', '1.4', '1.5', 
                   '1.6', '1.7', '1.8', '1.9', '2', '2.1', '2.2', '2.3', '2.4', 
                   '2.5', '2.7', '2.9', '3.2', '3.6', '4', '4.5', '5.1', '5.7', 
                   '6.4', '7.2', '8.1', '9', '10', '12', 'CPP', 'Sum 1.1-12']] = np.nan

print('\tRows removed: ', len(bad_rows))
# Update dataset length
length = length - len(bad_rows)

# 5) Compute summary statistics before and after phase 2 cleaning, if requested

choice = input('Print summary statistics? Enter Y or N: ')
if choice == 'Y' or choice == 'y':

    print('\n--Phase 1 Cleaning Results--')
    # Input the before & after CFA data into the summary statistics function
    summary_statistics(cfa_phase1)
    print('\n--Phase 2 Cleaning Results--')
    summary_statistics(cfa)
    
#
# ADD DATA INTERPOLATIO HERE, AFTER SUMMARY STATISTICS?
#
    
# 6) Export CFA file to CSV. Report final length.
print('\n\nFinished CFA outlier & contamination removal')
print('\n\tFinal CFA dataset length:', length)

cfa.to_csv('Cleaned_CFA_Phase2_' + str(date.today()) + '.csv')
bad_cfa.to_csv('Bad_CFA_Phase2_' + str(date.today()) + '.csv')
print('\n\tData exported to CSV. Bad data saved in separate file.')
print('-----------------------------------------------------------------------')





