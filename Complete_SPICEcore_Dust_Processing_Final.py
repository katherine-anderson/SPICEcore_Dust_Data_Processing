# --------------------------------------------------------------------------------------
#                     Master SPICEcore Dust Processing Script
#
#    - Dust processing has 2 phases
#      1) Mechanical error removal 
#      2) Removal of outliers & anomalies
#
#    - Gives user the option to run all of the data processing from this script
#      - Select 'Y' on 1st prompt to process all data from beginning to end (Phase 1 and Phase 2)
#      - Select 'N' if you ran the Phase 1 script already and have a CLEANED_CFA_Phase1_ dataset
#
# Phase 2 Dust Processing
#    - Cleans anomalies and outliers from the continuous flow analysis (CFA) data after Phase 1 processing
#      - Preserves data during known dust and volcanic events
#      - Saves 'bad' data into another file, labelled by error type
#      - NaNs 'bad' data in the CFA file and prints error counts
#      - Error types:
#        1) Median absolute deviation (MAD) outliers
#        2) Manually-identified issues which remain
#    - Prints summary statistics
#    - Saves cleaned and 'bad' data to two separate files
#
# Aaron Chesler and Katie Anderson, 7/16/20
# ---------------------------------------------------------------------------------------
#%%
# ---------------------------------------------------------------------------------------
#                               1: FILE PREPARATION
# ---------------------------------------------------------------------------------------
print('\n\n...................................................................')
print('                 SPICEcore DUST DATA PROCESSING')
print('...................................................................')
# Import modules and packages
import numpy  as np
import pandas as pd
import os
from   datetime import date

# Ask user whether to run phase 1 processing (melter error removal)
choice = input('Select from the following options: \n1) All data processing (Phase 1 & Phase 2)\n2) Phase 2 processing only\n\nChoice: ')
if choice == '1':
    
    # Ask user for directory where scripts are located
    directory = input('Enter path for SPICEcore dust scripts: ')
    os.chdir(directory)
    
    # Run Phase 1 processing script
    exec(open('SPICEcore_Dust_Phase1_Processing_Final.py').read())
    
# Ask user for directory where data are located
elif choice == '2': 
    
    # Ask user for directory where scripts are located
    directory = input('Enter path for SPICEcore dust scripts: ')
    os.chdir(directory)
    
    # Run function definitions script
    exec(open('SPICEcore_Dust_Processing_Functions_Final.py').read())
    
    # Get directory for data files
    directory = input('Enter path for SPICEcore dust data: ')
    os.chdir(directory)
else:
    
    print('Invalid choice.')
    # Stop running the program
    exit
    
# Print header for Phase 2 data processing
print('\n\n...................................................................')
print('  SPICEcore Dust Data Phase 2 Cleaning: Outliers and Contamination')
print('...................................................................')

# Load complete CFA file after Phase 2 processing
# Ask user for CFA file to use
file = input('Enter name of the SPICEcore dust file after Phase 1 processing with .csv extension: ')
cfa_phase1 = pd.read_csv(file, header = 0)
del cfa_phase1['Unnamed: 0']
# Make separate copies of the CFA data before and after phase 2 cleaning to compare summary statistics
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
print('\n\nRemoving outliers.')
print('CFA dataset length after error removal:', length)

# 1) Identify and remove particle concentration & CPP outliers, using MAD

# Set # of measurements to use for background medians
window = 500
# Set threshold for accepted Median Absolute Deviations (MAD) (e.g., 2 * MAD)
threshold = 2

# Remove overlapping concentration & CPP outliers
# Inputs: CFA data, dust event indices, volcanic event indices, background window size, and MAD threshold
bad_rows = remove_outliers_MAD(cfa, dust_rows, volc_rows, window, threshold)

# Create empty dataframe to hold bad data
bad_cfa = pd.DataFrame()
# bad_cfa['Error Type'] = 'MAD Outlier'

# Add bad data to the bad CFA dataframe
bad_cfa = bad_cfa.append(cfa.loc[bad_rows, :], sort = False)
# Label error type
bad_cfa['Error Type'] = 'MAD Outlier'

# NaN values in the bad rows, except depth, age, & boolean columns
cfa.loc[bad_rows, ['Flow Rate', 'ECM', '1', '1.1', '1.2', '1.3', '1.4', '1.5', 
                   '1.6', '1.7', '1.8', '1.9', '2', '2.1', '2.2', '2.3', '2.4', 
                   '2.5', '2.7', '2.9', '3.2', '3.6', '4', '4.5', '5.1', '5.7', 
                   '6.4', '7.2', '8.1', '9', '10', '12', 'CPP', 'Sum 1.1-12']] = np.nan

print('\tRows removed: ', len(bad_rows))

# Update dataset length
length = length - len(bad_rows)

# 3) Remove remaining manually-identified issues
print('\n Removing manually-identified issues.')
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

#%%
# 4) Compute summary statistics before and after Phase 2 processing, if requested
choice = input('Print summary statistics? Enter Y or N: ')
if choice == 'Y' or choice == 'y':

    print('\n--Results After Phase 1 Processing--')
    # Input the before & after CFA data into the summary statistics function
    summary_statistics(cfa_phase1)
    print('\n--Results After Phase 2 Processing--')
    summary_statistics(cfa)
        
# 5) Export CFA file to CSV. Report final length.
print('\n\nFinished SPICEcore dust data processing.')
print('\n\tFinal dataset length:', length)

cfa.to_csv('Cleaned_CFA_Phase2_' + str(date.today()) + '.csv')
bad_cfa.to_csv('Bad_CFA_Phase2_' + str(date.today()) + '.csv')
print('\n\tData exported to CSV [Cleaned_CFA_Phase2_...].\n\tBad data saved in separate file [Bad_CFA_Phase2_...].')
print('-----------------------------------------------------------------------')
