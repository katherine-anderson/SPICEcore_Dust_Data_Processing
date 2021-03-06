# --------------------------------------------------------------------------------------
#                     SPICEcore DUST PROCESSING FUNCTIONS

# Supplementary script with function definitions for functions called in other scripts
# No need to run this script on its own- the other scripts run it internally
#
# List of functions:
#
#  1) correct_meltday:           Correct time units during melt day 7/19/2016
#  2) label_core_breaks:         Get a list of indices for each continuous flow analysis (CFA) row near a core break
#  3) label_volc_events:         Get a list of indices for each row in a volcanic window (by age)
#  4) label_dust_events:         Get a list of indices for each row in a dust event (by depth)
#  5) find_cpp:                  Calculate CPP for a CFA dataframe
#  6) median_absolute_deviation: Calculate median absolute deviation (MAD) for one column of CFA data
#  7) remove_outliers_MAD:       Remove outliers from the CFA data, using MAD
#  8) select_cfa:                Subset CFA data for given depth or age range
#  9) summary_statistics:        Print summary statistics for dust concentration & CPP during data cleaning
    
# Katie Anderson, 7/16/20
# ---------------------------------------------------------------------------------------
#%%
# Function to correct time units during melt day 7/19/2019
# Inputs: CFA dataframe
# Output: Corrected CFA dataframe

def correct_meltday(cfa_data):
    # Correcting low Abakus values during melt day 7/19/2016 (302-312 m)
    # Multiply all Abakus values by 60. Data were recorded in /min, not /sec
       
    # Select data from Jul 19 (302-312 m), which need the correction
    cfa_Jul19 = cfa_data.loc[cfa_data['Depth (m)'] >= 302].copy()
    cfa_Jul19 = cfa_Jul19.loc[cfa_Jul19['Depth (m)'] <  312].copy()
    
    # Select only the Abakus columns
    cfa_Jul19 = cfa_Jul19.loc[:, '1':'12']
    
    # Multiply all Abakus values by 60 to convert to /sec
    cfa_Jul19 = cfa_Jul19.mul(60)
    
    # Update original CFA data with the new values in the corrected dataframe
    cfa.update(cfa_Jul19)
    
    # Return corrected CFA dataframe
    return(cfa)

#%%
# Function to get indices of all CFA measurements taken within a specified core break range
# Inputs: CFA dataframe, core breaks dataframe, specified +/- core break range (in meters)
# Output: List of rows within core breaks

def label_core_breaks(cfa_data, core_breaks, core_range):
    
    # Create an empty list to record the CFA measurements which occurred around a core break range
    break_true  = []
    # Create an empty list to record the first row within each core break range
    new_break = []
    
    # Subset the CFA data for depths within range of each core break
    for corebreak in core_breaks['Depth (m)']:
        new_cfa = cfa_data[(cfa_data['Depth (m)'] >= (corebreak - core_range)) & 
                           (cfa_data['Depth (m)'] <= (corebreak + core_range))]
        # Check that there are CFA measurements in the core break interval
        if new_cfa.empty: continue  
        # Add all depth values in the core break interval to a list
        else: 
            # Add all indices within core breaks to the list
            break_true.extend(new_cfa.index.values.tolist())
            # Add the first index of the CFA data for one core break to the list
            new_break.append(new_cfa.index[0])
            
    # Return list of indices occurring within core breaks
    return break_true, new_break
#%%
# Function to get indices of all CFA measurements taken within range of years around volcanic events
# Inputs: Holocene CFA, Holocene volcanic dates, before/after buffers, in years
# Output: List of rows within volcanic range

def label_volc_events(cfa_data, volc_record, start_buffer, end_buffer):
    
    # Create an empty list to record the CFA measurements which occurred near a volcanic event
    volc_true  = []
    # Create an empty list to record the first row within each volcanic event
    new_volc = []
        
    # Subset the CFA data for depths within range of other volcanic events
    for start_year in volc_record['Start Year (b1950)']:
        new_cfa = cfa_data[(cfa_data['AgeBP'] <= (start_year + start_buffer)) & 
                           (cfa_data['AgeBP'] >= (start_year - end_buffer  ))]
        # Check that there are CFA measurements in the interval around the volcanic events
        if new_cfa.empty: continue  
        else: 
            # Add all indices within volcanic events to the list
            volc_true.extend(new_cfa.index.values.tolist())
            # Add the first index of the CFA data for one volcanic event to the list
            new_volc.append(new_cfa.index[0])
            
    # Return list of rows within buffer dates of volcanic events
    return volc_true, new_volc
#%%
# Function to get a list of rows within dust events
# Inputs: CFA data, dust event dataframe with depth intervals
# Output: List of rows within dust events

def label_dust_events(cfa_data, dust_depths):
    
    # Create an empty list to record the CFA measurements within depth range of dust events
    dust_event_true = []
    
    # Subset the CFA data for depths within range of dust events
    for index, row in dust_events.iterrows():
        new_cfa = cfa_data[(cfa_data['Depth (m)'] >= row['Dust Event Start (m)']) &
                           (cfa_data['Depth (m)'] <= row['Dust Event End (m)'])]
        # Check that there are CFA measurements in the dust event depth intervals
        if new_cfa.empty: continue
        # Add all indices within dust events to the list
        else:
            dust_event_true.extend(new_cfa.index.values.tolist())
            
    # Return list of rows within dust events 
    return dust_event_true
#%%
# Function to calculate CPP per measurement
# Input: CFA data
# Output: List of CPP for each row

def find_cpp(cfa_data):
   
    # Create dataframe to record particle sums. Need to prevent dividing by 0
    cpp_df = pd.DataFrame(columns = ['Sum_All', 'Sum_Coarse'])
    
    col_list = ['1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9', 
                '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', '3.2', 
                '3.6', '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', 
                '10', '12']
    
    # Sum particle counts for each measurement using the above columns
    cpp_df['Sum_All'] = cfa_data[col_list].sum(min_count = 1, axis = 1)

    # Remake the column lists for only the coarse particles (>= 4.5 um)
    col_list = ['4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12']

    # Sum coarse particle counts for each measurement using the above columns
    cpp_df['Sum_Coarse'] = cfa_data[col_list].sum(min_count = 1, axis = 1)
    
    # Remove rows with 0 sum_all counts. Don't divide by 0
    cpp_df = cpp_df[cpp_df['Sum_All'] > 0]

    # Return a series of the percent of particles that are coarse per row
    return(cpp_df['Sum_Coarse'] / cpp_df['Sum_All'] * 100)

#%%
# Function to calculate Median Absolute Deviation (MAD)
# Inputs: One-dimensional dataset (like CFA particle concentration or CPP)
# Output: MAD

def median_absolute_deviation(x, axis = None):
    # Get the absolute value of every element minus the overall mean
    deviation = abs(x - x.median())
    
    # Return the median of that deviation
    return deviation.median()
#%%
# Function to remove outliers given different background & sensitivity conditions
# Inputs: CFA data, list of dust event rows, list of volcanic event rows, background window size, MAD threshold
# Outputs: Dataframe with overlapping outliers removed, list of outlier indices

def remove_outliers_MAD(cfa_data, dust_indices, volc_indices, background_interval, threshold):
    print('\nRemoving MAD outliers.')
    
    # Calculate rolling medians and overall median absolute deviation (MAD)
    # Will calculate if 3 measurements in the window that aren't NaN
    
    cpp_background  = cfa_data['CPP'].rolling(background_interval, min_periods = 3).median()
    conc_background = cfa_data['Sum 1.1-12'].rolling(background_interval, min_periods = 3).median()
    
    cpp_mad  = median_absolute_deviation(cfa_data['CPP'])
    conc_mad = median_absolute_deviation(cfa_data['Sum 1.1-12'])

    # Subsetting CFA data for the outliers
    # Point is an outlier if it exceeds threshold * MAD from the background
    cpp_peaks  = cfa_data[(cfa_data['CPP']        >= (cpp_background  + threshold * cpp_mad))]
    conc_peaks = cfa_data[(cfa_data['Sum 1.1-12'] >= (conc_background + threshold * conc_mad))]

    # Want to find when these outliers occur at the same time
    overlap = conc_peaks.index.intersection(cpp_peaks.index)
    # Prevent rows in real dust events from being removed
    overlap = overlap.difference(dust_indices)
    
    # Ask the user whether or not to preserve outliers at volcanic events
    choice1 = input('\tPreserve outliers at volcanic events? Enter Y or N: ')
    
    if choice1 == 'n' or choice1 == 'N':
        # Remove variable has the indices at which to NaN values
        remove = overlap
    elif choice1 == 'y' or choice1 == 'Y':
        # Subtract the volcanic event indices from the overlapping outlier indices
        remove = overlap.difference(volc_indices)
    else:
        print('Invalid entry. Defaulted to preserving outliers at volcanic events.')
        remove = overlap.difference(volc_indices)
       
    return remove
#%%
# Function to subset CFA data for given depth or age range
#     Inputs: CFA dataframe, starting value, ending value, how ('Depth' or 'Age')
#     Outputs: Subset of CFA data

def select_cfa(cfa_data, lower, upper, variable):
    
    if variable == 'Depth (m)':
        cfa_data = cfa_data[cfa_data['Depth (m)'] >= lower]
        cfa_data = cfa_data[cfa_data['Depth (m)'] <  upper]
        return cfa_data
    elif variable == 'AgeBP':
        cfa_data = cfa_data[cfa_data['AgeBP'] >= lower]
        cfa_data = cfa_data[cfa_data['AgeBP'] <  upper]
        return cfa_data
    else:
        print('Invalid entry')
        return cfa_data

#%%  
# Function to print summary statistics for dust concentration & CPP during data cleaning
#     Inputs: CFA dataframe with particle sum and CPP columns
#     Output: None. Prints summary statistics.

def summary_statistics(cfa_data):
    
    if 'Sum 1.1-12' in cfa_data.columns and 'CPP' in cfa_data.columns:  
        # Make local copies of particle concentration & CPP
        row_sums = cfa_data.loc[:, 'Sum 1.1-12'].copy()
        cpp      = cfa_data.loc[:, 'CPP'].copy()
        
        # Convert number concentration to # / mL
        number_conc = row_sums.mul(1000).copy()
        
        # Call function to calculate MAD for the particle concentration & CPP
        conc_mad = median_absolute_deviation(cfa_data['Sum 1.1-12'])
        cpp_mad  = median_absolute_deviation(cfa_data['CPP'])
    
        # Skip NaNs when calculating these
        print('Dust Number Concentration (#/mL):')
        print('    Mean:   %.2f' % (np.nanmean(number_conc)))
        print('    Median: %.2f' % (np.nanmedian(number_conc)))
        print('    Min:    %.2f' % (np.nanmin(number_conc)))
        print('    Max:    %.2f' % (np.nanmax(number_conc)))
        print('    StDev:  %.2f' % (np.nanstd(number_conc)))
        print('    MAD:    %.2f' % conc_mad)
    
        print('\nCoarse Particles (%):')
        print('    Mean:   %.2f' % (np.nanmean(cpp)))
        print('    Median: %.2f' % (np.nanmedian(cpp)))
        print('    Min:    %.2f' % (np.nanmin(cpp)))
        print('    Max:    %.2f' % (np.nanmax(cpp)))
        print('    StDev:  %.2f' % (np.nanstd(cpp)))
        print('    MAD:    %.2f' % cpp_mad)
        
    else: print('Input data with particle sum and CPP columns.')

#%%