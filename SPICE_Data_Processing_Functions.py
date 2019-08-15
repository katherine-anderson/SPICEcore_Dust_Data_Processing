# In[1]:


# Script with function definitions. This way, the commands for each calculation can just be written once.
#
# Katie Anderson, 8/13/19
#
# List of functions:
#
#  1) correct_meltday:           Corrects low Abakus values during melt day 7/19/2016
#  2) label_core_breaks:         Get a list of indices for each CFA row near a core break (by depth)
#  3) label_volc_events:         Get a list of indices for each row in a volcanic window (by age)
#  4) label_dust_events:         Get a list of indices for each row in a dust event (by depth)
#  5) find_cpp:                  Calculate CPP for a CFA dataframe
#  6) find_humps:                Find hump-shaped PSD anomalies in a CFA dataframe
#  7) median_absolute_deviation: Calculate MAD for one column of CFA data
#  8) remove_outliers_MAD:       Remove outliers from the CFA data, using MAD
#  9) remove_outliers_integrals: Get a list of indices with an integral outlier
# 10) summary_statistics:        Print summary statistics for dust concentration & CPP during data cleaning

#%%
# 1)
# Function to correct anomalously low Abakus values during melt day 7/19/2019
# Inputs: CFA dataframe
# Output: Corrected CFA dataframe

def correct_meltday(cfa_data):
    # Correcting low Abakus values during melt day 7/19/2016 (302-312 m)
    # Get median of the days before (290-302 m) and after (312-322)
    # Correct the 7/19 data so it has the same median as the neighboring days
    
    
    # Select data from the day before (290-302 m)
    cfa_BandA = cfa.loc[cfa['Depth (m)'] >= 290].copy()
    cfa_BandA = cfa_BandA.loc[cfa_BandA['Depth (m)'] <  302].copy()
    
    # Select data from the day after (312-322 m)
    cfa_A = cfa.loc[cfa['Depth (m)'] >= 312].copy()
    cfa_A = cfa_A.loc[cfa_A['Depth (m)'] <  322].copy()
    
    # Combine the data from the day before and the day after
    cfa_BandA = cfa_BandA.append(cfa_A)
    # Get median values in each Abakus column for the combined data
    med_BandA = cfa_BandA.median(axis = 0)
    med_BandA = med_BandA.drop(labels = ['Depth (m)', 'Flow Rate', 'ECM'])
    
    # Select data from Jul 19 (302-312 m), which need the correction
    cfa_Jul19 = cfa.loc[cfa['Depth (m)'] >= 302].copy()
    cfa_Jul19 = cfa_Jul19.loc[cfa_Jul19['Depth (m)'] <  312].copy()
    # Get median values in each Abakus column for the day to correct
    med_Jul19 = cfa_Jul19.median(axis = 0)
    med_Jul19 = med_Jul19.drop(labels = ['Depth (m)', 'Flow Rate', 'ECM'])
    
    # Divide the before/after medians by the Jul 19 medians
    # This is the correction to apply to Jul 19 to get it to the right levels
    # Only want positive, not-NaN correction values
    correction = med_BandA.div(med_Jul19)
    # Drop NaNs, where it divided by zero
    correction = correction.dropna()
    # Drop 0s-- don't want to multiply everything by 0
    correction = correction[(correction > 0)]
    
    # Get the list of the remaining columns to correct
    columns = correction.index.values
    # Select these columns in the Jul 19 data
    cfa_Jul19 = cfa_Jul19.loc[:, columns]
    # Multiply each value from Jul 19 by the correction for its column
    new_Jul19 = cfa_Jul19.mul(correction, axis = 1)
    
    # Update the corrected columns in the Jul 19 data
    cfa_Jul19.loc[:, columns] = new_Jul19
    # New medians now equal the before/after median
    
    # Update original CFA data with the new values in the corrected dataframe
    cfa.update(cfa_Jul19)
    # Return corrected CFA dataframe
    return(cfa)

#%%
# 2)
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
# 3)
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
# 4)
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
# 5) 
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
    
    # Select the non-NaN CFA rows
    depth_isnull = cfa_data['Depth (m)'].isnull()
    valid_depth = depth_isnull[depth_isnull == False]
    # Get the indices of these rows
    valid_depth = list(valid_depth.index.values)
    # Select these rows in the CFA data
    cfa_data = cfa_data.loc[valid_depth].copy() 
    
    # Sum particle counts for each measurement using the above columns
    cpp_df['Sum_All'] = cfa_data[col_list].sum(axis = 1)

    # Remake the column lists for only the coarse particles (>= 4.5 um)
    col_list = ['4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12']

    # Sum coarse particle counts for each measurement using the above columns
    cpp_df['Sum_Coarse'] = cfa_data[col_list].sum(axis = 1)
    
    # Remove rows with 0 sum_all counts. Don't divide by 0
    cpp_df = cpp_df[cpp_df['Sum_All'] > 0]

    # Return a series of the percent of particles that are coarse per row
    return(cpp_df['Sum_Coarse'] / cpp_df['Sum_All'] * 100)

#%%
# 6) 
# Identifies CFA measurements with hump PSD anomalies
# Inputs: CFA data and the min/max depths in which to find the humps
# Output: CFA dataframe with only the hump measurements
# Hump criteria: Measurements where all bins 3.5-10 um have higher counts than the
#    average count for bins 1.5-2.9 um.

def find_humps(cfa_data, min_depth, max_depth):
    print('\nRemoving PSD hump anomalies')

    # Subset the CFA data for the selected depth range
    cfa_data = cfa_data[(cfa_data['Depth (m)'] >= min_depth) & 
                        (cfa_data['Depth (m)'] <= max_depth)]
    
    # Make a copy of CFA data for bins 3.2-10
    humps_col_list = ['3.2', '3.6', '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10']
    
    cfa_humps = cfa_data[humps_col_list]
    
    # Make copy of CFA data for bins 1.5-2.9
    smalls_col_list = ['1.5', '1.6', '1.7', '1.8', '1.9', '2', '2.1', '2.2',
                       '2.3', '2.4', '2.5', '2.7', '2.9']
    
    cfa_smalls = cfa_data[smalls_col_list]  
    
    # Get mean concentration for bins 1.5-2.9 per row
    smalls_mean = cfa_smalls.mean(axis = 'columns')
    
    # Subtract the 1.5-2.9 bin mean from the 3.2-10 values
    cfa_humps = cfa_humps.subtract(smalls_mean, axis = 'index')
    
    # If all subtracted values are positive, it is a hump
    # Mark all positive differences as True
    criteria = cfa_humps > 0
    # Check for rows with only True values
    criteria = criteria.all(axis = 'columns')

    # Subset the CFA data for rows where all values are elevated above small particles
    cfa_data = cfa_data[criteria]
    
    # Count the number of discrete humps
    # Copy the first two columns of the CFA data into new dataframe
    humps_diff = cfa_data.loc[:, 'Depth (m)':'ECM'].copy()
    # Subtract each row from the one before (to get diff in depth)
    humps_diff = humps_diff.diff()
    # Subset the diff dataframe for rows where depth changes by 3+ cm
    # ~3 cm melt resolution. >3 cm diff = new hump event
    new_humps = humps_diff[(humps_diff['Depth (m)'] >= 0.03)]
    
    # Report the number of hump measurements and events
    choice = input('\tPrint detailed PSD hump anomaly counts? Enter Y or N: ')
    if choice == 'y' or choice == 'Y':
        print('\tNumber of measurements:    ', len(cfa_data))
        print('\tNumber of discrete events: ', len(new_humps))
    
    return cfa_data
#%%
# 7)  
# Function to calculate Median Absolute Deviation (MAD)
# Inputs: One-dimensional dataset (like CFA particle concentration or CPP)
# Output: MAD

def median_absolute_deviation(x, axis = None):
    # Get the absolute value of every element minus the overall mean
    deviation = abs(x - x.median())
    
    # Return the median of that deviation
    return deviation.median()
#%%
# 8) 
# Function to remove outliers given different background & sensitivity conditions
# Inputs: CFA data, list of dust event rows, list of volcanic event rows, background window size, MAD threshold
# Outputs: Dataframe with overlapping outliers removed, list of outlier indices

def remove_outliers_MAD(cfa_data, dust_indices, volc_indices, background_interval, threshold):
    print('\nRemoving MAD outliers')
    
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
    # Ask the user whether or not to display # of outliers found & removed
    choice2 = input('\tPrint detailed MAD outlier counts? Enter Y or N: ')
    
    if choice1 == 'n' or choice1 == 'N':
        # Remove variable has the indices at which to NaN values
        remove = overlap
    elif choice1 == 'y' or choice1 == 'Y':
        # Subtract the volcanic event indices from the overlapping outlier indices
        remove = overlap.difference(volc_indices)
    else:
        print('Invalid entry. Defaulted to preserving outliers at volcanic events.')
        remove = overlap.difference(volc_indices)
    
    # Count number of discrete volcanic events in the rows we're about to remove
    # Get a copy of the CFA data with only the 'New Volcanic Event?' column
    temp_cfa = cfa_data.loc[remove, :].copy()
    volc_event = temp_cfa[(temp_cfa['New Volcanic Event?'] == True)]
    if choice2 == 'y' or choice2 == 'Y':
        print('\tNumber of discrete volcanic events in removed rows:', len(volc_event))
        
    # Count number of core breaks in the rows we're about to remove
    temp_cfa = cfa_data.loc[remove, :].copy()
    break_outliers = temp_cfa[(temp_cfa['New Break?'] == True)]
    if choice2 == 'y' or choice2 == 'Y':
        print('\tNumber of discrete core breaks in removed rows:    ', len(break_outliers))
        
    return remove
#%%
# 9) 
# Function to remove outliers using rolling 2-pt integrals (Aaron)
# Inputs: CFA data, stdev threshold above median, list of dust event rows, list of volcanic rows
# Output: List of rows to remove from the CFA data

def remove_outliers_integrals(cfa_data, threshold, dust_indices, volc_indices):
    print('\nRemoving integral outliers')
    # Calculate thresholds for concentration & CPP outliers
    
    # Lists for integral values
    # y is the integral size
    y = 2
    conc_tpz = []
    cpp_tpz  = []
    
    for x in range(0, len(cfa_data), 2):
        # Concentration integrals
        # Make sure neither value is NaN. All NaN values in concentration column should be NaN in CPP column
        if cfa_data.loc[x, 'Sum 1.1-12'] == np.nan or cfa_data.loc[y, 'Sum 1.1-12'] == np.nan: 
            continue
        else:
            conc_tpz.append(trapz(cfa_data['Sum 1.1-12'][x: (x+y)]))
            # CPP integrals
            cpp_tpz.append(trapz(cfa_data['CPP'][x: (x+y)]))
               
    conc_stdev      = np.nanstd(conc_tpz)
    conc_median     = np.nanmedian(conc_tpz)
    # Threshold for error
    conc_error = (threshold * conc_stdev) + conc_median

    cpp_stdev      = np.nanstd(cpp_tpz)
    cpp_median     = np.nanmedian(cpp_tpz)
    # Threshold for error
    cpp_error  = (threshold * cpp_stdev) + cpp_median
                        
    # Concentration outliers
    y = 2
    conc_outlier_indices = []
    cpp_outlier_indices  = []
    
    for x in range(0,len(cfa_data), 2):
    
        # Concentration outliers
        if trapz(cfa_data['Sum 1.1-12'].iloc[x:y]) >= conc_error:  # If area increases outside threshold
            # Only append indices that aren't already in the list. Avoids duplicates.
            if x not in conc_outlier_indices:
                conc_outlier_indices.append(x)
            if y not in conc_outlier_indices:
                conc_outlier_indices.append(y)
            
        # CPP outliers        
        if trapz(cfa_data['CPP'].iloc[x:y]) >= cpp_error:  # If area increases outside threshold
            if x not in cpp_outlier_indices:
                cpp_outlier_indices.append(x)
            if y not in cpp_outlier_indices:
                cpp_outlier_indices.append(y)
        y = y + 2

    # Convert index lists to sets to find the overlap between them
    conc_outlier_set = set(conc_outlier_indices)
    cpp_outlier_set  = set(cpp_outlier_indices)

    # Find overlapping indices with both concentration & CPP outliers
    overlap = conc_outlier_set.intersection(cpp_outlier_set)
    # Convert to list again and sort
    overlap = list(overlap)
    overlap.sort()
    
    # Select those rows in the CFA data
    cfa_data = cfa_data.loc[overlap, :]
    # Prevent real dust & volc. events from being removed
    remove = cfa_data.index.difference(dust_indices)
    remove = cfa_data.index.difference(volc_indices)

    # Return list of outlier rows and the rows to remove
    return overlap, remove
  #%%  
# 10)    
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
        print('Dust number concentration (#/mL):')
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
        
    else: print('Input data with particle sum and CPP columns')

