#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Script with function definitions. This way, the commands for each calculation can just be written once.
#
# Katie Anderson, 5/28/19
#
# List of functions:
#
# 1) label_core_breaks:   Get a list of indices for each CFA row near a core break (by depth)
# 2) label_volc_events:   Get a list of indices for each row in a volcanic window (by age)
# 3) find_cpp:            Calculate CPP for a CFA dataframe
# 4) label_dust_events:   Get a list of indices for each row in a dust event (by depth)
# 5) find_humps:          Find hump-shaped PSD anomalies in a CFA dataframe
# 6) test_scenario:       Create an outlier removal scenario for a CFA dataframe
# 7) plot_single_psd:     Create histogram of particle counts around a given depth


# In[1]:


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


# In[ ]:


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
        new_cfa = cfa_data[(cfa_data['Age b 1950'] <= (start_year + start_buffer)) & 
                           (cfa_data['Age b 1950'] >= (start_year - end_buffer  ))]
        # Check that there are CFA measurements in the interval around the volcanic events
        if new_cfa.empty: continue  
        else: 
            # Add all indices within volcanic events to the list
            volc_true.extend(new_cfa.index.values.tolist())
            # Add the first index of the CFA data for one volcanic event to the list
            new_volc.append(new_cfa.index[0])
            
    # Return list of rows within buffer dates of volcanic events
    return volc_true, new_volc


# In[1]:


# Function to calculate CPP per measurement. Asks the user to keep/exclude smallest and largest bins.
# Input: CFA data
# Output: List of CPP for each row

def find_cpp(cfa_data):
    # Create dataframe to record particle sums. Needed to prevent dividing by 0
    cpp_df = pd.DataFrame(columns = ['Sum_All', 'Sum_Coarse'])
    
    # Ask the user to include/exclude smallest & largest bins
    choice = input('-->Use smallest and largest bins for CPP? Enter Y or N: ')

    # Create column lists for either option
    if choice == 'y' or choice == 'Y':
        col_list = ['1', '1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', 
                    '1.9', '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', 
                    '3.2', '3.6', '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', 
                    '9', '10', '12']
    if choice == 'n' or choice == 'N':
           col_list = ['1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', 
                    '1.9', '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', 
                    '3.2', '3.6', '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', 
                    '9', '10']
    # Sum particle counts for each measurement using the above columns
    cpp_df['Sum_All'] = cfa_data[col_list].sum(axis = 1)
    # Check for negative counts
    if min(cpp_df['Sum_All']) < 0: 
        print('CPP function found negative sum of all particles.')

    # Remake the column lists for only the coarse particles (>= 4.5 um)
    if choice == 'y' or choice == 'Y':
        col_list = ['4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12']

    if choice == 'n' or choice == 'N':
        col_list = ['4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10']

    # Sum coarse particle counts for each measurement using the above columns
    cpp_df['Sum_Coarse'] = cfa_data[col_list].sum(axis = 1)
    # Check for negative counts
    if min(cpp_df['Sum_Coarse']) < 0: 
        print('CPP function found negative sum of coarse particles.')
    
    # Remove rows with 0 sum_all counts. Don't divide by 0
    cpp_df = cpp_df[cpp_df['Sum_All'] > 0]

    # Return a series of the percent of particles that are coarse per row
    return(cpp_df['Sum_Coarse'] / cpp_df['Sum_All'] * 100)


# In[ ]:


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


# In[2]:


# Identifies CFA measurements with hump PSD anomalies
# Inputs: CFA data and the min/max depths in which to find the humps
# Output: CFA dataframe with only the hump measurements
# Hump criteria: Measurements where all bins 3.5-10 um have higher counts than the
#    average count for bins 1.5-2.9 um.

def find_humps(cfa_data, min_depth, max_depth):
    
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
    print('Number of measurements:    ', len(cfa_data))
    print('Number of discrete events: ', len(new_humps))
    
    return cfa_data


# In[ ]:


# Function to remove outliers given different background & sensitivity conditions
# Inputs: CFA data, list of dust event rows, # of measurements for background levels, # of stdevs. accepted
# Output: Outlier counts, dataframe with overlapping outliers removed

def test_threshold(cfa_data, dust_indices, background_interval, threshold_stdev):
    print('\nNEW OUTLIER REMOVAL SCENARIO')
    
    # Calculate rolling standard deviations
    # Will calculate std. if only 1 measurement in the window that isn't NaN
    cpp_background  = cfa_data['CPP'].rolling(background_interval, min_periods = 1).std()
    conc_background = cfa_data['Sum 1.1-10'].rolling(background_interval, min_periods = 1).std()   

    # Subsetting CFA data for the outliers
    conc_peaks =  cfa_data[(cfa_data['CPP'] >= ( threshold_stdev *  conc_background))]
    cpp_peaks  =  cfa_data[(cfa_data['CPP'] >= ( threshold_stdev *  cpp_background)) ]

    # Want to find when these outliers occur at the same time
    overlap = conc_peaks.index.intersection(cpp_peaks.index)
    # Prevent rows in real dust events from being removed
    overlap = overlap.difference(dust_indices)
    
    # Ask the user whether or not to display # of outliers found & removed
    choice1 = input('\n-->Print outlier counts? Enter Y or N: ')
    if choice1 == 'y' or choice1 == 'Y':    
        print('\nConc. outliers:              ', len(conc_peaks))
        print('CPP outliers:                ', len(cpp_peaks))
        print('Overlap:                     ', len(overlap))
    
    if 'Volcanic Event?' in cfa_data.columns:
        # Subset the outlier dataframes for rows within volcanic events
        volc_conc_peaks = conc_peaks[(conc_peaks['Volcanic Event?'] == True)]
        volc_cpp_peaks  = cpp_peaks [(cpp_peaks ['Volcanic Event?'] == True)]

        # Count the number of overlapped outliers which occur within a volcanic event
        volc_overlap = volc_conc_peaks.index.intersection(volc_cpp_peaks.index)
        
        if choice1 == 'y' or choice1 == 'Y':
            print('\nConc. outliers w/in volcanic event:', len(volc_conc_peaks))
            print('CPP outliers w/in volcanic event:  ', len(volc_cpp_peaks))
            print('Overlap w/in volcanic event:       ', len(volc_overlap))
    
        # Want to have the option of preserving volcanic events or not
        choice2 = input('\n-->Remove outliers at volcanic events? Enter Y or N: ')
        if choice2 == 'y' or choice2 == 'Y':
            # Remove variable is the indices at which to NaN values
            remove = overlap
            if choice1 == 'y' or choice1 == 'Y':
                print('\nSUMMARY OF OUTLIER REMOVAL')
                print('1) Rows Removed: ', len(remove))
        if choice2 == 'n' or choice2 == 'N':
            # Subtract the volc_overlap indices from the overlap indices
            remove = overlap.difference(volc_overlap)
            if choice1 == 'y' or choice1 == 'Y':
                print('\nSUMMARY OF OUTLIER REMOVAL')
                print('1) Rows Removed: ', len(remove))
    
        # Count number of discrete volcanic events in the rows we're about to remove
        # Get a copy of the CFA data with only the 'New Volcanic Event?' column
        temp_cfa = cfa_data.loc[remove, :].copy()
        volc_event = temp_cfa[(temp_cfa['New Volcanic Event?'] == True)]
        print('2) Number of discrete volcanic events in removed rows:', len(volc_event))
        
        # Count number of core breaks in the rows we're about to remove
        temp_cfa = cfa_data.loc[remove, :].copy()
        break_outliers = temp_cfa[(temp_cfa['New Break?'] == True)]
        print('3) Number of discrete core breaks in removed rows:    ', len(break_outliers))
    else: 
        # Just remove the rows w/ overlapping CPP and conc. outliers
        remove = overlap
        
        # Count number of core breaks in the rows we're about to remove
        temp_cfa = cfa_data.loc[remove, :].copy()
        break_outliers = temp_cfa[(temp_cfa['New Break?'] == True)]
        
        print('\nSUMMARY OF OUTLIER REMOVAL')
        print('1) Rows Removed:', len(remove))
        print('2) Number of discrete core breaks in removed rows:', len(break_outliers))
        
    # Change all rows with overlapping outliers to NaN
    # Don't NaN the 'Break?', 'New Break?', 'Volcanic Event?', and 'New Volcanic Event?' columns
    cfa_data.loc[remove, ['Depth (m)', 'Age b 1950', 'Flow Rate', 'ECM', '1', '1.1',
                            '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9', '2', 
                            '2.1', '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', '3.2', '3.6', 
                            '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12',
                            'Original Depth (m)', 'CPP', 'Sum 1.1-10']] = np.nan

    return cfa_data


# In[ ]:


# Function to create 1 bar plot of particle counts within x-centimenters of a point
#     Inputs: CFA dataframe, depth of interest, depth range around that point
#     Output: Bar plot of particle counts per bin, with option to save

def plot_single_psd(cfa_data, point, depth_range):

    # The range around the main point of interest
    point_min = point - depth_range
    point_max = point + depth_range

    # Subset the CFA data to range around given point
    point_cfa = cfa_data[(cfa_data['Depth (m)'] >= point_min) 
                        & (cfa_data['Depth (m)'] <= point_max)]
    
    # Check for and remove columns
    col_list = ['1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', 
                    '1.9', '2', '2.1', '2.2', '2.3', '2.4', '2.5', '2.7', '2.9', 
                    '3.2', '3.6', '4', '4.5', '5.1', '5.7', '6.4', '7.2', '8.1', 
                    '9', '10']
    point_cfa = point_cfa[col_list]
    
    # Sum particles by column around main & comp. points
    point_count = point_cfa.sum(axis = 0)
    
    # Make figures
    fig, ax = plt.subplots(figsize = (5,5))
    
    ax.bar(col_list, point_count, width = 1, color = 'black');
    ax.set_xticks([0,10,20,29])
    ax.tick_params(labelsize = 14)
    ax.set_ylabel('Counts (#/uL)', fontsize = 16)
    ax.set_title(str(point) + ' Meters +/- ' + str(depth_range) + ' Meters', fontsize = 18)    
    
    choice = input('Save Figure? Enter Y or N: ')
    if choice == 'y' or choice == 'Y':
        os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\SPICE\\PSD Plots')
        plt.savefig('PSD_1_' + str(point) + '.png')

