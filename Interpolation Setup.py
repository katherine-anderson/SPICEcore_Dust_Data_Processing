# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 14:21:12 2019

@author: katie
"""
#cfa.loc[288710:288730, 'AgeBP':'CPP']
#%%
# Import modules
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# Load function definition script
os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\Git Projects\\SPICE_Data_Processing')
# Need to read this script for the 'select_cfa' function
exec(open('SPICE_Data_Processing_Functions.py').read())

# Change to SPICE data directory
os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\SPICE\\Data')
#%%
# Load cleaned CFA data, delete extra columns
cfa = pd.read_csv('Cleaned_CFA_Phase2_2019-10-22.csv')
del cfa['Unnamed: 0']
#%%
# Function to create interpolation setup
# Inputs: Cleaned CFA data
# Output: Dataframe with time intervals and years/rows to interpolate over

def setup_interp(cfa_data):
    # Set up interpolation dataframe
    interp = pd.DataFrame()
    
    # Need to get lists of incrementing age ranges
    start_age = []
    end_age = []
    # Get youngest age value to start
    years = cfa_data.loc[0, 'AgeBP']
    
    # Increment ages until the max age of the data
    # Column 33 is AgeBP
    while years <= cfa_data.iloc[-1, 33]:
        # Add the start age to the list
        start_age.append(years)
        # Increment by 500 years
        years = years + 500
        # Add the end age to the list
        end_age.append(years)
    
    # Add lists to new interp dataframe columns
    interp['Starting Age'] = start_age
    interp['Ending Age']   = end_age
    # Create empty column to indicate the # of years to interpolate over in that interval
    interp['Years to Interp'] = np.nan
    
    # Subset the data for different time intervals (arbitrary)
    # Set different limits on the # of years to interpolate for each interval
    holocene = interp[interp['Starting Age'] <= 12000].copy()
    holocene['Years to Interp'] = 0.5
    
    lgm = interp[interp['Starting Age'] > 12000].copy()
    lgm = lgm[lgm['Starting Age'] <= 22000]
    lgm['Years to Interp'] = 1
    
    glacial = interp[interp['Starting Age'] > 22000].copy()
    glacial['Years to Interp'] = 2
    
    # Merge all 3 dataframes, now with the 'years to interp' data
    dfs = [holocene, lgm, glacial]
    interp = pd.concat(dfs)

    # Create empty list to store the # of rows to interpolate over for each interval
    interp_rows = []
    
    # Loop through each time interval in the interpolation dataframe
    for start, end, years in zip(interp['Starting Age'], interp['Ending Age'], interp['Years to Interp']):
        # Subset the CFA data for this age interval
        temp = select_cfa(cfa_data, start, end, 'AgeBP')
        # Only keep the age values
        temp = temp.loc[:, 'AgeBP'].copy()
        # Get the difference in age between each row
        age_diff = temp.diff()
        # Get the mean age difference between rows
        mean_diff = np.nanmean(age_diff)
        # Find how many rows are needed to reach the specified number of years
        number_rows = years / mean_diff
        # Round to nearest integer. This is the number of rows to interpolate over for this interval
        number_rows = np.round(number_rows, decimals = 0)
        # Change to integer
        number_rows = number_rows.astype(np.int64)
        # Append the number of rows
        interp_rows.append(number_rows)
    
    # Add list to new column in the interpolation dataframe
    interp['Rows to Interp'] = interp_rows    

    return interp
#%%
# INTERPOLATION

def interpolate_cfa(cfa_data, interp_df):
    # Create empty series for final concentration & CPP columns
    final_conc = pd.Series()
    final_cpp  = pd.Series()
    
    # Ask user how to interpolate
    choice = input('Method:\n1) Linear\n2) Polynomial\n3) Spline\n\nChoice: ')
    if choice == '1': how = 'linear'
    elif choice == '2':
        how = 'polynomial'
        order_choice = input('Enter integer for order: ')
        # Make sure it's an integer
        order_choice = int(order_choice)
    elif choice == '3':
        how = 'spline'
        order_choice = input('Enter integer for order: ')
        # Make sure it's an integer
        order_choice = int(order_choice)
    else:
        print('Invalid choice.')
        return
    
    # If 'linear'
    if choice == '1':
        # Loop through each interpolation interval
        for start, end, rows in zip(interp_df['Starting Age'], interp_df['Ending Age'], interp_df['Rows to Interp']):
            # Subset the CFA data for this age interval
            temp = select_cfa(cfa_data, start, end, 'AgeBP')
            # Interpolate concentration up to x number of rows
            # Limit_direction = backward, so it interpolates forward through time
            # No order specified for linear interpolation
            interp_conc = temp['Sum 1.1-12'].interpolate(method = how, limit = rows, limit_direction = 'backward')
            # Interpolate CPP
            interp_cpp = temp['CPP'].interpolate(method = how, limit = rows, limit_direction = 'backward')    
            # Append interpolated values to final lists before moving to next interval
            final_conc = final_conc.append(interp_conc)
            final_cpp  = final_cpp.append(interp_cpp)
            
        # Return series of interpolated concentration & CPP data
        return final_conc, final_cpp 
        
    # If 'polynomial' or 'spline'        
    elif choice == '2' or choice == '3':
        # Loop through each interpolation interval
        for start, end, rows in zip(interp_df['Starting Age'], interp_df['Ending Age'], interp_df['Rows to Interp']):
            # Subset the CFA data for this age interval
            temp = select_cfa(cfa_data, start, end, 'AgeBP')
            # Interpolate concentration up to x number of rows
            # Limit_direction = backward, so it interpolates forward through time
            # Use the user-defined order for polynomial & spline interpolation
            interp_conc = temp['Sum 1.1-12'].interpolate(method = how, order = order_choice, limit = rows, limit_direction = 'backward')
            # Interpolate CPP
            interp_cpp = temp['CPP'].interpolate(method = how, order = order_choice, limit = rows, limit_direction = 'backward')    
            # Append interpolated values to final lists before moving to next interval
            final_conc = final_conc.append(interp_conc)
            final_cpp  = final_cpp.append(interp_cpp)
            
        # Return series of interpolated concentration & CPP data
        return final_conc, final_cpp 
#%%
# Get interpolation setup dataframe
interp = setup_interp(cfa)
#%%
# Interpolation #1
linear_conc, linear_cpp = interpolate_cfa(cfa, interp)
# Remove interpolated values which are too small or too large
linear_conc = linear_conc[linear_conc >= 0]
linear_conc = linear_conc[linear_conc <  300]
linear_cpp  = linear_cpp[linear_cpp >= 0]
linear_cpp  = linear_cpp[linear_cpp <  70]
# Copy interpolated data into new dataframe
cfa_linear = cfa.copy()
cfa_linear['Sum 1.1-12'] = linear_conc
cfa_linear['CPP'] = linear_cpp
# Delete extra columns
cfa_linear = cfa_linear.loc[:, 'AgeBP':'CPP'].copy()

# Interpolation #2
poly2_conc, poly2_cpp = interpolate_cfa(cfa, interp)
# Remove interpolated values which are too small or too large
poly2_conc = poly2_conc[poly2_conc >= 0]
poly2_conc = poly2_conc[poly2_conc <  300]
poly2_cpp  = poly2_cpp[poly2_cpp >= 0]
poly2_cpp  = poly2_cpp[poly2_cpp <  70]
# Copy interpolated data into new dataframe
cfa_poly2 = cfa.copy()
cfa_poly2['Sum 1.1-12'] = poly2_conc
cfa_poly2['CPP'] = poly2_cpp
# Delete extra columns
cfa_poly2 = cfa_poly2.loc[:, 'AgeBP':'CPP'].copy()

# Interpolation #3
poly3_conc, poly3_cpp = interpolate_cfa(cfa, interp)
# Remove interpolated values which are too small or too large
poly3_conc = poly3_conc[poly3_conc >= 0]
poly3_conc = poly3_conc[poly3_conc <  300]
poly3_cpp  = poly3_cpp[poly3_cpp >= 0]
poly3_cpp  = poly3_cpp[poly3_cpp <  70]
# Copy interpolated data into new dataframe
cfa_poly3 = cfa.copy()
cfa_poly3['Sum 1.1-12'] = poly3_conc
cfa_poly3['CPP'] = poly3_cpp
# Delete extra columns
cfa_poly3 = cfa_poly3.loc[:, 'AgeBP':'CPP'].copy()
#%%
# Comparing linear, 2nd-order polynomial, 3rd-order spline
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax1)
fig.subplots_adjust(hspace = 0.05)

ax1.plot(cfa_linear['AgeBP'], cfa_linear['Sum 1.1-12'], color = 'orange')
ax1.plot(cfa_poly2['AgeBP'],  cfa_poly2['Sum 1.1-12'],  color = 'blue')
ax1.plot(cfa_poly3['AgeBP'],  cfa_poly2['Sum 1.1-12'],  color = 'green')
ax1.plot(cfa['AgeBP'], cfa['Sum 1.1-12'], color = 'black')
ax1.set_ylabel('Conc')
mylabels = ['Linear', '2nd-Order Poly', '3rd-Order Poly', 'None']
ax1.legend(mylabels, loc = 'upper left')
ax1.axes.get_xaxis().set_visible(False)  

ax2.plot(cfa_linear['AgeBP'], cfa_linear['CPP'], color = 'orange')
ax2.plot(cfa_poly2['AgeBP'],  cfa_poly2['CPP'],  color = 'blue')
ax2.plot(cfa_poly3['AgeBP'],  cfa_poly3['CPP'],  color = 'green')
ax2.plot(cfa['AgeBP'], cfa['CPP'], color = 'black')
ax2.set_ylim(0, 80)
ax2.set_ylabel('CPP')
ax2.set_xlabel('Age BP')

#%%
# Interpolation setup 2

#poly3_conc, poly3_cpp = interpolate_cfa(cfa, interp)
#cfa_poly3 = cfa.copy()
#cfa_poly3['Sum 1.1-12'] = poly3_conc
#cfa_poly3['CPP'] = poly3_cpp
## Delete extra columns
#cfa_poly3 = cfa_poly3.loc[:, 'AgeBP':'CPP'].copy()
#
#spline2_conc, spline2_cpp = interpolate_cfa(cfa, interp)
#cfa_spline2 = cfa.copy()
#cfa_spline2['Sum 1.1-12'] = spline2_conc
#cfa_spline2['CPP'] = spline2_cpp
## Delete extra columns
#cfa_spline2 = cfa_spline2.loc[:, 'AgeBP':'CPP'].copy()

#%%
# Comparing 2nd- and 3rd- order poly and spline interpolations
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#ax2 = fig.add_subplot(212, sharex=ax1)
#fig.subplots_adjust(hspace = 0.05)
#
#ax1.plot(cfa_poly['AgeBP'], cfa_poly['Sum 1.1-12'], color = 'blue')
#ax1.plot(cfa_poly3['AgeBP'], cfa_poly3['Sum 1.1-12'], color = 'orange')
#ax1.plot(cfa_spline2['AgeBP'], cfa_spline2['Sum 1.1-12'], color = 'red')
#ax1.plot(cfa_spline['AgeBP'], cfa_spline['Sum 1.1-12'], color = 'green')
#ax1.plot(cfa['AgeBP'], cfa['Sum 1.1-12'], color = 'black')
#ax1.set_ylabel('Conc')
#mylabels = ['2nd-Order Poly', '3rd-Order Poly', '2nd-Order Spline', '3rd-Order Spline', 'None']
#ax1.legend(mylabels)
#ax1.axes.get_xaxis().set_visible(False)  
#
#ax2.plot(cfa_poly['AgeBP'], cfa_poly['CPP'], color = 'blue')
#ax2.plot(cfa_poly3['AgeBP'], cfa_poly3['CPP'], color = 'orange')
#ax2.plot(cfa_spline2['AgeBP'], cfa_spline2['CPP'], color = 'red')
#ax2.plot(cfa_spline['AgeBP'], cfa_spline['CPP'], color = 'green')
#ax2.plot(cfa['AgeBP'], cfa['CPP'], color = 'black')
#ax2.set_ylim(0, 80)
#ax2.set_ylabel('CPP')
#ax2.set_xlabel('Age BP')
#%%