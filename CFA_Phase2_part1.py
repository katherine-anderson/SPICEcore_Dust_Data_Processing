# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 14:14:30 2019

@author: Aaron
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz

os.chdir(r'C:\Users\Aaron\Desktop\PhD Stuff\SPICE\Cleaning Procedure')


cfa = pd.read_csv('..\Data\Cleaned_CFA_Phase1_2019-08-07.csv')
cfa = cfa.drop(['Unnamed: 0', 'Flow Rate', 'ECM'], axis = 1)

y = 2 #interval of data points per integral

tpz = []
index = []
cpp_tpz = []

for x in range(0,len(cfa),2):

    tpz.append(trapz(cfa['Sum 1.1-12'][x: (x+y)]))
    cpp_tpz.append(trapz(cfa['CPP'][x: (x+y)]))
    index.append(x)
  
integrals = pd.DataFrame({'index': index, 'integral':tpz, 'cpp_integral': cpp_tpz})   

stdev = np.nanstd(tpz)
median = np.nanmedian(tpz)
error = (2 * stdev) + median

cpp_stdev = np.nanstd(cpp_tpz)
cpp_median = np.nanmedian(cpp_tpz)
cpp_error = (2 * cpp_stdev) + cpp_median

contamination=[] #make this a dataframe 
cpp_contamination=[] #make this a dataframe 

integrals_bad = integrals[(integrals['integral'] >= error) & (integrals['cpp_integral'] >= cpp_error)]

other = pd.DataFrame()
other['index'] = integrals_bad['index'] + 1

all_contam_int = pd.concat([other, integrals_bad], axis = 0)
all_contam_int = all_contam_int.sort_values('index')


dust_rows = pd.Series(cfa[cfa['Dust Event?'] == True].index.values.astype('float'))
volc_rows = pd.Series(cfa[cfa['Volcanic Event?'] == True].index.values.astype('float'))

a = all_contam_int.loc[all_contam_int['index'].isin(dust_rows)]
b = all_contam_int.loc[all_contam_int['index'].isin(volc_rows)]

all_contam_int = all_contam_int.set_index('index')
a = a.set_index('index')
b = b.set_index('index')

to_drop = pd.concat([all_contam_int, a, b], axis = 0)
to_drop = to_drop.reset_index()
bad_data = to_drop.drop_duplicates('index', keep = False)



















#for x in range(0,len(cfa), 2):
#
#    if trapz(cfa1['conc'].iloc[x:y]) >= error:  # if area increases outside standard deviation of all intergrals
#         contamination.append([trapz(cfa['conc'].iloc[x:y]),cfa['Depth(m)'].iloc[x], cfa['Depth(m)'].iloc[y]]) #need to get integrals associated depth
#    #        cfa1['Concentration'][x] = np.nan
#    y = y + 2



#contam_depths = pd.DataFrame(contamination)
#contam_depths.columns = ['Integral', 'top_Depth(m)', 'bot_Depth(m)']
#
#integral_contam_depths = contam_depths.copy()

#integral_contam_depths.to_csv('integral_contam_depths.csv')