# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 18:55:33 2018

@author: Aaron and Heather

MASTER SPICE FILE

"""

import numpy as np
import pandas as pd
import scipy.interpolate as interp1d
import matplotlib.pyplot as plt
import scipy as sp
import os
import math
from datetime import date
import statistics as stats

os.chdir('C:\\Users\\Aaron\\Desktop\\PhD Stuff\\SPICE\\Data')

cfa = pd.read_excel('RAW_CFA.xlsx', header = 0)
#cfa = pd.read_excel('CFA_734_removed_breaks_questions.xlsx',header=0)
annual_depths = pd.read_excel('SPICE_Timescale_Feb2018.xlsx', header=0)
annual_depths['Years'] = annual_depths.index+1
AR = pd.read_excel('SPICE_accumulation.xlsx')
SPICE_CFA_Dates = pd.read_csv('SPICE_CFA_Dates.csv')



cfa1 =SPICE_CFA_Dates.copy()

del cfa1['date_time1']
#SPICE734a = SPICE734.copy()
#ECM = SPICE734a[['Depth(m)','ECM']]
error_dict ={}
error2=[]

'''
machine error 2 (highlights areas where there is a verticle line in the data and 
there are two depth measurements assigned.  Values are nan'd and highlightd - 
checked and works for top few meters)
'''
#error_m2top = []
#error_m2bot = []

fig,ax = plt.subplots()
ax.plot(cfa1['Depth(m)'],cfa1['1'],c='k')
ax.semilogy()
ax1=ax.twinx()
ax2=ax.twinx()
ax3=ax.twinx()
ax4=ax.twinx()
ax5=ax.twinx()
ax1.semilogy()
ax2.semilogy()
ax3.semilogy()
ax4.semilogy()
ax5.semilogy()


for i in range(1,len(cfa1['Depth(m)'])-1):
    print(i)
    if cfa1['Depth(m)'][i] == cfa1['Depth(m)'][i-1]:
#        print(cfa1['Depth(m)'][i])
#        error_m2top.append(cfa1['Depth(m)'][i])
#        error_m2bot.append(cfa1['Depth(m)'][i-1])
        cfa1.loc[i,2:-2]=np.nan
        cfa1.loc[i-1,2:-2]=np.nan

        
#error_dict['machine_error_2']=error2
print(cfa1)
ax1.plot(cfa1['Depth(m)'],cfa1['1'],c='r')
'''
removing bubbles - removes dust values that occur during a bubble (where the slope of the ECM decrease and increase between
three points - sent to Bess)
'''

#bubble=[]
for i in range(1,len(cfa1['Depth(m)'])-1):
    print(i)
    x2 = cfa1['Depth(m)'][i]
    x1 = cfa1['Depth(m)'][i-1]
    y2 = cfa1['ECM'][i]
    y1 = cfa1['ECM'][i-1]
    x3 = cfa1['Depth(m)'][i+1]
    y3 = cfa1['ECM'][i+1]
    slope = (y2-y1)/(x2-x1)
    slope2 = (y3-y2)/(x3-x2)
    if (slope <= -25) & (slope2 >= 25):#make sure this is the correct slope for bubbles!!!!!!!!!!!!!!
#        bubble.append(cfa1['Depth(m)'][i])
        cfa1.iloc[i,2:-2]=np.nan
print(cfa1)
ax2.plot(cfa1['Depth(m)'],cfa1['1'],c='g')
#error_dict['bubble']=bubble
#top =7
#bot =25
#cfa_tp = cfa.loc[(cfa['Depth(m)'] > top) & (cfa['Depth(m)'] < bot)]
#cfa1_tp = cfa1.loc[(cfa1['Depth(m)'] > top) & (cfa1['Depth(m)'] < bot)] #making depth range smaller for example

'''
machine error 3 (Nans data where depth value decrease - checked and works)
'''
#error3=[]
#error_top = []
#error_bot = []

for i in range(1,len(cfa1['Depth(m)'])-1):
    print(i)
    if cfa1['Depth(m)'][i] <= cfa1['Depth(m)'][i-1]:
        print(cfa1['Depth(m)'][i])
#        error_top.append(cfa1['Depth(m)'][i])
#        error_bot.append(cfa1['Depth(m)'][i-1])
        cfa1.loc[i,2:-2]=np.nan
        cfa1.loc[i-1,2:-2]=np.nan

print(cfa1)
ax3.plot(cfa1['Depth(m)'],cfa1['1'],c='m')
#error_dict['machine_error_3'] = error3


'''
machine error 4 (Finds data that is too close together)
'''

#error_m4top = []
#error_m4bot = []

for i in range(1,len(cfa1['Depth(m)'])-1):
    print(i)
    if cfa1['Depth(m)'][i] - cfa1['Depth(m)'][i-1] <= 0.00002: #change this value to change threshold for error
        print(cfa1['Depth(m)'][i])
#        error_m4top.append(cfa1['Depth(m)'][i])
#        error_m4bot.append(cfa1['Depth(m)'][i-1])
        cfa1.loc[i,2:-2]=np.nan
        cfa1.loc[i-1,2:-2]=np.nan

print(cfa1)
ax4.plot(cfa1['Depth(m)'],cfa1['1'],c='orange')

'''
removing negative values, Flow Rate, ECM and Particle Size bins.  Any row with
a negative value has been nan'd to not influence the PSD - checked and works
'''

cfa1 = cfa1.set_index(['Depth(m)','melt_date1'])

#nan_df_breaks = cfa1.loc[cfa1.isnull().any(axis=1),:]

cfa1[cfa1< 0] = np.nan

cfa1.loc[cfa1.isnull().any(axis=1),:] = np.nan
#nan_df_errors = cfa1.loc[cfa1.isnull().any(axis=1),:]
print(cfa1)
ax5.plot(cfa1.index,cfa1['1'],c='grey')

'''
Contamination
Need to highlight areas that have area associated with a core break, spike in 
dust data, and followed by a decreasing back to background concentrations

'''

'''
Highlight data that is too far apart (ie. 195m)
'''




'''
Writing cleaned cfa1 file to a text file.  Date of cleaning to be included in
title
'''
del cfa1['Unnamed: 0']
cfa1.to_csv('cleaned_CFA_4SEPT2018.csv')
        
