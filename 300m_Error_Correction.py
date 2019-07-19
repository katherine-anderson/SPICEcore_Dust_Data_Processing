# -*- coding: utf-8 -*-
"""
Created on Sun May 13 11:32:00 2018

@author: Aaron
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
#cfa = pd.read_excel('CFA_734_removed_breaks_questions.xlsx',header=0)
#annual_depths = pd.read_excel('SPICE_Timescale_Feb2018.xlsx', header=0)
#annual_depths['Years'] = annual_depths.index+1
#AR = pd.read_excel('SPICE_accumulation.xlsx')
SPICE_CFA_Dates = pd.read_csv('cleaned_CFA_11June2018.csv')

cfa1 =SPICE_CFA_Dates.copy()
#SPICE734a = SPICE734.copy()
#ECM = SPICE734a[['Depth(m)','ECM']]
error_dict ={}
error2=[]

'''
Wierd Section that returned low values at 300m
'''
#moving section up to neighboring values
#1) take median of run day from good areas? Melted on 7/19/2016 take median of 7/18 and 7/20
#2) multiple by x orders of magnitude?

#upload melt_log and append dates.  Then take mean of the median values of the neighboring 
# days and multiple by the section below.

cfa1['datetime']=pd.to_datetime(cfa1['melt_date1'])
cfa1=cfa1.set_index('datetime')
fig,ax = plt.subplots(2,sharex=True)

Day1 = cfa1.loc['2016-07-18 00:00:00']
ax[0].plot(Day1['Depth(m)'],Day1['1'], c = 'b')
ax[1].plot(Day1['Depth(m)'],Day1['1'], c = 'b')
Day1 = Day1.iloc[:,4:].median()
Day2 = cfa1.loc['2016-07-19 00:00:00']
ax[0].plot(Day2['Depth(m)'],Day2['1'], c = 'orange')
Day2m = Day2.iloc[:,4:].median()
Day3 = cfa1.loc['2016-07-20 00:00:00']
ax[0].plot(Day3['Depth(m)'],Day3['1'], c = 'g')
ax[1].plot(Day3['Depth(m)'],Day3['1'], c = 'g')
Day3 = Day3.iloc[:,4:].median()
ax[0].semilogy()
ax[1].semilogy()


for i in Day1.index:
    d1=Day1[i]/Day2m[i]
    d3=Day3[i]/Day2m[i]
    Day2[i]=Day2[i]*stats.mean([d1,d3])
    
    print(stats.mean([d1,d3]))
ax[1].plot(Day2['Depth(m)'],Day2['1'], c = 'orange')


#cfa1.loc['2016-07-19 00:00:00']*stats.mean([d1,d3])
cfa1.loc['2016-07-19 00:00:00']=Day2
fig,ax = plt.subplots()
ax.plot(cfa1['Depth(m)'],cfa1['1'])
ax.semilogy()

'''
Output file with date in file name
'''

cfa1.to_csv('Master_CFA_11June2018.csv')
