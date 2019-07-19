# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 10:13:42 2018

@author: Aaron

Highlights and deals with errors that are associated with Ethanol-140 contamination

Ethanol-140 contains particles in the fluid and is found in the core breaks (stick ends, breaks and 
annealed fractures).

The goal of this script is to highlight those areas of likely contamination in the data and remove them

Steps
1) Upload machine error cleaned dataset after going through the Master File and the 300 error removal
for SPICE

2) Highlight core break (stick ends, breaks) regions to see if Ethanol-140 matches with those sections

3) Find areas that peak have a local maximum above a certain value (3 std) with stratled by local
minimums (background, median values)
"""
from __future__ import print_function

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.integrate import simps
import scipy.signal 
import scipy.stats as ss
from numpy import trapz
import plotly.plotly as py
import plotly.graph_objs as go


#os.chdir('C:\\Users\\Aaron\\Desktop\\PhD Stuff\\SPICE\\Python_Programs\sympy-master\sympy-master\sympy\calculus')
#
#import sympy

os.chdir('C:\\Users\\Aaron\\Desktop\\PhD Stuff\\SPICE\\Data')

cfa = pd.read_csv('Master_CFA_11June2018.csv')
annual_depths = pd.read_excel('SPICE_Timescale_Feb2018.xlsx', header=0)
annual_depths['Years'] = annual_depths.index+1
breaks = pd.read_excel('SPICE_core_breaks_734.xlsx', sheetname = 'half cm')

cfa1 = cfa.copy()
breaks1 = breaks.copy()

breaks1 = breaks1/100 #converting from cm to m
breaks1 = breaks1 + 5.15 #setting breaks to top of core

cfa1['Concentration']=np.sum(cfa1.iloc[:, 5:], axis=1)*1000

#years_interp = pd.Series(np.interp(cfa1['Depth(m)'],annual_depths['Depth'],annual_depths['Age'] ))
#cfa1.insert(1,'AgeBP',years_interp)
#cfa1 = cfa1.set_index(cfa1['AgeBP'])
cfa2 = cfa1.drop(['datetime','melt_date1','Flow Rate','ECM'], axis = 1)

'''
Using numpy gradient to identify signal changes
'''
cfa2['gradient'] = np.gradient(cfa2['Concentration'])
cfa2 = cfa2.set_index('Depth(m)')

#cfa3 = cfa2.dropna(axis = 0)
#
#cfa3['peaks'] = scipy.signal.find_peaks(cfa3['gradient'])


'''
Highlight core breaks
'''

fig,ax = plt.subplots(figsize = (16,12))
ax.plot(cfa1['Depth(m)'], cfa1['Concentration'])
       
for i in range(len(breaks1)):
    ax.axvspan(breaks1['Lower'][i],breaks1['Upper'][i],facecolor='r', alpha=0.5)

ax.set_title('South Pole Concentration', fontsize = 20)
ax.set_xlabel('Depth(m)', fontsize = 16)
ax.set_ylabel('Concentration (particles/mL)', fontsize = 16)
plt.savefig('Core_breaks.png')
#        
#        
#
#d = 10000
#with PdfPages('Core_Breaks.pdf') as pdf:
#    
#    for i in range(0,len(cfa1),d):
#        print(i)
#        fig,ax = plt.subplots(figsize = (20,12))
#        ax.plot(cfa1['Depth(m)'][i:d], cfa1['Concentration'][i:d])
#        try:
#            ax.set_title('{}-{}'.format(round(cfa1['Depth(m)'][i]),round(cfa1['Depth(m)'][d])), fontsize = 16)
#        except:
#            ax.set_xlim([round(cfa1['Depth(m)'][i]),round(cfa1['Depth(m)'].max())])
#        
#        
##               
#        for j in range(len(breaks1)):
#            ax.axvspan(breaks1['Lower'][j],breaks1['Upper'][j],facecolor='r', alpha=0.5)
#        try:    
#            ax.set_xlim([cfa1['Depth(m)'][i],cfa1['Depth(m)'][d]])
#        except:
#            ax.set_xlim([cfa1['Depth(m)'][i],cfa1['Depth(m)'].max()])
#        d = d+10000
##    ax.set_title('South Pole Concentration', fontsize = 20)
#        ax.set_xlabel('Depth(m)', fontsize = 16)
#        ax.set_ylabel('Concentration (particles/mL)', fontsize = 16)
#        ax.tick_params(axis ='both', labelsize = 16)
#        pdf.savefig(fig)
#        plt.close()
##        
#     
         
        
'''
Using python built integrals for idenfyng sections of core with large integrals

Script that identifies areas within y number data points
'''

##Compute the area using the comppsite Trapezoidal rule
#area = trapz(cfa1['Concentration'], dx=5)
#print("area =", area)
#
## Compute the area using the composite Simpson's rule.
#area = simps(cfa1['Concentration'], dx=5)
#print("area =", area)


y = 10 #interval of data points per integral
tpz = []

for x in range(0,len(cfa1),10):
    tpz.append(trapz(cfa1['Concentration'][x:y]))
    y = y+10



stdev = np.std(tpz)
median = np.median(tpz)
error = stdev + median

y = 10
contamination=[] #make this a dataframe 
integrals = [] #making this into a list of a list and then transfer to Dataframe (contamination)
for x in range(0,len(cfa1),10):
    print('{} {}'.format(x,y))
    integrals.append(trapz(cfa1['Concentration'][x:y]))
    if trapz(cfa1['Concentration'][x:y]) >= error:  # if area increases outside standard deviation of all intergrals
         contamination.append([trapz(cfa1['Concentration'][x:y]),cfa1['Depth(m)'][x], cfa1['Depth(m)'][y]]) #need to get integrals associated depth
#        cfa1['Concentration'][x] = np.nan
    y = y + 10



contam_depths = pd.DataFrame(contamination)
contam_depths.columns = ['Integral', 'top_Depth(m)', 'bot_Depth(m)']
#contam_depths.to_csv('Contam_Integrals.csv')


'''
Highlighting sections of the core where the integral is outside of the standard deviation of median
'''
props = props = dict(boxstyle='round', facecolor='wheat', alpha= 1.0)

d = 10000
with PdfPages('Core_Breaks_Integrals.pdf') as pdf:
    
    for i in range(0,len(cfa1),d):
        print(i)
        fig,ax = plt.subplots(figsize = (20,12))
        ax.plot(cfa1['Depth(m)'][i:d], cfa1['Concentration'][i:d])
        try:
            ax.set_title('{}-{}'.format(round(cfa1['Depth(m)'][i]),round(cfa1['Depth(m)'][d])), fontsize = 16)
        except:
            ax.set_xlim([round(cfa1['Depth(m)'][i]),round(cfa1['Depth(m)'].max())])
        
        
#               
        for j in range(len(breaks1)):
            ax.axvspan(breaks1['Lower'][j],breaks1['Upper'][j],facecolor='r', alpha=0.5)
        try:    
            ax.set_xlim([cfa1['Depth(m)'][i],cfa1['Depth(m)'][d]])
        except:
            ax.set_xlim([cfa1['Depth(m)'][i],cfa1['Depth(m)'].max()])
        
        
        for k in range(len(contam_depths)):
            ax.axvspan(contam_depths['top_Depth(m)'][k],contam_depths['bot_Depth(m)'][k],facecolor='g', alpha=0.5)
        try:    
            ax.set_xlim([cfa1['Depth(m)'][i],cfa1['Depth(m)'][d]])
        except:
            ax.set_xlim([cfa1['Depth(m)'][i],cfa1['Depth(m)'].max()])
            
        d = d+10000
#   ax.set_title('South Pole Concentration', fontsize = 20)
        ax.set_xlabel('Depth(m)', fontsize = 16)
        ax.set_ylabel('Concentration (particles/mL)', fontsize = 16)
        ax.tick_params(axis ='both', labelsize = 16)
        text_input='Red = Stick Breaks \nGreen = Contamination'
        ax.text(0.95, 0.95, text_input, verticalalignment='top', transform=ax.transAxes, bbox = props, size = 16)
        pdf.savefig(fig)
        plt.close()


'''
Removing Core Breaks and Integral Spikes

'''

#
for a in range(len(cfa1['Depth(m)'])):
    for b,c in zip(breaks1['Lower'], breaks1['Upper']):
        if (cfa1['Depth(m)'][a] >= b) & (cfa1['Depth(m)'][a] <= c):
            print(cfa1['Depth(m)'][a])
            cfa1['Depth(m)'][a] = np.nan
#
#
#for a in range(len(cfa1['Depth(m)'])):
#    for b,c in zip(contam_depths['top_Depth(m)'], contam_depths['top_Depth(m)']):
#        if (cfa1['Depth(m)'][a] >= b) & (cfa1['Depth(m)'][a] <= c):
#            print(cfa1['Depth(m)'][a])
#            cfa1['Depth(m)'][a] = np.nan


'''
Identifying Regions by finding peaks, with width and prominence
'''

median = np.median(cfa2['Concentration'])
std = np.std(cfa2['Concentration'])
error = median + std

peaks = scipy.signal.find_peaks(cfa2['Concentration'], prominence = error, width = (1,10))

peaks_index = peaks[0][:]
dictionary = peaks[1:6][:]
contam_dictionary = dictionary[0]
contam_peaks = pd.DataFrame.from_dict(contam_dictionary)
contam_peaks = contam_peaks.set_index(peaks_index)

fig,ax = plt.subplots(2,figsize = (20,12))
ax[0].plot(cfa1['Depth(m)'], cfa1['Concentration'])
ax[1].plot(cfa1.index, cfa1['Concentration'])
ax[0].tick_params(labelsize = 20)
ax[1].tick_params(labelsize = 20)

for i in range(len(cfa1.index)):
    for m,n in zip(contam_peaks['left_bases'], contam_peaks['right_bases']):
        if (cfa1.index[i] >= m) & (cfa1['Depth(m)'][i] <= n):
            ax.axvspan(m,n, facecolor = 'red', alpha = 0.5)




