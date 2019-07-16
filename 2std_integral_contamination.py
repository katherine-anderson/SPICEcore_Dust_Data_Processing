# -*- coding: utf-8 -*-
"""
Created on Sun May 19 17:08:36 2019

@author: Aaron
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 15:55:43 2019

@author: Aaron
"""

import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz


os.chdir(r'C:\Users\Aaron\Desktop\PhD Stuff\SPICE\Data')

cfa = pd.read_csv('SPICE_raw_26FEB2019.csv')
volc = pd.read_excel('SPICE_Holocene_volcanic_events.xlsx')
measured_breaks = pd.read_excel('core_breaks_full.xlsx')
annual_depths = pd.read_excel('SPICE_Timescale_2_25_2019.xlsx', header=0)
annual_depths['Years'] = annual_depths.index+1

cfa = cfa.replace('#NAME?', np.nan)
cfa = cfa.astype('float')
cfa = cfa.drop(['conc', 'Cond', 'Flow Rate'], axis = 1)
cfa = cfa.set_index('Depth(m)')


cfa = cfa.replace([np.inf, -np.inf], np.nan) #removing infinites 
cfa = cfa.dropna() #dropping all nan values
cfa = cfa[cfa >= 0] #Getting rid of any concentrations below 0

conc = cfa.sum(axis = 1)
coursep = cfa[['4.5', '5.1', '5.7', '6.4', '7.2', '8.1', '9', '10', '12']].sum(axis = 1)
cpp = (coursep/conc) *100

cfa['conc'] = conc
cfa['cpp'] = cpp

'''
Removing errors from core
'''

md = np.nanmedian(cfa['conc'])
std = np.nanstd(cfa['conc'])
error = md + std

cfa = cfa[cfa['conc'] < error]

cfa = cfa[cfa.index <= 800]
measured_breaks = measured_breaks[measured_breaks['Depth (m)'] <= 800]



'''
Creating 5cm window of volcanic tie points as measure
of where the volcanoes are
'''

volc['tr'] = volc['depth(m)']-0.05
volc['br'] = volc['depth(m)']+0.05

#volc = volc[volc['Bot D (m)'] <= 800]

'''
Index is reset in order to calculate integrals underneath the concentration
profile
'''

cfa1 = cfa.reset_index()


y = 2 #interval of data points per integral
tpz = []

for x in range(0,len(cfa1),2):
    tpz.append(trapz(cfa1['conc'][x: (x+y)]))

stdev = np.nanstd(tpz)
median = np.nanmedian(tpz)
error = (2 * stdev) + median


y = 2
contamination=[] #make this a dataframe 
integrals = [] #making this into a list of a list and then transfer to Dataframe (contamination)

for x in range(0,len(cfa1), 2):
    print('{} {}'.format(x,y))
    integrals.append(trapz(cfa1['conc'].iloc[x:y]))
    if trapz(cfa1['conc'].iloc[x:y]) >= error:  # if area increases outside standard deviation of all intergrals
         contamination.append([trapz(cfa1['conc'].iloc[x:y]),cfa1['Depth(m)'].iloc[x], cfa1['Depth(m)'].iloc[y]]) #need to get integrals associated depth
    #        cfa1['Concentration'][x] = np.nan
    y = y + 2



contam_depths = pd.DataFrame(contamination)
contam_depths.columns = ['Integral', 'top_Depth(m)', 'bot_Depth(m)']

integral_contam_depths = contam_depths.copy()

integral_contam_depths.to_csv('integral_contam_depths.csv')


'''
Calculating integrals for CPP
'''




y = 2 #interval of data points per integral
cpp_tpz = []

for x in range(0,len(cfa1),2):
    cpp_tpz.append(trapz(cfa1['cpp'][x: (x+y)]))

cpp_stdev = np.nanstd(cpp_tpz)
cpp_median = np.nanmedian(cpp_tpz)
cpp_error = (2 * cpp_stdev) + cpp_median



y = 2
cpp_contamination=[] #make this a dataframe 
cpp_integrals = [] #making this into a list of a list and then transfer to Dataframe (contamination)

for x in range(0,len(cfa1),2):
#    while y <= 800:
    print('{} {}'.format(x,y))
    cpp_integrals.append(trapz(cfa1['cpp'].iloc[x:y]))
    if trapz(cfa1['cpp'].iloc[x:y]) >= cpp_error:  # if area increases outside standard deviation of all intergrals
         cpp_contamination.append([trapz(cfa1['cpp'].iloc[x:y]),cfa1['Depth(m)'].iloc[x], cfa1['Depth(m)'].iloc[y]]) #need to get integrals associated depth
#        cfa1['Concentration'][x] = np.nan
    y = y + 2



cpp_contam_depths = pd.DataFrame(cpp_contamination)
cpp_contam_depths.columns = ['Integral', 'top_Depth(m)', 'bot_Depth(m)']


cpp_contam_depths.to_csv('cpp_contam_depths.csv')


#'''
#Measureing integrals with ECM
#'''
#
#bubble=[]
#
#for i in range(1,len(cfa1['Depth(m)'])-1):
#
#    x2 = cfa1['Depth(m)'][i]
#    x1 = cfa1['Depth(m)'][i-1]
#    y2 = cfa1['ECM'][i]
#    y1 = cfa1['ECM'][i-1]
#    x3 = cfa1['Depth(m)'][i+1]
#    y3 = cfa1['ECM'][i+1]
#    slope = (y2-y1)/(x2-x1)
#    slope2 = (y3-y2)/(x3-x2)
#    
#    if (slope <= -25) & (slope2 >= 25):#make sure this is the correct slope for bubbles!!!!!!!!!!!!!!
#        bubble.append(cfa1['Depth(m)'][i])
#
#'''
#Check for bubbles
#'''
#
#
#fig,ax = plt.subplots(figsize = (25,12))
#ax.plot(cfa1['Depth(m)'], cfa1['ECM'], color = 'k')
#for i in bubble:
#    ax.axvspan(i - 0.01, i + 0.01, facecolor = 'orange', alpha = 0.75)#highlighting 1cm region around bubbles that could be affected
# 
#   


'''
Plotting Integrals with Depths
'''

fig,ax1 = plt.subplots(figsize = (20,12))
ax1.plot(cfa1['Depth(m)'], cfa1['conc'], color = 'k')
ax2 = ax1.twinx()
ax2.plot(cfa1['Depth(m)'], cfa1['cpp'], color = 'steelblue', alpha = 0.25)

for a,b in zip(contam_depths['top_Depth(m)'], contam_depths['bot_Depth(m)']):
    ax1.axvspan(a,b, facecolor = 'r', alpha = 0.5)

for c,d in zip(volc['tr'], volc['br']):
    ax1.axvspan(c,d, facecolor = 'b', alpha = 0.5) 

for e,f in zip(cpp_contam_depths['top_Depth(m)'], cpp_contam_depths['bot_Depth(m)']):
    ax1.axvspan(e,f, facecolor = 'green', alpha = 0.25)
    
for g,h,i in zip(measured_breaks['Minus 3 cm'], measured_breaks['Plus 3 cm'], measured_breaks['Depth (m)']):
    ax1.axvspan(g,h, facecolor = 'grey')
    ax1.axvline(i, color = 'k', linestyle = '--')
    
    

ax1.tick_params(labelsize = 20)
ax2.tick_params(labelsize = 20)
ax2.set_ylabel('CPP', fontsize = 25)
ax1.set_xlabel('Depth(m)', fontsize = 25)
ax1.set_ylabel('Conc. (part/uL)', fontsize = 25)
ax1.set_xlim(0,850)




'''
isolating just cpp and integral values hat overlap with each other
'''


uol = []
ulol = []
lol = []
llol = []

for a,b in zip(contam_depths['top_Depth(m)'], contam_depths['bot_Depth(m)']):
    for num,idx in enumerate(cpp_contam_depths['top_Depth(m)']):
        if (idx >= a) and (idx <= b):
#            print(cpp_contam_depths['top_Depth(m)'][num])
            uol.append(cpp_contam_depths['top_Depth(m)'][num])
            ulol.append(cpp_contam_depths['bot_Depth(m)'][num])


for a,b in zip(contam_depths['top_Depth(m)'], contam_depths['bot_Depth(m)']):
    for num,idx in enumerate(cpp_contam_depths['bot_Depth(m)']):
        if (idx >= a) and (idx <= b):
            print(cpp_contam_depths['top_Depth(m)'][num])
            lol.append(cpp_contam_depths['top_Depth(m)'][num])
            llol.append(cpp_contam_depths['bot_Depth(m)'][num])



upper_overlap = pd.DataFrame({'top_depth' : uol, 'bot_depth' : ulol})
lower_overlap = pd.DataFrame({'top_depth' : lol, 'bot_depth' : llol})

all_contam = pd.concat([upper_overlap, lower_overlap], axis = 0)
all_contam = all_contam.sort_values('bot_depth')

all_contam.to_csv('all_contam.csv')


fig,ax = plt.subplots(figsize = (20,12))
ax.plot(cfa1['Depth(m)'], cfa1['conc'], color = 'k')
ax.tick_params(labelsize = 20)
ax.set_ylabel('Conc µm/L', fontsize = 25)

for a,b in zip(upper_overlap['top_depth'], upper_overlap['bot_depth']):
    ax.axvspan(a,b, facecolor = 'r', alpha = 0.5)

for c,d in zip(lower_overlap['top_depth'], lower_overlap['bot_depth']):
    ax.axvspan(c,d, facecolor = 'g', alpha = 0.5)
    
#for e,f in zip(volc['tr'], volc['br']):
#    ax.axvspan(e,f, facecolor = 'b', alpha = 0.5)


'''
Identifying volcanic events so those depths can be kept by identifying a six yr
window were events could have had an impact
''' 

years_interp = pd.Series(np.interp(cfa1['Depth(m)'],annual_depths['Depth (m)'],annual_depths['Age (yr b 1950)'] ))
cfa1.insert(1,'AgeBP',years_interp)

cfa1 = cfa1.set_index('AgeBP')

volc['volc_end'] = 1950 - volc['volc_start'] + 6
volc['volc_start'] = 1950 - volc['volc_start']

volc['volc_end'] = 1950 - volc['volc_end']
volc['volc_start'] = 1950 - volc['volc_start']

volc = volc.astype(dtype = 'float')


upper_overlap['u_contam_yr'] = pd.Series(np.interp(upper_overlap['top_depth'],annual_depths['Depth (m)'],annual_depths['Age (yr b 1950)']))
upper_overlap['ul_contam_yr'] = pd.Series(np.interp(upper_overlap['bot_depth'],annual_depths['Depth (m)'],annual_depths['Age (yr b 1950)']))

lower_overlap['l_contam_yr'] = pd.Series(np.interp(lower_overlap['top_depth'],annual_depths['Depth (m)'],annual_depths['Age (yr b 1950)']))
lower_overlap['ll_contam_yr'] = pd.Series(np.interp(lower_overlap['bot_depth'],annual_depths['Depth (m)'],annual_depths['Age (yr b 1950)']))

#cfa2 = cfa1.reset_index()


volc_uoverlap = []
volc_uloverlap = []

for e,f in zip(upper_overlap['u_contam_yr'], upper_overlap['ul_contam_yr']):
    for num, g in enumerate (volc['volc_start']):
        if (g >= e) and (g <= f):
            print(g)
            volc_uoverlap.append(volc['volc_start'][num])
            volc_uloverlap.append(volc['volc_end'][num])

volc_loverlap = []
volc_lloverlap = []

for j,k in zip(lower_overlap['l_contam_yr'], lower_overlap['ll_contam_yr']):
    for num, g in enumerate (volc['volc_end']):
        if (g >= j) and (g <= k):
            print(g)
            volc_loverlap.append(volc['volc_start'][num])
            volc_lloverlap.append(volc['volc_end'][num])
            
volc_upper_overlap = pd.DataFrame({'end_yr' : volc_uoverlap, 'start_yr' : volc_uloverlap})
volc_lower_overlap = pd.DataFrame({'end_yr' : volc_loverlap, 'start_yr' : volc_lloverlap})

volc_upper_overlap1 = volc_upper_overlap.drop_duplicates(keep = 'first')
volc_lower_overlap1 = volc_lower_overlap.drop_duplicates(keep = 'first')


'''
Finding number of data points between volcanic events

combine volcanic dataframes

reseting index and will highlight start and end of volc year and save index value
use that value to find the number of data points in a volcanic time period
'''

volc_events = pd.concat([volc_upper_overlap1, volc_lower_overlap1], axis = 0)
volc_events = volc_events.sort_values(by=['end_yr'])

volc_events.to_csv('volc_events.csv')



volc_data = []


for t,y in zip(volc_events['end_yr'], volc_events['start_yr']):
    for r in (cfa1.index):
        if (r <= t) and (r >= y):
            volc_data.append(r)
print(len(volc_data))




'''
Plotting volcanic events
'''

fig,ax = plt.subplots(figsize = (20,12))
ax.plot(cfa1['conc'], color = 'k')
ax.tick_params(labelsize = 20)
ax.set_ylabel('Conc µm/L', fontsize = 25)
ax.set_xlabel('Years BP', fontsize = 25)

for aa,bb in zip(upper_overlap['u_contam_yr'], upper_overlap['ul_contam_yr']):
    ax.axvspan(aa,bb, facecolor = 'r', alpha = 0.5)

for cc,dd in zip(lower_overlap['l_contam_yr'], lower_overlap['ll_contam_yr']):
    ax.axvspan(cc,dd, facecolor = 'g', alpha = 0.5)
    
for ee,ff in zip(volc_upper_overlap1['start_yr'], volc_upper_overlap1['end_yr']):
    ax.axvspan(ee,ff, facecolor = 'b', alpha = 0.5)

for gg,hh in zip(volc_lower_overlap1['start_yr'], volc_lower_overlap1['end_yr']):
    print(gg,hh)
    ax.axvspan(gg,hh, facecolor = 'grey', alpha = 0.5)

volc_periods = pd.concat([volc_upper_overlap1, volc_lower_overlap1], axis = 0)
volc_periods = volc_periods.set_index('end_yr')
volc_periods = volc_periods.sort_index()
volc_periods = volc_periods.reset_index()




'''
Removing Contaminated Depths

Steps:

1) create new dataframe of contaminated depths from cfa and drop those rows

2) In contaminated depths in new dataframe create another dataframe of depths
    that overlap with volcanic depths
    
3) Concat those rows back into cleaned dataframe and sort index
    
Process will remove contaminated data and then put back contamination depths
associated with volcanic events.
'''


cfa['AgeBP'] = cfa1.index

contaminated_depths = []
contaminated_depths1 = []
contaminated_depths_u = []
contaminated_depths_l = []


#for depth in cfa.index:
for uc,ulc in zip(all_contam['top_depth'], all_contam['bot_depth']):
    cfalist = cfa.index[cfa.index.to_series().between(uc,ulc)].tolist()
#    contaminated_depths1.extend(cfalist)
    for i in cfalist:
        print(i)
        contaminated_depths.append(i)
#        if (depth >= uc) and (depth <= ulc):
#            contaminated_depths.append(depth)
            
contam_series = pd.Series(contaminated_depths)

contam_df = []

#for index, depth in enumerate(cfa.index):
#for contam in contam_series:

#    print(cfalist)
##        if depth == contam:
##            contam_df.append(cfa.loc[[depth]])
contam_df1 = cfa[cfa.index.to_series().isin(contam_series)]
#contam_df1 = pd.concat(contam_df)
contam_df2 = contam_df1.drop_duplicates(keep = 'first')
contam_df2 = contam_df2.reset_index()
contam_df2 = contam_df2.set_index('AgeBP')

not_bad = []

for idx1, year in enumerate(contam_df2.index):
    print(idx1)
    for v_sy, v_ey in zip(volc_upper_overlap['start_yr'], volc_upper_overlap['end_yr']) :
        if (year >= v_sy) and (year <= v_ey):
            not_bad.append(contam_df2.loc[[year]])
            
not_bad_yrs = pd.concat(not_bad)
not_bad_yrs1 = not_bad_yrs.drop_duplicates(keep = 'first') 

'''
Dropping all contaminated depth (including volcanic years) from cfa and
concating the volcanic years back in
'''

contam_df2 = contam_df2.set_index('Depth(m)')

trial = cfa.copy()


#for depth_1 in trial.index:
##    print(depth_1)
#    for contaminated in contam_df2.index:
#        if depth_1 == contaminated:
#            print(cfa.loc[[depth_1]])
#            cleaned = trial.drop([depth_1])
# 
#for contaminated in contam_df2.index:
##    trial1 = trial[trial.index != contaminated]
#    trial1 = trial[trial.index != contaminated]
    
trial1 = trial[~trial.index.isin(contam_df2.index)].dropna()
trial1 = trial1.reset_index()
trial1 = trial1.set_index('AgeBP')

clean_cfa = pd.concat([trial1, not_bad_yrs1], axis = 0)#not bad years needs to be by depth not YEAR!
clean_cfa1 = clean_cfa.sort_index()
clean_cfa1 = clean_cfa1.reset_index()
clean_cfa1 = clean_cfa1.set_index('Depth(m)')

fig,ax = plt.subplots(figsize = (20,12))
ax.plot(clean_cfa1['conc'], color = 'k', linewidth = 2)
ax.tick_params(labelsize = 20)
ax.set_xlabel('Years BP', fontsize = 25)
ax.set_ylabel('Conc um/L', fontsize = 25)

for a,b in zip(upper_overlap['top_depth'], upper_overlap['bot_depth']):
    ax.axvspan(a,b, facecolor = 'r', alpha = 0.5)

for c,d in zip(lower_overlap['top_depth'], lower_overlap['bot_depth']):
    ax.axvspan(c,d, facecolor = 'g', alpha = 0.5)

clean_cfa1.to_csv('contam_removed_volc_preserved_cfa.csv')

'''
Run Stats on cleaned data
'''


stats = clean_cfa1['conc'].describe()
'''
Stat Figures
'''

roll_std = clean_cfa1['conc'].rolling(5).std()

fig,ax = plt.subplots(2, figsize = (20,12), sharex = True)
ax[0].plot(roll_std, color = 'steelblue')
ax[0].set_ylabel('Rolling STDev', fontsize = 25)
ax[1].plot(clean_cfa1['conc'], 'k' )
ax[1].set_ylabel('Conc (particles/uL)', fontsize = 25)
ax[1].set_xlabel('Depth(m)', fontsize = 25)
ax[0].tick_params(labelsize = 20)
ax[1].tick_params(labelsize = 20)
ax[0].axvline(160, color = 'k', linestyle = '--', linewidth = 2)
ax[1].axvline(160, color = 'k', linestyle = '--', linewidth = 2)


u = roll_std[roll_std.index >= 160]
l = roll_std[roll_std.index <= 160]






