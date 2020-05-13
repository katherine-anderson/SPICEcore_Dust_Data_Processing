# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:38:11 2020

@author: Katie
"""


import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt
os.chdir('C:\\Users\\katie\\OneDrive\\Documents\\SPICE\\Data')

wout_hump = pd.read_csv('Bad_CFA_Phase2_NoHumps_4-28.csv', index_col = 'Unnamed: 0')
with_hump = pd.read_csv('Bad_CFA_Phase2_2020-03-13.csv', index_col = 'Unnamed: 0')

with_hump = with_hump.sort_index()
wout_hump = wout_hump.sort_index()

#%%

hump_rows = with_hump[(with_hump['Error Type'] == 'PSD Hump Anomaly')].index.values.tolist()

test = wout_hump[wout_hump.index.isin(hump_rows) == True].copy()

#%% Plot

which = 'Sum 1.1-12'

plt.figure()
plt.scatter(with_hump.loc[hump_rows, 'AgeBP'], with_hump.loc[hump_rows, which])
# plt.figure()
plt.scatter(test['AgeBP'], test[which])
plt.xlabel('AgeBP')
plt.ylabel(which)
plt.legend(['Humps', 'Humps Identified as Outliers'], fontsize = 15)

