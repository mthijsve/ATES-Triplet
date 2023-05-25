#%%
import os
import sys
import shutil 
import time
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import flopy
import flopy.modflow as mf
import flopy.mt3d as mt3
import flopy.seawat as swt
import flopy.utils.binaryfile as bf

from functions_triplet import *
from flow_function_Triplet import *
from pathlib import Path
#%%

WD = os.getcwd()
out = 'output'
imagesOut = 'Images'
path = os.path.join(WD,out)
files = glob.glob(path + "/**/*.csv")
impath = os.path.join(path,imagesOut)
#%%
def AllDataFromFiles(files):
    df = pd.read_csv(files)
    file = os.path.basename(files)
    string = os.path.splitext(file)[0]
    string = string.replace('Run_output__test','')
    data = string.split('_')
    list = ['t','x','z']
    for i in range(len(list)):
        data[i] = data[i].replace(list[i],'')
        if data[i] == '05':
            data[i] = 0.5
        else:
            data[i] = float(data[i])
            
    dfexp = pd.DataFrame({'dt':[data[0]],'x':[data[1]],'z':[data[2]]})
    dfexp['dnm'] = df['dnmh'].sum()
    dfexp['dnm'] = dfexp['dnm'].multiply(data[0])
    dfexp['V_h_in'] = df['W0_Vin'].sum()
    dfexp['V_h_out']= df['W0_Vout'].sum()
    dfexp['V_c'] = df['W1_Vin'].sum()+df['W1_Vout'].sum()
    dfexp['V_m'] = df['W2_Vin'].sum()+df['W2_Vout'].sum()
    dfexp['T_h'] = df['W0_T_mf_out'].mean()
    dfexp['Tempsum'] = df['W0_T_mf_out'].sum()*dfexp['dt']
    dfexp['T_c'] = df['W1_T_mf_out'].mean()
    dfexp['T_m'] = df['W2_T_mf_out'].mean()
    return dfexp
df = pd.DataFrame()
for i in files:
    df = pd.concat([df,AllDataFromFiles(i)])
df.sort_values(by=['dt'], axis = 0, inplace=True)
df.set_index(keys=['dt'],drop=False,inplace=True)
timestep = df['dt'].unique().tolist()
vls = ['dnm','V_h_in','V_h_out','V_c','V_m','T_h','Tempsum','T_c','T_m']
for i in timestep:
    v='dnm'
    x = [0.5,1,2,5,10]
    z = [0.5,1,2,5,10]
    pv = pd.pivot_table(df.loc[df['dt']==i],values=v,index='z',columns='x')
    plt.figure()
    ax = sns.heatmap(pv,vmin = 41,vmax = 45,annot=True,fmt=".2f")
    ax.set_xticks([0.5,1.5,2.5,3.5,4.5],['0.5','1','2','5','10'])
    ax.set_yticks([0.5,1.5,2.5,3.5,4.5],['0.5','1','2','5','10'])
    ax.set_title('timestep:'+str(i)+' - '+v)
    ax.invert_yaxis()
    plt.savefig(os.path.join(imagesOut,str(i)+str(v)+'.png'))
    
#%% [E] plots
"""days = np.zeros(rl)
for i in range(rl):
     days[i] = perlen*(i+1)
colors = ['r','k','b', 'g', 'c', 'y','grey', 'violet', 'dodgerblue','firebrick','coral','yellow', 'lightgreen', 'orange','firebrick','cyan']
                   
# ''' Temperatures and discharge'''
fig, (ax0, ax1, ax2)  = plt.subplots(nrows=3, sharex=True)
ax0.set_title('Discharge of wells',fontsize=10)
im0 = ax0.plot(days,Run_output.loc[:,'W0_Vin'] + Run_output.loc[:,'W0_Vout'], color=colors[0], label=well_obj_list[0].type)
im1 = ax1.plot(days,Run_output.loc[:,'W1_Vin'] + Run_output.loc[:,'W1_Vout'], color=colors[1], label=well_obj_list[1].type)
im2 = ax2.plot(days,Run_output.loc[:,'W2_Vin'] + Run_output.loc[:,'W2_Vout'], color=colors[2], label=well_obj_list[2].type)
ax1.set_ylabel('well discharge [m3/day]',fontsize=10)

plt.grid(True, 'minor', lw=1, c='grey')   
plt.legend() 

fig, (ax1) = plt.subplots(nrows=1, sharex=True)
ax1.set_title('Temperature  wells [C]',fontsize=10)
for i in range(nW):
    im1 = ax1.plot(days,Run_output.loc[:,'W'+str(i)+'_T_mf_out'], color=colors[i], label=well_obj_list[i].type),
ax1.set_ylabel('Temperature filter screens [C]',fontsize=10) 
ax1.set_facecolor('lightgrey')
plt.grid(True, 'major', lw=1, c='w') 
plt.show()
"""