import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import flopy
from pathlib import Path
import os
import re

def extract_variables(filename):
    # Extract the numeric part of the string after 'Qh', 'Qc', 'injectionT', 'CP', and 'RP' (github copilot helped me with this one)
    Qh = int(re.search('Qh(\d+)', filename).group(1))/10
    Qc = int(re.search('Qc(\d+)', filename).group(1))/10
    injectionT = float(re.search('injectionT(\d+)', filename).group(1))
    CP = float(re.search('CP(\d+)', filename).group(1))/10
    RP = float(re.search('RP(\d+)', filename).group(1))/10

    return Qh, Qc, injectionT, CP, RP

WD = Path.cwd()
output_path = os.path.join(WD, "Paper1_output")

filenames = []
for root,dirs,files in os.walk(output_path):
    for file in files:
        if file.endswith(".csv"):
            filenames.append(os.path.join(root,file))

#empty list to store results
Qh_list = []
Qc_list = []
injectionT_list = []
Cutoffper_list = []
Returnper_list = []
ehta_list = []
V_h_list = []
V_c_list = []
V_b_list = []
T_b_max_list = []
T_b_avg_list = []

for filename in filenames:

    Qyh,Qyc,injectionT,CP,RP = extract_variables(filename)
    df = pd.read_csv(filename)
    if Qyc != 0:
        ehta = 1 / df['Efficiency_h'].iloc[2]
        V_h = df['W0_Vin'].abs().sum() + df['W0_Vout'].abs().sum()
        V_c = df['W1_Vin'].abs().sum() + df['W1_Vout'].abs().sum()
        V_b = df['W2_Vin'].abs().sum() + df['W2_Vout'].abs().sum()
        T_b_max = df['W2_T_mf_out'].max()
        T_b_avg = df['W2_T_mf_out'].mean()
        V_h_list.append(V_h)
        V_c_list.append(V_c)
        V_b_list.append(V_b)
        Qh_list.append(Qyh)
        Qc_list.append(Qyc)
        injectionT_list.append(injectionT)
        Cutoffper_list.append(CP)
        Returnper_list.append(RP)

        ehta_list.append(ehta)
        T_b_max_list.append(T_b_max)
        T_b_avg_list.append(T_b_avg)



data = {'Qyh':Qh_list, 'Qyc':Qc_list, 'injectionT':injectionT_list, 'CP':Cutoffper_list, 'RP':Returnper_list, 'ehta':ehta_list, 'V_h':V_h_list, 'V_c':V_c_list, 'V_b':V_b_list, 'T_b_max':T_b_max_list, 'T_b_avg':T_b_avg_list}
df_results = pd.DataFrame(data)
df_results['Qyh'] = df_results['Qyh']/1e12  #J to TJ
df_results['V_h'] = df_results['V_h']/1e6   #m3 to million m3
df_results['V_b'] = df_results['V_b']/1e6   #m3 to million m3
df_results['V_c'] = df_results['V_c']/1e6   #m3 to million m3
df_results['Qyc'] = df_results['Qyc']/1e12  #J to TJ
df_results['V_h_norm']=df_results['V_h']/df_results['Qyh']
df_results['V_c_norm']=df_results['V_c']/df_results['Qyc']
df_results['Q_in'] = df_results['Qyh']/df_results['ehta']
df_results['m^2'] = df_results['Q_in']/0.0018
df_results['T_b_max'] = T_b_max_list
df_results['T_b_avg'] = T_b_avg_list

df_results.to_csv('results.csv')

cmap = plt.cm.get_cmap('viridis')

# Define a dictionary of markers based on the Thmin values
markers = {0.5:'o', 0.7:'*',0.9:'^'}

plt.figure(0)
for cp, group in df_results.groupby('CP'):
    plt.scatter(group['V_h'], group['ehta'], c=group['injectionT'], cmap=cmap, marker=markers[cp], label=f'Cut-off ={cp}')
# plt.xlabel('solar heat collector area [m$^2$]')
# plt.ylabel('efficiency [-]')
plt.ylabel(r'total pumped volume hot well [million m$^3$]')
plt.title('The efficiency for all combinations of different parameters')
plt.colorbar().set_label('Injection Temperature [°C]')
plt.legend()
plt.show()

plt.figure(1)
plt.plot((df_results['Qyh']/df_results['Qyc']),df_results['T_b_max'],'o')
plt.xlabel('Qyh/Qyc')
plt.ylabel('Maximum temperature Buffer [°C]')
plt.show()

plt.figure(2)
plt.plot((df_results['Qyh']/df_results['Qyc']),df_results['T_b_avg'],'o')
plt.xlabel('Qyh/Qyc')
plt.ylabel('Average temperature Buffer [°C]')
plt.show()