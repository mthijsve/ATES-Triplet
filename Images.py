import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import flopy
from pathlib import Path
import os

def extract_variables(filename):

    #used chatgpt to write this function
    filename = os.path.basename(filename)
    start = filename.find("Qh")
    filename = filename[start:]
    parts = filename.split("_")

    Qh = int(parts[0][2:])/10
    Qc = int(parts[1][2:])/10
    injectionT = float(parts[2][10:])
    Thmin = float(parts[3].split(".")[0][5:])/10
    return Qh, Qc, injectionT, Thmin

WD = Path.cwd()
output_path = os.path.join(WD, "output")

filenames = []
for root,dirs,files in os.walk(output_path):
    for file in files:
        if file.endswith(".csv"):
            filenames.append(os.path.join(root,file))

#empty list to store results
Qh_list = []
Qc_list = []
injectionT_list = []
Thmin_list = []
ehta_list = []
V_h_list = []
V_c_list = []
V_b_list = []

for filename in filenames:
    Qyh,Qyc,injectionT,Thmin = extract_variables(filename)
    df = pd.read_csv(filename)
    ehta = 1 / df['Efficiency_h'].iloc[2]
    V_h = df['W0_Vin'].abs().sum() + df['W0_Vout'].abs().sum()
    V_c = df['W1_Vin'].abs().sum() + df['W1_Vout'].abs().sum()
    V_b = df['W2_Vin'].abs().sum() + df['W2_Vout'].abs().sum()
    
    # Only save lines where Qyh is equal to Qyc
    if Qyh == Qyc:
        V_h_list.append(V_h)
        V_c_list.append(V_c)
        V_b_list.append(V_b)
        Qh_list.append(Qyh)
        Qc_list.append(Qyc)
        injectionT_list.append(injectionT)
        Thmin_list.append(Thmin)
        ehta_list.append(ehta)

# df_good = df_results[(df_results['ehta'] > 0.70) & (df_results['ehta'] < 0.75)]
# df_medium = df_results[(df_results['ehta'] > 0.45) & (df_results['ehta'] < 0.55)]
# df_bad = df_results[(df_results['ehta'] > 0.25) & (df_results['ehta'] < 0.35)]

#Define a colormap based on the injectionT values


data = {'Qyh':Qh_list, 'Qyc':Qc_list, 'injectionT':injectionT_list, 'Thmin':Thmin_list, 'ehta':ehta_list, 'V_h':V_h_list, 'V_c':V_c_list, 'V_b':V_b_list}
df_results = pd.DataFrame(data)
df_results['Qyh'] = df_results['Qyh']/1e12
df_results['V_h'] = df_results['V_h']/1e6
df_results['V_c'] = df_results['V_c']/1e6
df_results['Qyc'] = df_results['Qyc']/1e12
df_results['V_h_norm']=df_results['V_h']/df_results['Qyh']
df_results['V_c_norm']=df_results['V_c']/df_results['Qyc']
cmap = plt.cm.get_cmap('viridis')

# Define a dictionary of markers based on the Thmin values
markers = {0.1:'.', 0.3:'o', 0.5:'^', 0.7:'*',0.9:'P'}

plt.figure(0)
for thmin, group in df_results.groupby('Thmin'):
    plt.scatter(group['Qyh'], group['ehta'], c=group['injectionT'], cmap=cmap, marker=markers[thmin], label=f'Cut-off ={thmin}')
plt.xlabel('yearly demand [TJ]')
plt.ylabel('efficiency [-]')
#plt.ylabel(r'total pumped volume hot well [million m$^3$]')
#plt.title('The efficiency for all combinations of different parameters')
plt.colorbar().set_label('Injection Temperature [°C]')
plt.legend()
plt.show()

plt.figure(1)
for thmin, group in df_results.groupby('Thmin'):
    plt.scatter(group['V_h'], group['ehta'], c=group['injectionT'], cmap=cmap, marker=markers[thmin], label=f'Cut-off ={thmin}')
plt.xlabel(r'total pumped volume hot well [million m$^3$]')
plt.ylabel('efficiency [-]')
plt.colorbar().set_label('Injection Temperature [°C]')
plt.legend()
plt.show()

plt.figure(2)
for thmin, group in df_results.groupby('Thmin'):
    plt.scatter(group['V_h_norm'], group['ehta'], c=group['injectionT'], cmap=cmap, marker=markers[thmin], label=f'Cut-off ={thmin}')
plt.xlabel(r'total pumped volume hot well [million m$^3$]')
plt.ylabel('efficiency [-]')
plt.colorbar().set_label('Injection Temperature [°C]')
plt.legend()
plt.show()

plt.figure(3)
for thmin, group in df_results.groupby('Thmin'):
    plt.scatter(group['V_c'], group['Qyc'], c=group['Thmin'], cmap=cmap, label=f'Cut-off ={thmin}')
plt.xlabel(r'total pumped volume cold well [million m$^3$]')
plt.ylabel('efficiency [-]')
plt.colorbar().set_label('yearly demand cold well [TJ]')
plt.legend()
plt.show()

