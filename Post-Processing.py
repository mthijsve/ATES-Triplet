import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import flopy
from pathlib import Path
import os
import re
import seaborn as sns


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
Qinc_list = []
Qinh_list = []
checklist = []

T_cinj = 5
T_crp = 18

for filename in filenames:

    Qyh,Qyc,injectionT,CP,RP = extract_variables(filename)
    df = pd.read_csv(filename)
    CuT = 22 + (injectionT - 22) * CP
    ReT = 22 + (CuT - 22) * RP
    if Qyc != 0:
        
        # Calculate Wh_T_r
        df['Wh_T_r'] = np.where(
            (df['W0_Vin'] > 0) & (df['W2_Vout'] > 0),
            (df['W1_Vout']*18 + df['W2_Vout']*df['W2_T_mf_out'])/(df['W1_Vout']+df['W2_Vout']),
            np.where(
                df['W0_Vin'] > 0,
                18,
                np.nan
            )
        )
        
        # Calculate Wc_T_r
        df['Wc_T_r'] = np.where(
            (df['W1_Vin'] > 0) & (df['W2_Vout'] > 0),
            (df['W0_Vout']*ReT + df['W2_Vout']*df['W2_T_mf_out'])/(df['W0_Vout']+df['W2_Vout']),
            np.where(
                df['W1_Vin'] > 0,
                ReT,
                np.nan
            )
        )
        
        df['Qinh'] = 4183000 * df['W0_Vin'] * (df['Wh_T_r'] - injectionT)
        df['Qinc'] = 4183000 * df['W1_Vin'] * (df['Wc_T_r'] - T_cinj)
        #for each timestep add the correct flow and temperature for regeneration.
        #calculate the net flow from the return and buffer mixture. 
        ehta = 1 / df['Efficiency_h'].iloc[2]
        V_h = df['W0_Vin'].abs().sum() + df['W0_Vout'].abs().sum()
        V_c = df['W1_Vin'].abs().sum() + df['W1_Vout'].abs().sum()
        V_b = df['W2_Vin'].abs().sum() + df['W2_Vout'].abs().sum()
        Qinh = df['Qinh'].abs().sum()
        Qinc = df['Qinc'].abs().sum()
        T_b_max = df['W2_T_mf_out'].max()
        T_b_avg = df['W2_T_mf_out'].mean()
                # Select the last five years of data (5 years * 365 days/year)
        last_five_years = df.tail(5 * 30*12)

        # Count the number of 'True' statements in the 'check' column
        num_true = last_five_years['dnmh'].sum()
        checklist.append(num_true)
        V_h_list.append(V_h)
        V_c_list.append(V_c)
        V_b_list.append(V_b)
        Qh_list.append(Qyh)
        Qc_list.append(Qyc)
        injectionT_list.append(injectionT)
        Cutoffper_list.append(CP)
        Returnper_list.append(RP)
        Qinc_list.append(Qinc)
        Qinh_list.append(Qinh)
        ehta_list.append(ehta)
        T_b_max_list.append(T_b_max)
        T_b_avg_list.append(T_b_avg)

data = {'Qyh':Qh_list, 'Qyc':Qc_list, 'injectionT':injectionT_list, 'CP':Cutoffper_list, 'RP':Returnper_list, 'ehta':ehta_list, 'V_h':V_h_list, 'V_c':V_c_list, 'V_b':V_b_list, 'T_b_max':T_b_max_list, 'T_b_avg':T_b_avg_list, 'Qinc':Qinc_list, 'Qinh':Qinh_list,'check':checklist}

df_results = pd.DataFrame(data)
df_results['Qyh'] = df_results['Qyh']/1e12  #J to TJ
df_results['Qinc'] = df_results['Qinc']/8/1e12  #J to TJ
df_results['Qinh'] = df_results['Qinh']/8/1e12  #J to TJ
df_results['V_h'] = df_results['V_h']/1e6   #m3 to million m3
df_results['V_b'] = df_results['V_b']/1e6   #m3 to million m3
df_results['V_c'] = df_results['V_c']/1e6   #m3 to million m3
df_results['Qyc'] = df_results['Qyc']/1e12  #J to TJ
df_results['V_tot']=df_results['V_h']+df_results['V_c']+df_results['V_b']
df_results['Qpump'] = df_results['V_tot']*0.54/8 #[TJ]
df_results['V_h_norm']=df_results['V_h']/df_results['Qyh']
df_results['V_c_norm']=df_results['V_c']/df_results['Qyc']
df_results['V_b_norm']=df_results['V_b']/df_results['Qyh']
df_results['Q_in'] = df_results['Qyh']/df_results['ehta']
df_results['m^2'] = df_results['Q_in']/0.0018
df_results['C_m^2'] = df_results['Qinh']/0.0018
df_results['DC_cap'] = df_results['Qyc']*1e12/(803*3600)/1e6 #MW
df_results['C_DC_cap'] = df_results['Qinc']/8*1e12/(803*3600)/1e6 #MW
df_results['T_b_max'] = T_b_max_list

#df_results = df_results[df_results['check'] == 0]
#df_results = df_results.drop(columns=['check'])
#df_results = df_results[df_results['ehta'] > 0.1]
#df_results = df_results[df_results['ehta'] < 0.87]


# Count the number of unique entries in the 'Qyc' column
coldcheck = df_results['Qyc'].value_counts()
heatcheck = df_results['Qyh'].value_counts()
cpcheck = df_results['CP'].value_counts()
rpcheck = df_results['RP'].value_counts()
injcheck = df_results['injectionT'].value_counts()

print('cold:',coldcheck)
print('heat:',heatcheck)
print('cp:',cpcheck)
print('rp:',rpcheck)
print('inj:',injcheck)

#df_results.to_csv('results.csv')
#COP DC 20f
#space .1 m2/kw

# Create five plots for total pumped volume V_h
plt.figure(figsize=(12, 8))


sns.set_palette('colorblind')
colorblind_palette = sns.color_palette()
cyan_colorblind = colorblind_palette[0]  # cyan color in colorblind palette

# Box Plot 1
plt.subplot(2, 3, 1)
sns.boxplot(x='Qyh', y='C_m^2', data=df_results, color=cyan_colorblind)
plt.xlabel('yearly heating demand (TJ)')
plt.ylabel('Solar Collectors (m²)')

# Box Plot 2
plt.subplot(2, 3, 2)
sns.boxplot(x='Qyc', y='C_m^2', data=df_results, color=cyan_colorblind)
plt.xlabel('Yearly cooling demand (TJ)')
plt.ylabel('Solar Collectors (m²)')

# Box Plot 3
plt.subplot(2, 3, 3)
sns.boxplot(x='injectionT', y='C_m^2', data=df_results, color=cyan_colorblind)
plt.xlabel('injection temperature (°C)')
plt.ylabel('Solar Collectors (m²)')

# Box Plot 4
plt.subplot(2, 3, 4)
sns.boxplot(x='CP', y='C_m^2', data=df_results, color=cyan_colorblind)
plt.xlabel('Cutoff percentage (%)')
plt.ylabel('Solar Collectors (m²)')

# Box Plot 5
plt.subplot(2, 3, 5)
sns.boxplot(x='RP', y='C_m^2', data=df_results, color=cyan_colorblind)
plt.xlabel('return percentage (%)')
plt.ylabel('Solar Collectors (m²)')
plt.tight_layout()
plt.savefig(r'C:\Users\msvanesch\Documents\06-geschreven_werk\Paper1\Images\SCResultbox.png')
plt.show()