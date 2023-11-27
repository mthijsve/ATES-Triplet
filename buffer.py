import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')

def extract_variables(filename):
    
    filename = os.path.basename(filename)
    start = filename.find("Qh")
    filename = filename[start:]
    parts = filename.split("_")

    Qh = int(parts[0][2:])/10
    Qc = int(parts[1][2:])/10
    injectionT = float(parts[2][10:])
    Thmin = float(parts[3].split(".")[0][5:])/10
    return Qh, Qc, injectionT, Thmin

WD = os.getcwd()
output_path = os.path.join(WD, "output")
output_B_path = os.path.join(WD, "Buffercheck output")

filenames = []
filenames_B = []

def get_filenames(path,list):
    for root,dirs,files in os.walk(path):
        for file in files:
            if file.endswith(".csv"):
                list.append(os.path.join(root,file))
    return list

filenames_B = get_filenames(output_B_path,filenames_B)
filenames = get_filenames(output_path,filenames)

# Extract the endings from filenames_B
endings_B = [os.path.basename(filename).split('Qh')[1] for filename in filenames_B if 'Qh' in filename]

# Keep only filenames in filenames that have a matching ending in filenames_B
filenames = [filename for filename in filenames if any(os.path.basename(filename).endswith(ending) for ending in endings_B)]

# Initialize an empty DataFrame
def get_results(filenames):
    dfs = [] #list to store Dataframes
    for filename in filenames:
        Qyh, Qyc, injectionT, Thmin = extract_variables(filename)
        df = pd.read_csv(filename)
        ehta = 1 / df['Efficiency_h'].iloc[2]
        V_h = df['W0_Vin'].abs().sum() + df['W0_Vout'].abs().sum()
        V_c = df['W1_Vin'].abs().sum() + df['W1_Vout'].abs().sum()
        V_b = df['W2_Vin'].abs().sum() + df['W2_Vout'].abs().sum()
        dfs.append(pd.DataFrame({'Qyh': [Qyh], 'Qyc': [Qyc], 'injectionT': [injectionT], 'Thmin': [Thmin], 'ehta': [ehta], 'V_h': [V_h], 'V_c': [V_c], 'V_b': [V_b]}))

    df_results = pd.concat(dfs, ignore_index=True)
    # Convert units
    df_results[['Qyh', 'Qyc']] /= 1e12  # J to TJ
    df_results[['V_h', 'V_c', 'V_b']] /= 1e6  # m3 to million m3

    # Calculate additional columns
    df_results['V_h_norm'] = df_results['V_h'] / df_results['Qyh']
    df_results['V_c_norm'] = df_results['V_c'] / df_results['Qyc']
    df_results['Q_in'] = df_results['Qyh'] / df_results['ehta']
    df_results['m^2'] = df_results['Q_in'] / 0.0018
    df_results['Cutoff_T'] = 22 + (df_results['injectionT'] - 22) * df_results['Thmin']

    return df_results

df_NB = get_results(filenames)
df_NB['New'] = 0
df_B = get_results(filenames_B)
df_B['New'] = 1

df_Combi = pd.concat([df_NB,df_B], ignore_index=True)
df_Combi.to_csv('combicheck.csv')

y_min = df_Combi['m^2'].min()
y_max = df_Combi['m^2'].max()
unique_injectionT = df_Combi['injectionT'].unique()

fig, axs = plt.subplots(len(unique_injectionT), 1, figsize=(10, 15), sharex=True)

for ax, injectionT in zip(axs, unique_injectionT):
    df_injectionT = df_Combi[df_Combi['injectionT'] == injectionT]
    grouped = df_injectionT.groupby(['Qyh', 'New'])
    for (Qyh, New), group in grouped:
        ax.plot(group['Cutoff_T'], group['m^2'], marker='o', linestyle='-', label=f'Qyh={Qyh}, New={New}')
    ax.set_xlabel('Cutoff Temperature [°C]')
    ax.set_ylabel('m^2 [-]')
    ax.legend()
    ax.set_ylim(y_min, y_max)  # Set y-axis limits

plt.tight_layout()
plt.show()