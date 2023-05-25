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
path = os.path.join(WD,out)
files = glob.glob(path + "/**/*.csv")


#%%
for file in files:
    df = pd.read_csv(file)
    file2 = os.path.basename(file)
    string = os.path.splitext(file2)[0]
    string = string.replace('Run_output__Test','')    
    Vol= abs(df['W0_Vin'].sum())+abs(df['W0_Vout'].sum())+abs(df['W1_Vin'].sum())+abs(df['W1_Vout'].sum())+abs(df['W2_Vin'].sum())+abs(df['W2_Vout'].sum())
    print(string,':',Vol)
    
