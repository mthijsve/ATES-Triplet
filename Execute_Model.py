import os
from pathlib import Path

from functions_triplet import *
from flow_function_Triplet import *

import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import flopy
import flopy.modflow as mf
import flopy.mt3d as mt3
import flopy.seawat as swt
import flopy.utils.binaryfile as bf
import PySeawaTriplet as pst
import time

st = time.time()
WD = os.getcwd()

Qyh =[2e12] #[0, 0.5e12, 1e12, 2e12, 5e12, 10e12]
Qyc =[2e12] #[0, 0.5e12, 1e12, 2e12, 5e12, 10e12]
injectionT =[45] #[40, 50, 60, 70, 80, 90]
Thmin=[.656] #[0.5, 0.6, 0.7, 0.8, 0.9, .95]

combinations = itertools.product(Qyh, Qyc, injectionT, Thmin)

def objective(x,Qyh,Qyc,injectionT,Thmin):
    print(x)
    Z = abs(pst.Modelrun(x[0],x[1],Qyh,Qyc,injectionT,Thmin))
    return Z


for Qyh, Qyc, injectionT, Thmin in combinations:
    #minimize triplet model

    mini = minimize(objective,[1.0,1.0],args=(Qyh,Qyc,injectionT,Thmin),method = 'Nelder-Mead',bounds=[(1,5),(1,4)])
    print(mini)

print('calibraten took', (time.time()-st)/60, 'minutes' )

