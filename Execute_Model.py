import os
import itertools
from functions_triplet import *
from flow_function_Triplet import *
from calibration import model_run, parallel_model_run
import numpy as np

WD = os.getcwd()

Qyh = [0, 0.5e12, 1e12, 2e12, 5e12, 10e12]
Qyc = [0, 0.5e12, 1e12, 2e12, 5e12, 10e12]
injectionT = [40, 50, 60, 70, 80, 90]
Thmin= [0.5, 0.6, 0.7, 0.8, 0.9, .95]

combinations = itertools.product(Qyh, Qyc, injectionT, Thmin)


for Qyh, Qyc, injectionT, Thmin in combinations:
    # Set name for run
    Qyhls = [Qyh]
    Qycls = [Qyc]
    T_h = [injectionT] 
    Thminls = [Thmin]
    results = parallel_model_run(Qycls, Qyhls, T_h, Thminls)
    #write code that calibrates the model so that dnmhs = 0 and corr_ws is the lowest it can be, as well as dnmh = 0 and corr_h is the lowest it can be
    



