import os
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import flopy.modflow as mf
import flopy.mt3d as mt3
import flopy.seawat as swt
import flopy.utils.binaryfile as bf
import PySeawaTriplet as pst
import time
from concurrent.futures import ProcessPoolExecutor

from functions_triplet import *
from flow_function_Triplet import *


st = time.time()
WD = os.getcwd()

Qyh =[0.5e12, 1e12]#, 2e12, 5e12, 10e12]
Qyc =[0.5e12, 1e12]#, 2e12, 5e12, 10e12]
injectionT =[50, 60]#, 70, 80, 90]


combinations = itertools.product(Qyh, Qyc, injectionT)

def calibration(bounds, Qyh, Qyc, injectionT, corr_w, corr_c):
    def find_zero(Z, start, end, step):
        temp = 1
        values = np.arange(start, end, step)
        for i in values:
            if Z(i) == 0:
                temp = i
                break
        return temp

    def model_run_wrapper_dnmc(corr_c):
        return abs(pst.Modelrun(corr_w, corr_c, Qyh, Qyc, injectionT)[1])

    def model_run_wrapper_dnmh(corr_w):
        return abs(pst.Modelrun(corr_w, corr_c, Qyh, Qyc, injectionT)[0])
    # First calibrate to dnmc
    corr_c = find_zero(model_run_wrapper_dnmc, bounds[0], bounds[1], 0.5)
    if corr_c is not None:
        corr_c = find_zero(model_run_wrapper_dnmc, max(1,corr_c - 0.4), corr_c + 0.5, 0.1)
    if corr_c is not None:
        corr_c = find_zero(model_run_wrapper_dnmc, max(1,corr_c - 0.09), corr_c + 0.1, 0.01)

    # Then calibrate to dnmh
    corr_w = find_zero(model_run_wrapper_dnmh, bounds[0], bounds[1], 0.5)
    if corr_w is not None:
        corr_w = find_zero(model_run_wrapper_dnmh, max(1,corr_w - 0.4), corr_w + 0.5, 0.1)
    if corr_w is not None:
        corr_w = find_zero(model_run_wrapper_dnmh, max(1,corr_w - 0.09), corr_w + 0.1, 0.01)

    print(corr_w, corr_c)
    return corr_w, corr_c

def run_calibration_combination(combination):
    Qyh, Qyc, injectionT = combination
    mini = calibration([1, 5], 2e12, 2e12, 60, 1,1)
    print(mini)

if __name__ == '__main__':
    with ProcessPoolExecutor() as executor:
        executor.map(run_calibration_combination, combinations)

    print('Calibration took', (time.time() - st) / 60, 'minutes')

