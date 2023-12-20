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
Cutoffper=[0.5, 0.7]#, 0.9]
Returnper=[0.5, 0.7]#, 0.9]

combinations = itertools.product(Qyh, Qyc, injectionT, Cutoffper, Returnper)

def calibration(bounds, Qyh, Qyc, injectionT, Cutoffper, Returnper):
    def find_zero(Z, start, end, step):
        temp = 1
        values = np.arange(start, end, step)
        for i in values:
            if Z(i) == 0:
                temp = i
                break
        return temp

    def model_run_wrapper(i):
        return abs(pst.Modelrun(i, Qyh, Qyc, injectionT, Cutoffper, Returnper))

    temp = find_zero(model_run_wrapper, bounds[0], bounds[1], 0.5)
    if temp is not None:
        temp = find_zero(model_run_wrapper, max(1,temp - 0.4), temp + 0.5, 0.1)
    if temp is not None:
        temp = find_zero(model_run_wrapper, max(1,temp - 0.09), temp + 0.1, 0.01)
    print(temp)
    return temp


def run_calibration_combination(combination):
    Qyh, Qyc, injectionT, Cutoffper, Returnper = combination
    mini = calibration([1, 15], Qyh, Qyc, injectionT, Cutoffper, Returnper)
    print(mini)

if __name__ == '__main__':
    with ProcessPoolExecutor() as executor:
        executor.map(run_calibration_combination, combinations)

    print('Calibration took', (time.time() - st) / 60, 'minutes')

