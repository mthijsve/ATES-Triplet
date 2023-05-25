
from pathos.threading import ThreadPool as Pool
import os
from functions_triplet import *
from flow_function_Triplet import *


def model_run(Qyh, Qyc, T_h, Thmin):
    WD = os.getcwd()
    exec(open('setup.py').read())

    exec(open('PySeawaTriplet.py').read())

    loss2 = dnmh
    
    return loss2

def parallel_model_run( Qyc, Qyh,T_h, Thmin, threads=None):

     # Set number of threads (cores) used for parallel run and map threads
    if threads is None:
        pool = Pool()
    else:
        pool = Pool(nodes=threads)
    # Run parallel models
    results = pool.map(
        model_run,
        Qyc,Qyh,T_h,Thmin,
    )
    return results