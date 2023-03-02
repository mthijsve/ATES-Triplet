import os

from functions_triplet import *
from flow_function_Triplet import *

WD = os.getcwd()
print(WD)

perlent = [0.5,1,2,5,10,20,30]    # (DAYS)
dmint = [0.5,1,2,5,10]              # smallest cel size at well [m]
dzt = [0.5,1,2,5,10]                # vertical gridlayer thickness [m]  important to syncronize with layer thicknesses in csv file!

for i in range(len(perlent)):
    for k in range(len(dmint)):
        for j in range(len(dzt)):
            perlen = perlent[i]
            dmin = dmint[k]
            dz = dzt[j]           
            dmin_bound = 200      # total distance from well with 'dmin' sized cells [m]
            dmax = 200            # largest cell size at model boundary [m]
            aroundAll = 1000      # normal=1500 [m] size of grid around well.
            nstep = 20            # minimum number of steps that gridfunctions must add  to acquire the dmax size requirement [-]
            '''Set name for run''' 
            naam = 'test'
            case = str(perlen + dmin + dz)
            name = naam + case
            exec(open('setup.py').read())
            exec(open('PySeawaTriplet.py').read())