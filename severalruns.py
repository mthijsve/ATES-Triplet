import os

from functions_triplet import *
from flow_function_Triplet import *

WD = os.getcwd()
print(WD)

perlent = [10]            # (DAYS)
dmint = [0.5,1,2,5,10]              # smallest cel size at well [m]
dzt = [0.5,1,2,5,10]               # vertical gridlayer thickness [m]  important to syncronize with layer thicknesses in csv file!

dmin_bound = 200      # total distance from well with 'dmin' sized cells [m]
dmax = 200            # largest cell size at model boundary [m]
aroundAll = 1500      # normal=1500 [m] size of grid around well.
nstep = 20            # minimum number of steps that gridfunctions must add  to acquire the dmax size requirement [-]

"""
perlen = 1
dmin = 1
dz = 1
name = 'testt1_x1_z1'

exec(open('setup.py').read())
exec(open('contourplot.py').read())

"""
for i in range(len(perlent)):
    perlen = perlent[i]

    for j in range(len(dmint)):
        dmin = dmint[j]

        for k in range(len(dzt)):
            dz = dzt[k]        

            '''Set name for run''' 
            p = str(perlen).replace('.','')
            x = str(dmin).replace('.','')
            z = str(dz).replace('.','')

            naam = 'test'
            case = 't' + p +'_x' + x +'_z' + z
            name = naam + case
            print(name)
            exec(open('setup.py').read())
            exec(open('PySeawaTriplet.py').read())

