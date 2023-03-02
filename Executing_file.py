'''



'''
#%%import packages and set directories
import os

from functions_triplet import *
from flow_function_Triplet import *


WD = r'C:\Users\msvanesch\Documents\08-Python_Scripts\PySeawaTriplet'
os.chdir(WD)       # Set the proper working directory (WD)

'''Define input files'''
swtexe_name = 'swt_v4x64.exe' 
wellfile = 'wells.csv'
demandfile = 'demand.csv'     
subsurffile = 'sub-surface.csv'     
monitoringfile = 'monitoring.csv'

#%%create the objects needed
well_obj_list = createobjfromCSV(PyWell,wellfile) #creates well object from  data in csv file
form_obj_list = createobjfromCSV(PySystem,subsurffile)
mon_obj_list = createobjfromCSV(PySystem,monitoringfile) # creates list of monitoring points in your grid where you can track the temperature

nmon = len(mon_obj_list)
nform = len(form_obj_list)
nW=len(well_obj_list)         # numer of wells
T_amb = form_obj_list[0].s1   # Ambient Temperature of the subsurface
PEFF = form_obj_list[0].por   # porosity [-]

#%% B: settings for model run

'''Set name for run''' 
naam = 'test'
case = '2'
name = naam + case

dirs = WD+'/output/'+name                                                      
if os.path.exists(dirs)==False:
    os.makedirs(dirs)

'''main time setting parameters'''
perlen = 1    # (DAYS)
years = 5        # set the minimum #years that a system will run
ppy = 365/perlen                                                       

startwinter = 1 # simulation starts in winter (1) or summer (0). for distribution of flows of wells and surface temperature
        
nper = 1      # Number of Modflow stress periods per coupled simulation period
nstp = 1      # Number of Modflow time steps per stress period
timprs = np.array([perlen])
nprs = len(timprs)

repeatrun = 0     # repeatrun starts the run with the last temperature file of the previous run as initial conditions
finaltimestep = 365/perlen*years-1

'''Model flow input'''
flowtype = 2     # choise between 0 and 4, see readme
gwflow_x = -0.   # Groundwater flowvelocity in x direction [m/y]
gwflow_y = -0.   # Groundwater flowvelocity in y direction [m/y]

densflow = 1     # turn densityflow off or on by switching the VDF packages off or on  
WaterDensON = 1  # Calculation of parameters standard =1

'''correction factors determine the recharge of the wells by multiplying the needed flow to meet the demand by the factor'''
corr_w = 1.5      # correction factor hot well
corr_1 = corr_w   # startup correction factor, only applied to hot well.
corr_c = 1.0      # correction factor cold well

Thmin = 5          #Thmin is the temperature difference in the hot well between the injection temperature and the cutoff temperature

for i in well_obj_list:       # Update each active Python well object with the temperature and head at its grid location
    if i.type == 'warm':
        print(i.T_inj)
        cutofftemp_h = i.T_inj-Thmin
        
cutofftemp_c = T_amb-2   # Temperature when the heating turns off, thus no more discharge from the hot well

T_room = 20                 # Room temperature in building, used for return temperature groundwater after heat exchanger
T_loss_building = 2         # degrees that the return temperature deviates from room temperature. The higher, the bigger the losses in the HVAC system. 



'''grid settings '''
AXI = 1              # axial symmetric grid or not. | 1=AXI, 0 = 3D |
LIN = 0              # Linear (1) or logarithmic (0) cell sizes buildup around the center of the modelgrid 
ICBUND_Up = -1       # TOP:: -1 = boundary has constant temperature (NO Flow, YES Temperature boundary), 1 = NO Flow, No Tempertature boundary
ICBUND_Down = -1     # BOTTOM: -1 = boundary has constant temperature (NO Flow, YES Temperature boundary), 1 = NO Flow, No Tempertature boundary
OutsideAirBound = 0  # 1 = ON, 0=OFF, if ON: the temperature boundary at model top is adjusted following outside air temperature 
steady = True        # Switch for steady-state/transient Modflow simulation

dmin = 1              # smallest cel size at well [m]
dmin_bound = 200      # total distance from well with 'dmin' sized cells [m]
dmax = 200            # largest cell size at model boundary [m]
aroundAll = 1000      # normal=1500 [m] size of grid around well.
nstep = 20            # minimum number of steps that gridfunctions must add  to acquire the dmax size requirement [-]
grid_extents = None   # example:[-300,300,-300,300] #set model boundaries [m]
dz = 1                # vertical gridlayer thickness [m]  important to syncronize with layer thicknesses in csv file!
#%%

exec(open('PySeawaTriplet.py').read())

