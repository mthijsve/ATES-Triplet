import numpy as np
from functions_triplet import createobjfromCSV,PySystem,PyWell

#%%import packages and set directories
'''Define input files'''
swtexe_name = 'swt_v4x64.exe' 
wellfile = 'wells.csv'
demandfile = 'demand.csv'     
subsurffile = 'sub-surface.csv'     
monitoringfile = 'monitoring.csv'

'''main time setting parameters'''
perlen = 1    # (DAYS)
startup_years = 1 # set the minimum #years that a system will run in startup mode

# startup_years is the time that the system is in startup mode.
# At the end of the startup time the wells should be at their maximum capacity to run in a steady state mode.

years = 2        # set the minimum #years that a system will run
ppy = 365/perlen
sl = int(round (ppy * startup_years/perlen, 0)) # startup length                                                       
rl = int(round (365 * years / perlen, 0)) # run length

startwinter = 1 # simulation starts in winter (1) or summer (0). for distribution of flows of wells and surface temperature
        
nper = 1      # Number of Modflow stress periods per coupled simulation period
nstp = 1      # Number of Modflow time steps per stress period
timprs = np.array([perlen])
nprs = len(timprs)
savefiles = 0 # save head and temperature files (1) or not (0)
finaltimestep = 365/perlen*years-1
 
Tmin = 5        
T_room = 20                 # Room temperature in building, used for return temperature groundwater after heat exchanger
T_loss_building = 2         # degrees that the return temperature deviates from room temperature. The higher, the bigger the losses in the HVAC system. 

'''grid settings '''
dmin =  2             # smallest cel size at well [m]
dz =   5           # vertical gridlayer thickness [m]  important to syncronize with layer thicknesses in csv file!

dmin_bound = 200      # total distance from well with 'dmin' sized cells [m]
dmax = 200            # largest cell size at model boundary [m]
aroundAll = 1500      # normal=1500 [m] size of grid around well.
nstep = 20            # minimum number of steps that gridfunctions must add  to acquire the dmax size requirement [-]

AXI = 1              # axial symmetric grid or not. | 1=AXI, 0 = 3D |
LIN = 0              # Linear (1) or logarithmic (0) cell sizes buildup around the center of the modelgrid 
ICBUND_Up = -1       # TOP:: -1 = boundary has constant temperature (NO Flow, YES Temperature boundary), 1 = NO Flow, No Tempertature boundary
ICBUND_Down = -1     # BOTTOM: -1 = boundary has constant temperature (NO Flow, YES Temperature boundary), 1 = NO Flow, No Tempertature boundary
OutsideAirBound = 0  # 1 = ON, 0=OFF, if ON: the temperature boundary at model top is adjusted following outside air temperature 
steady = True        # Switch for steady-state/transient Modflow simulation
grid_extents = None   # example:[-300,300,-300,300] #set model boundaries [m]


#%% [B] Detailed MODFLOW/geohydr. inputs. Default conductivities and porosity - assumed constant for  model, can be handled as arrays like temperature/head

''' basic parameter settings'''
al = 0.5                                 # dispersivity [m]                                 
trpv = 0.005                             # trans dispersivity [m]                           
trpt = 0.05                              # transv. dispersivity [m]

rho_s = 2640.                            # density solids [kg/m3]
rho_f = 1000.                            # density fluids [kg/m3]
rho_b = rho_s * (1-PEFF) + rho_f * PEFF  # bulk density [kg/m3]

'''thermal properties'''
Cp_s = 710                                               # specific. heatcapacity [J/kg/K]
Cp_f = 4183.                                             # specific. heatcapacity [J/kg/K]
Cw = Cp_f * rho_f                                        # volumetric heatcapacity [J/m3/K]
Cs = Cp_s * rho_s                                        # volumetric heatcapacity [J/m3/K]
Caq =  Cs * (1-PEFF) + Cw * PEFF                         # volumetric aqufier heatcapacity
kT_s = 2.                                                # thermal conductivity sandy solids [W/m/K]
kT_clay = 1.7                                            # thermal conductivity clayy solids [W/m/K]
kT_f = 0.58                                              # thermal conductivity water [W/m/K]
kT_aq = kT_s * (1-PEFF) + kT_f * PEFF                    # thermal conductivity aquifer bulk
kT_aqt = kT_clay * (1 - PEFF) + kT_f * PEFF              # thermal conductivity aquitard bulk
Tdif_aq = kT_aq / (PEFF * rho_f * Cp_f) * 24 * 3600      # thermal diffusivity, deals with conduction in dsp package
Tdif_aqt = kT_aqt / (PEFF * rho_f * Cp_f) * 24 * 3600 
Kdist=Cp_s / (rho_f * Cp_f)                              # thermal distribution coeff, deals with thermal retardation in rct package

'''Seawat parameters'''
crho_ref = T_amb
visc_ref = 0.00002394*10**(248.37/(crho_ref+133.15))  # reference viscosity (VOSS 1982)

Tbh  = T_room + T_loss_building          # return temperature from HVAC
Tbc  = T_room - T_loss_building          # return temperature cooling from HVAC


for i in range(len(well_obj_list)):
    if well_obj_list[i].T_inj > Tmax:
        Tmax = well_obj_list[i].T_inj
    if well_obj_list[i].T_inj < Tmin:
        Tmin = well_obj_list[i].T_inj
denseref = 1000 - (((Tmin - 4)**2) / 207)      # reference density
densemin = 1000 - (((Tmax - 4)**2) / 207)      # minimum density (maximum temperature)
drhodT = (denseref - densemin) / (Tmin - Tmax) # linear relationship between T and p


'''Model flow input'''
flowtype = 2     # choise between 0 and 4, see readme
gwflow_x = -0.   # Groundwater flowvelocity in x direction [m/y]
gwflow_y = -0.   # Groundwater flowvelocity in y direction [m/y]

densflow = 0     # turn densityflow off or on by switching the VDF packages off or on  
WaterDensON = 0  # Calculation of parameters standard =1

'''correction factors determine the recharge of the wells by multiplying the needed flow to meet the demand by the factor'''
corr_w = 1.5     # correction factor hot well
corr_c = 1.0      # correction factor cold well

Thmin = 15          #Thmin is the temperature difference in the hot well between the injection temperature and the cutoff temperature

for i in well_obj_list:       # Update each active Python well object with the temperature and head at its grid location
    if i.type == 'warm':
        cutofftemp_h = i.T_inj-Thmin
        Tmax =  i.T_inj  # To calculate an approximation of the Temperature|density relation, the minimum(Tmax) and maximum density (4 C)
        # please note: seawat uses a linear approximation of this relation
cutofftemp_c = T_amb   # Temperature when the heating turns off, thus no more discharge from the hot well
 
Tmin = 5        
T_room = 20                 # Room temperature in building, used for return temperature groundwater after heat exchanger
T_loss_building = 2         # degrees that the return temperature deviates from room temperature. The higher, the bigger the losses in the HVAC system. 

'''grid settings '''
dmin =  1             # smallest cel size at well [m]
dz =   2              # vertical gridlayer thickness [m]  important to syncronize with layer thicknesses in csv file!

dmin_bound = 200      # total distance from well with 'dmin' sized cells [m]
dmax = 200            # largest cell size at model boundary [m]
aroundAll = 1000      # normal=1500 [m] size of grid around well.
nstep = 20            # minimum number of steps that gridfunctions must add  to acquire the dmax size requirement [-]

AXI = 1              # axial symmetric grid or not. | 1=AXI, 0 = 3D |
LIN = 0              # Linear (1) or logarithmic (0) cell sizes buildup around the center of the modelgrid 
ICBUND_Up = -1       # TOP:: -1 = boundary has constant temperature (NO Flow, YES Temperature boundary), 1 = NO Flow, No Tempertature boundary
ICBUND_Down = -1     # BOTTOM: -1 = boundary has constant temperature (NO Flow, YES Temperature boundary), 1 = NO Flow, No Tempertature boundary
OutsideAirBound = 0  # 1 = ON, 0=OFF, if ON: the temperature boundary at model top is adjusted following outside air temperature 
steady = True        # Switch for steady-state/transient Modflow simulation

grid_extents = None   # example:[-300,300,-300,300] #set model boundaries [m]
