from functions_triplet import *
from flow_function_Triplet import *

import pandas as pd
import numpy as np
import flopy
import flopy.modflow as mf
import flopy.mt3d as mt3
import flopy.seawat as swt
import flopy.utils.binaryfile as bf
import os
from pathlib import Path
import time

def Modelrun(corr_w,Qyh, Qyc, injectionT, Thmin):
    print('start_modelrun with corr_w = ', corr_w,'. Qyh=', Qyh, 'Qyc=', Qyc, 'injectionT=', injectionT, 'Thmin=', Thmin)
    corr_ws = corr_w

    '''Define input files'''
    swtexe_name = 'swt_v4x64.exe' 
    wellfile = 'wells.csv'  
    subsurffile = 'sub-surface.csv'     
    monitoringfile = 'monitoring.csv'

    p = str(Qyh).replace('.', '')
    x = str(Qyc).replace('.', '')
    y = str(injectionT).replace('.', '')
    z = str(Thmin).replace('.', '')
    name = 'sensitivity_Qh' + p + '_Qc' + x + '_injectionT' + y + '_Thmin' + z

    well_obj_list = createobjfromCSV(PyWell,wellfile) #creates well object from  data in csv file
    form_obj_list = createobjfromCSV(PySystem,subsurffile)
    mon_obj_list = createobjfromCSV(PySystem,monitoringfile) # creates list of monitoring points in your grid where you can track the temperature

    nmon = len(mon_obj_list)
    nform = len(form_obj_list)
    nW=len(well_obj_list)         # numer of wells
    T_amb = form_obj_list[0].s1   # Ambient Temperature of the subsurface
    PEFF = form_obj_list[0].por   # porosity [-]
    WD = os.getcwd()
    dirs = Path(WD + '/output/' + name)                                                      
    if os.path.exists(dirs)==False:
        os.makedirs(dirs)

    '''main time setting parameters'''
    perlen = 1    # (DAYS)
    startup_years = 3 # set the minimum #years that a system will run in startup mode

    # startup_years is the time that the system is in startup mode.
    # At the end of the startup time the wells should be at their maximum capacity to run in a steady state mode.

    years = 8        # set the minimum #years that a system will run
    ppy = 365/perlen
    sl = int(round (ppy * startup_years/perlen, 0)) # startup length                                                       
    rl = int(round (365 * years / perlen, 0)) # run length

    startwinter = 0 # simulation starts in winter (1) or summer (0). for distribution of flows of wells and surface temperature
            
    nper = 1      # Number of Modflow stress periods per coupled simulation period
    nstp = 1      # Number of Modflow time steps per stress period
    timprs = np.array([perlen])
    nprs = len(timprs)
    savefiles = 0 # save head and temperature files (1) or not (0)
    finaltimestep = 365/perlen*years-1
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

    '''room parameters'''
    Tmin = 5        
    T_room = 20                 # Room temperature in building, used for return temperature groundwater after heat exchanger
    T_loss_building = 2         # degrees that the return temperature deviates from room temperature. The higher, the bigger the losses in the HVAC system. 
    Tbh  = T_room + T_loss_building          # return temperature from HVAC
    Tbc  = T_room - T_loss_building          # return temperature cooling from HVAC

    for i in well_obj_list:       # Update each active Python well object with the temperature and head at its grid location
        if i.type == 'warm':
            i.T_inj = injectionT
            cutofftemp_h = i.T_inj*Thmin
            Tmax =  i.T_inj  # To calculate an approximation of the Temperature|density relation, the minimum(Tmax) and maximum density (4 C)
            # please note: seawat uses a linear approximation of this relation

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
    corr_c = 1.0      # correction factor cold well
    corr_cs = 1.0     # correction factor cold well startup
    cutofftemp_c = T_amb   # Temperature when the heating turns off, thus no more discharge from the hot well
    
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


    '''run period info'''
    Run_info = pd.DataFrame()                 # create run info dataframe
    Run_info.loc[0,'Period'] = 0
    Run_info.loc[0,'Day'] = 0
    for i in range(rl):    
        Run_info.loc[i,'Period'] = int(i*1)
        Run_info.loc[i,'Day'] = (i+1)*(perlen)
    Run_output = Run_info.copy()              # output dataframe

    zbot = form_obj_list[nform-1].zBot   # bottom of the modelgrid [m]
    ztop = form_obj_list[0].zTop         # top of the modelgrid [m]

    if AXI == 1:
        grid_obj = PyGridAxi()    # for axi-sym, geo-propterties are automatically adjusted see grid_functions for details
        if LIN == 0:      # Axisymmetrical grid with logarithmic increasing cell sizes from well 
            grid_obj.make_grid_axi_Log(well_obj_list, dmin=dmin, dmax=dmax, dz=dz,dmin_bound=dmin_bound,
                            nstep=int(nstep), ztop=ztop, zbot=zbot, aroundAll=aroundAll, grid_extents=grid_extents, ICBUND_Up=ICBUND_Up, ICBUND_Down=ICBUND_Down)
        else:                                                                       # Axisymmetrical grid with a zone with constant cell size at well, logarithmic increasing cell sizes from LIN boundary
            grid_obj.make_grid_axi_Lin(well_obj_list, dmin=dmin, dmax=dmax, dz=dz,dmin_bound=dmin_bound,
                            nstep=int(nstep), ztop=ztop, zbot=zbot, aroundAll=aroundAll, grid_extents=grid_extents, ICBUND_Up=ICBUND_Up, ICBUND_Down=ICBUND_Down)
    else:
        grid_obj = PyGrid()  
        if LIN ==0:                                                                 # 3D grid with logarithmic increasing cell sizes from well(s) 
            grid_obj.make_grid_Log(well_obj_list, dmin=dmin, dmax=dmax, dz=dz,dmin_bound=dmin_bound,
                                nstep=int(nstep), ztop=ztop, zbot=zbot, aroundAll=aroundAll, grid_extents=grid_extents, ICBUND_Up=ICBUND_Up, ICBUND_Down=ICBUND_Down)
        else:                                                                       # 3D grid with a constant cell size basis around the well(s). logarithmic increasing cell sizes from set boundary.
            grid_obj.make_grid_Lin(well_obj_list, dmin=dmin, dmax=dmax, dz=dz,dmin_bound=dmin_bound,
                                nstep=int(nstep), ztop=ztop, zbot=zbot, aroundAll=aroundAll, grid_extents=grid_extents, ICBUND_Up=ICBUND_Up, ICBUND_Down=ICBUND_Down)

    '''set geo properties for each gridcell'''
    set_geo_prop(grid_obj, form_obj_list, gwflow_x=gwflow_x, gwflow_y=gwflow_y, dz=dz, AqL=2, Tdif_aqt=Tdif_aqt, Tdif_aq=Tdif_aq, rho_b=rho_b,rho_f=rho_f,Cp_f=Cp_f,Kdist=Kdist)
    if OutsideAirBound == 1:
        set_Tboundary(grid_obj,  perlen=perlen, run_length=rl, Tmin=5, Tmax=20, startwinter = startwinter)  # sets varying temperature at surface level. 

    '''monitoring Layer, Row, Column number in grid'''
    mon_LRC_list = init_monitoring(grid_obj, mon_obj_list, dz, nmon, AXI)
    RES = np.zeros((rl,len(mon_LRC_list)))                                          # results file for monitoring points

    '''initialize well properties''' 
    for i in well_obj_list: 
        i.calc_LRC(grid_obj)                                                        # Locate each well in the simulation grid
        i.S_inj = 0
        
    itype = mt3.Mt3dSsm.itype_dict()
    well_LRCQ_list = {}
    well_LRCQ_list[0] = [[0, 0, 0, 0]]
    ssm_data = {}
    ssm_data[0] = [[0, 0, 0, 0, itype['WEL']]]
    laytyp = np.zeros([grid_obj.nlay])                                              # All confined layers

    '''calculate flows '''
    calc_demand(well_obj_list, Qyh , Qyc, perlen=perlen, flowtype=flowtype, run_length=rl, years=years, startwinter=startwinter) # calculates flows during simulation
    ''' Create SEAWAT model '''                                                     # set all the fixed conditions that don't change while running etc.
    mswtf = swt.Seawat(name, 'nam_swt',exe_name=swtexe_name,
                    model_ws=dirs)
    #to change to hours, add itmuni here, see documentation for options
    discret = mf.ModflowDis(mswtf, nrow=int(grid_obj.nrow), ncol=int(grid_obj.ncol), nlay=grid_obj.nlay,
                    delr=grid_obj.delr, delc=grid_obj.delc, laycbd=0., top=grid_obj.top, 
                    botm=grid_obj.botm, nper=nper, perlen=perlen, nstp=nstp, steady=steady)

    #%% [D] model itteration

    start_time = time.time()
    T_a_h = T_amb
    T_a_c = T_amb
    T_a_b = T_amb 

    for period in range(rl):

        if OutsideAirBound == 1:
            grid_obj.temp[0,:,:]= grid_obj.SurfaceT[period]# sets temperature at surface level according to time of year
        temp_QH=0
        temp_QC=0
        if well_obj_list:  # Create well and temperature lists following Modflow/MT3DMS format. Each timestep, flows and infiltration temp are assigned to the wells  
            for i in well_obj_list:
                if i.type == 'warm':   
                    if period < ppy*startup_years:           
                        if T_a_h < cutofftemp_h:
                            i.Q = corr_ws*pumpingrate(i.charge[period], Tbh, i.T_inj,Cw)/perlen
                            Qh = 0
                        else:
                            i.Q = (pumpingrate(i.flow[period] ,T_a_h ,Tbh ,Cw) + pumpingrate(i.charge[period], Tbh, i.T_inj,Cw)*corr_ws)/perlen
                            Qh = abs(pumpingrate(i.flow[period] ,T_a_h ,Tbh ,Cw))
                        temp_QH = i.Q
                    else:
                        if T_a_h < cutofftemp_h:
                            i.Q = corr_w*pumpingrate(i.charge[period], T_a_h, Tbh,Cw)/perlen
                            Qh = 0
                        else:
                            i.Q = (pumpingrate(i.flow[period] ,T_a_h ,Tbh ,Cw) + pumpingrate(i.charge[period], Tbh, i.T_inj,Cw)*corr_w)/perlen
                            Qh = abs(pumpingrate(i.flow[period] ,T_a_h ,Tbh ,Cw))
                        temp_QH = i.Q
                if i.type == 'cold':
                    if period < ppy*startup_years:
                        if T_a_c <cutofftemp_c:
                            i.Q = (pumpingrate(i.flow[period],T_a_c, Tbc,Cw) + pumpingrate(i.charge[period], Tbc,i.T_inj , Cw)*corr_cs)/perlen
                            
                            Qc = abs(pumpingrate(i.flow[period],T_a_c, Tbc,Cw))
                        else:
                            Qc = 0
                            i.Q = (pumpingrate(i.charge[period], Tbc, i.T_inj, Cw)*corr_cs)/perlen
                        temp_QC = i.Q
                    
                    else:
                        if T_a_c <cutofftemp_c:
                            i.Q = (pumpingrate(i.flow[period],T_a_c, Tbc,Cw) + pumpingrate(i.charge[period], Tbc, i.T_inj, Cw)*corr_c)/perlen
                            Qc = abs(pumpingrate(i.flow[period],T_a_c, Tbc,Cw))
                        else:
                            Qc = 0
                            i.Q = (pumpingrate(i.charge[period], Tbc, i.T_inj, Cw)*corr_c)/perlen
                            
                        temp_QC = i.Q
                    
                if i.type == 'buffer':
                    i.Q = -temp_QH -temp_QC
                
                    if Qc+Qh >0:
                        i.T_inj = (Qh*Tbh + Qc*Tbc)/(Qc+Qh)

        
        well_LRCQ_list = create_LRCQ_list(well_obj_list, grid_obj)                  # Set the settings of the wells for that timestep
        ssm_data = create_conc_list(well_obj_list, attrib='T_inj') 
        
        '''Initialize MODFLOW Packages'''
        bas = mf.ModflowBas(mswtf, ibound=grid_obj.IBOUND, strt=grid_obj.head)         # Basemodel Modflow
        wel = mf.ModflowWel(mswtf, stress_period_data = well_LRCQ_list)                  # Put in Wells
        words = ['head','drawdown','budget', 'phead', 'pbudget']                
        save_head_every = 1
        lpf = mf.ModflowLpf(mswtf, hk=grid_obj.HK, vka=grid_obj.VK, ss=grid_obj.ss, sy=0.15, laytyp=laytyp, layavg=0.,ipakcb=53) 

        oc = mf.ModflowOc(mswtf)                                                       # Output control package class --> moved (p3.7 iso p3.6)                   
        pcg = mf.ModflowPcg(mswtf, mxiter=200, iter1=200, npcond=1,                    # Preconditioned Conjugate-Gradient Package  --> solves the finite differences equations
                            hclose=0.001, rclose=0.001, relax=1.0, nbpol=0)
        
        '''Initialize MT3DMS packages''' 
        adv = mt3.Mt3dAdv(mswtf, mixelm=0, percel=0.8, nadvfd=1, #Particle based methods
                        nplane=0,mxpart=250000, dceps=1e-4, 
                        npl=5, nph=8, npmin=1, npmax=16)
        btn = mt3.Mt3dBtn(mswtf, cinact=-100., icbund=grid_obj.ICBUND, prsity=grid_obj.PEFF,
                        sconc=grid_obj.temp, #sconc2=grid_obj.salinity,
                        ifmtcn=-1, chkmas=False, nprobs=0, nprmas=1, dt0=0, ttsmult=1,
                        ttsmax=1., ncomp=1, mcomp=1, nprs=nprs, timprs=timprs, mxstrn=9999)
        dsp = mt3.Mt3dDsp(mswtf, al=al, trpt=trpt, trpv=trpv, multiDiff=True, dmcoef = grid_obj.Tdif) 
        rct = mt3.Mt3dRct(mswtf, isothm=1, ireact=0, igetsc=0, rhob=grid_obj.rho_b, sp1=grid_obj.Kdist)
        gcg = mt3.Mt3dGcg(mswtf, mxiter=50, iter1=50, isolve=1, cclose=1e-9, iprgcg=0)
        ssm = mt3.Mt3dSsm(mswtf, stress_period_data = ssm_data)
            
        '''Initialize SEAWAT packages''' #use this package when density/viscosity dependency is needed
        if densflow == 1:
            vdf = swt.SeawatVdf(mswtf, mtdnconc=-1, nsrhoeos=1, nswtcpl=1, iwtable=0, densemin=0, densemax=0, denseref=denseref,
                                mtrhospec=1,denseslp=drhodT,crhoref=crho_ref)
            vsc = swt.SeawatVsc(mswtf, mt3dmuflg=-1, viscmin=0.0, viscmax=0.0, viscref=visc_ref, #Viscref must be set to the reference viscosity at T=12 --> 
                                nsmueos=1, mtmutempspec=1, mutempopt=1, amucoeff=(2.39e-5, 10, 248.4, 133.2),  # temp is used to calc visc. according to Eq. 18 langevin et al 2008
                                mtmuspec=2, dmudc=1.923e-06, cmuref=0.0,                # solute influence on viscocity
                                invisc=-1, visc=-1, extension='vsc')
        
        mswtf.write_input()
        m = mswtf.run_model(silent=True)   #Run SEAWAT!      (silent=FALSE gives details for each timestep, silent=TRUE only gives end signal in the console)
        

        if savefiles == 1:
            '''Copy Modflow/MT3DMS output to new files so they wont be overwritten in next timestep.'''
            shutil.copyfile(os.path.join(dirs, name+'.hds'),
                            os.path.join(dirs, name+str(period)+'.hds'))
            shutil.copyfile(os.path.join(dirs, 'MT3D001.UCN'),
                            os.path.join(dirs, name+str(period)+'S1'+'.UCN'))
            '''Create head & concentrations file object and read head & concentrations arrays for next simulation period'''
            h_obj = bf.HeadFile(os.path.join(dirs, name+str(period)+'.hds'))
            grid_obj.head = h_obj.get_data(totim=perlen)
            t_obj = bf.UcnFile(os.path.join(dirs, name+str(period)+'S1'+'.UCN'))
            grid_obj.temp = t_obj.get_data(totim=perlen)                
        else:
            '''Create head & concentrations file object and read head & concentrations arrays for next simulation period'''
            h_obj = bf.HeadFile(os.path.join(dirs, name+'.hds'))
            grid_obj.head = h_obj.get_data(totim=perlen)
            t_obj = bf.UcnFile(os.path.join(dirs, 'MT3D001.UCN'))
            grid_obj.temp = t_obj.get_data(totim=perlen) 

        if well_obj_list:
            for i in well_obj_list:       # Update each active Python well object with the temperature and head at its grid location
                i.H_modflow = grid_obj.head[i.L[-1],i.R,i.C]
                i.T_modflow = np.average(grid_obj.temp[i.start_idx:i.stop_idx,i.R,i.C]) #the average of all the cells of the injection well! (start.idx and stop.idx)
                if i.type == 'warm':
                    T_a_h = i.T_modflow
                    
                if i.type == 'cold':
                    T_a_c = i.T_modflow
                    
                if i.type == 'buffer':
                    T_a_b = i.T_modflow
        '''Save temp monitoring pointsdata to results array'''
        for m in range(len(mon_LRC_list)):
            RES[period,m] = grid_obj.temp[int(mon_LRC_list[m,0]),int(mon_LRC_list[m,1]),int(mon_LRC_list[m,2])]    
        
        '''save the info the the Run_output file'''
        for j in range(len(well_obj_list)):
            if well_obj_list[j].Q <0:
                Run_output.loc[period,'W'+str(j)+'_Vin'] =  0
                Run_output.loc[period,'W'+str(j)+'_Vout'] =  well_obj_list[j].Q*perlen
                Ttemp = well_obj_list[j].T_modflow
            else:
                Run_output.loc[period,'W'+str(j)+'_Vin'] = well_obj_list[j].Q*perlen
                Run_output.loc[period,'W'+str(j)+'_Vout'] = 0  
                Ttemp = well_obj_list[j].T_modflow
                
            if well_obj_list[j].type == 'warm':
                Run_output.loc[period,'W'+str(j)+'_T_sys_in'] = well_obj_list[j].T_inj
                Run_output.loc[period,'W'+str(j)+'_T_mf_out'] = well_obj_list[j].T_modflow  
                Run_output.loc[period,'dnmh'] = Run_output.loc[period,'W'+str(j)+'_T_mf_out'] < cutofftemp_h            
                Run_output.loc[period,'Dens_water_in'] = denseref + drhodT * well_obj_list[j].T_inj
                Run_output.loc[period,'Dens_water_out'] = denseref + drhodT * well_obj_list[j].T_modflow
                if period < ppy*startup_years:
                    Run_output.loc[period,'Efficiency_h'] = corr_ws
                else:
                    Run_output.loc[period,'Efficiency_h'] = corr_w
            elif well_obj_list[j].type == 'cold':
                Run_output.loc[period,'W'+str(j)+'_T_sys_in'] = well_obj_list[j].T_inj
                Run_output.loc[period,'W'+str(j)+'_T_mf_out'] = well_obj_list[j].T_modflow
                Run_output.loc[period,'dnmc'] = Run_output.loc[period,'W'+str(j)+'_T_mf_out'] > cutofftemp_c 
                Run_output.loc[period,'Dens_water_in'] = denseref + drhodT * well_obj_list[j].T_inj
                Run_output.loc[period,'Dens_water_out'] = denseref + drhodT * well_obj_list[j].T_modflow
                if period < ppy*startup_years:
                    Run_output.loc[period,'Efficiency_h'] = corr_cs
                else:
                    Run_output.loc[period,'Efficiency_h'] = corr_c
            else:
                Run_output.loc[period,'W'+str(j)+'_T_sys_in'] = Ttemp
                Run_output.loc[period,'W'+str(j)+'_T_mf_out'] = well_obj_list[j].T_modflow
                Run_output.loc[period,'Dens_water_in'] = denseref + drhodT * well_obj_list[j].T_inj
                Run_output.loc[period,'Dens_water_out'] = denseref + drhodT * well_obj_list[j].T_modflow
        
        '''Density water influenced by Temperature'''                               # Update parameters based on rho_f per gridcell
        if WaterDensON ==1:
            grid_obj.rho_f = (1000 - ((grid_obj.temp-4)**2 / 207))                  # calculate updated water density for each cell       
            grid_obj.Kdist = Cp_s/(grid_obj.rho_f*Cp_f)                             # calculate updated Kdist for each cell (with updated rho_f
            for j in range (len(form_obj_list)):                                    # calculate updated Tdif for each layer (or Aqt or Aq)
                for k in range(int(form_obj_list[j].lbot - form_obj_list[j].ltop)):
                    if form_obj_list[j].type == 'aquitard':
                        grid_obj.Tdif[k,:,:] =  kT_aqt / (PEFF*grid_obj.rho_f[k,:,:]*Cp_f) * 24 * 3600
                    else:
                        grid_obj.Tdif[k,:,:] = kT_aq / (PEFF*grid_obj.rho_f[k,:,:]*Cp_f) * 24 * 3600 
    
        print (str(period)+' (of '+str(rl-1)+')') 
        h_obj.file.close()
        t_obj.file.close()

        elapsed_time = time.time() - start_time
        proxy = elapsed_time/60 * rl/(period+1) - elapsed_time/60
        print ('run time =', int(elapsed_time/60),' min')
        #print ('expected remaining run time =', int(proxy),' min')


    end_time = time.time()
    elapsed_time = (end_time- start_time) / 60
    print ('Simulation completed, Runtime= '+str(round(elapsed_time,1))+' min') 

    '''save run_output file'''
    Run_output.to_csv(os.path.join(dirs)+'/Run_output__'+name+'.csv')                
    mask2 = (Run_output['Period']> round(startup_years*ppy))
    #values dnm.. stands for demand not met (cold or hot)(startup period or not). This is used to calibrate the model. It should be 0 in the last year of startup and 0 in the years after
    Run_output['dnmh'].replace({False: 0, True: 1}, inplace=True)
    dnmh = Run_output['dnmh'].loc[mask2].sum()
    print('dmnm:',dnmh)
    return dnmh
