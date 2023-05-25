'''
flow calculations for PySeawaTriplet
adapted from grid_functions and agent_functions from PySeawATES
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def calc_demand(well_obj_list, Qyh, Qyc, run_length = 12, perlen=30, flowtype=1, demand=None, climate = None,years=1, startwinter = 1 ):
    
    '''  Calculate the average daily flow of a well object, as required by
    the Modflow WEL package.
    :param T_period: pandas series of daily temperatures over which the flows
                     should be computed
    flowtype 0: block function heat/cold (J), charge well next season
    flowtype 1: block function energy needed, charge constant throughtout year
    flowtype 2: sine function heat/cold (J), charge well next season with sine function
    flowtype 3: sine function heat/cold, charge well according to sine function through whole year based on
    ... (source for yearly charge pattern needed)
    flowtype 4: demand provided (J)
'''

    PerPerYear= int(round(365/perlen, 0))
    flowwarm = np.zeros(run_length)
    flowcold = np.zeros(run_length)
    chargewarm = np.zeros(run_length)
    chargecold = np.zeros(run_length)
    
    for i in well_obj_list:
        i.flow = np.zeros(run_length)
        i.charge = np.zeros(run_length)

    if flowtype == 0: #block function heat/cold (J), charge well next season

        PerPerSeason = PerPerYear/2
        start_winter = np.ones(years)
        start_summer = np.ones(years)
        
        if startwinter == 1:
            
            for a in range(years):
                start_winter[a] = a * PerPerYear
                start_summer[a] = a * PerPerYear + PerPerYear/2
    
            for i in well_obj_list:
                for j in range(years):
                    k=int(start_winter[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = -Qyh/PerPerSeason
                            i.charge[k+l] = 0

                        if i.type =='cold':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyc/PerPerSeason
                            
                    k=int(start_summer[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyh/PerPerSeason
                            
                        if i.type =='cold':
                            i.flow[k+l] = -Qyc/PerPerSeason
                            i.charge[k+l] = 0
                            
                if i.type =='warm':
                    flowwarm = i.flow
                    chargewarm = i.charge
                        
                if i.type =='cold':
                    flowcold = i.flow
                    chargecold = i.charge
                    
        if startwinter == 0:
            
            for a in range(years):
                start_winter[a] = a * PerPerYear
                start_summer[a] = a * PerPerYear + PerPerYear/2
    
            for i in well_obj_list:
                for j in range(years):
                    k=int(start_winter[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = -Qyh/PerPerSeason
                            i.charge[k+l] = 0

                        if i.type =='cold':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyc/PerPerSeason
                            
                    k=int(start_summer[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyh/PerPerSeason
                            
                        if i.type =='cold':
                            i.flow[k+l] = -Qyc/PerPerSeason
                            i.charge[k+l] = 0            

                if i.type =='warm':
                    flowwarm = i.flow
                    chargewarm = i.charge
                        
                if i.type =='cold':
                    flowcold = i.flow
                    chargecold = i.charge
                    
        # plt.figure()
        # plt.plot(-flowwarm,label='heating demand',color='#EC6842')
        # plt.plot(-flowcold,label='cooling demand',color='#00B8C8')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()
        
        # plt.figure()
        # plt.plot(-chargewarm,label='heating charge',color='#EC6842')
        # plt.plot(-chargecold,label='cooling charge',color='#00B8C8')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show() 
        
    if flowtype == 1: #block function energy needed, charge block function throughout the year
        
        PerPerSeason = PerPerYear/2
        start_winter = np.ones(years)
        start_summer = np.ones(years)
        
        if startwinter == 1:
            
            for a in range(years):
                start_winter[a] = a * PerPerYear
                start_summer[a] = a * PerPerYear + PerPerYear/2
    
            for i in well_obj_list:
                for j in range(years):
                    k=int(start_winter[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = -Qyh/PerPerSeason
                            i.charge[k+l] = Qyh/PerPerSeason/2

                        if i.type =='cold':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyc/PerPerSeason/2
                            
                    k=int(start_summer[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyh/PerPerSeason/2
                            
                        if i.type =='cold':
                            i.flow[k+l] = -Qyc/PerPerSeason
                            i.charge[k+l] = Qyc/PerPerSeason/2
                            
                if i.type =='warm':
                    flowwarm = i.flow
                    chargewarm = i.charge
                        
                if i.type =='cold':
                    flowcold = i.flow
                    chargecold = i.charge
                    
        if startwinter == 0:
            
            for a in range(years):
                start_winter[a] = a * PerPerYear
                start_summer[a] = a * PerPerYear + PerPerYear/2
    
            for i in well_obj_list:
                for j in range(years):
                    k=int(start_winter[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = -Qyh/PerPerSeason
                            i.charge[k+l] = Qyh/PerPerSeason/2

                        if i.type =='cold':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyc/PerPerSeason/2
                            
                    k=int(start_summer[j])
                    for l in range(int(PerPerSeason)):
                        
                        if i.type =='warm':
                            i.flow[k+l] = 0
                            i.charge[k+l] = Qyh/PerPerSeason/2
                            
                        if i.type =='cold':
                            i.flow[k+l] = -Qyc/PerPerSeason
                            i.charge[k+l] = Qyc/PerPerSeason/2  
                            
                if i.type =='warm':
                    flowwarm = i.flow
                    chargewarm = i.charge
                        
                if i.type =='cold':
                    flowcold = i.flow
                    chargecold = i.charge
                    
        # plt.figure()
        # plt.plot(-flowwarm,label='heating demand',color='#EC6842')
        # plt.plot(-flowcold,label='cooling demand',color='#00B8C8')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()
        
        # plt.figure()
        # plt.plot(-chargewarm,label='heating charge',color='#EC6842')
        # plt.plot(-chargecold,label='cooling charge',color='#00B8C8')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()

                    
    if flowtype == 2: #sine function heat/cold (J), charge well next season with sine function

        for i in well_obj_list:
            if i.type == 'warm':
                for j in range(run_length):
                    if startwinter == 1:
                        i.flow[j] = -np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyh/PerPerYear
                    else:
                        i.flow[j] = np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyh/PerPerYear   
                        
                i.flow = i.flow.clip(max=0)
                flowwarm = i.flow
                
            if i.type == 'cold':
                for j in range(run_length):
                    if startwinter == 1:
                        i.flow[j] = np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyc/PerPerYear
                    else:
                        i.flow[j] = -np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyc/PerPerYear
                        
                i.flow = i.flow.clip(max=0)
                flowcold = i.flow
        
        # plt.figure()
        # plt.plot(-flowwarm,label='heating demand',color='#EC6842')
        # plt.plot(-flowcold,label='cooling demand',color='#00B8C8')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()
        
        for i in well_obj_list:
            if i.type == 'warm':
                for j in range(run_length):
                    if startwinter == 1:
                        i.charge[j] = -np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyh/PerPerYear
                    else:
                        i.charge[j] = np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyh/PerPerYear
                        
                i.charge = i.charge.clip(min=0)
                chargewarm = i.charge
                
            if i.type == 'cold':
                for j in range(run_length):
                    if startwinter == 1:
                        i.charge[j] = np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyc/PerPerYear
                    else:
                        i.charge[j] = -np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyc/PerPerYear
                        
                i.charge = i.charge.clip(min=0)
                chargecold = i.charge
                
        # plt.figure()
        # plt.plot(-chargewarm,label='heating charge',color='#EC6842')
        # plt.plot(-chargecold,label='cooling charge',color='#00B8C8')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()            
                
    if flowtype == 3: #sine function heat/cold, charge well according to sine function through whole year
        for i in well_obj_list:
            if i.type == 'warm':
                for j in range(run_length):
                    if startwinter == 1:
                        i.flow[j] = -np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyh/PerPerYear
                    else:
                        i.flow[j] = np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyh/PerPerYear   
                        
                i.flow = i.flow.clip(max=0)
                flowwarm = i.flow
                
            if i.type == 'cold':
                for j in range(run_length):
                    if startwinter == 1:
                        i.flow[j] = np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyc/PerPerYear
                    else:
                        i.flow[j] = -np.sin(2*np.pi*j/PerPerYear)*np.pi*Qyc/PerPerYear
                        
                i.flow = i.flow.clip(max=0)
                flowcold = i.flow
        
        # plt.figure()
        # plt.plot(-flowwarm,label='heating',color='#EC6842')
        # plt.plot(-flowcold,label='cooling',color='#00B8C8')
        # plt.title('Yearly heating and cooling demand')
        # plt.ylabel('J / day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()
        
        for i in well_obj_list:
            if i.type == 'warm':
                for j in range(run_length):
                    if startwinter == 1:
                        i.charge[j] = -np.sin(2*np.pi*j/PerPerYear)*Qyh/PerPerYear + Qyh/PerPerYear
                        
                    else:             
                        i.charge[j] = np.sin(2*np.pi*j/PerPerYear)*Qyh/PerPerYear + Qyh/PerPerYear
                chargewarm = i.charge
            if i.type == 'cold':
                for j in range(run_length):
                    if startwinter == 1:
                        i.charge[j] = np.sin(2*np.pi*j/PerPerYear)*Qyc/PerPerYear + Qyc/PerPerYear
                        
                    else:
                        i.charge[j] = -np.sin(2*np.pi*j/PerPerYear)*Qyc/PerPerYear + Qyc/PerPerYear
                chargecold = i.charge
                
        # plt.figure()
        # plt.plot(chargewarm,label='hot charging',color='#EC6842')
        # plt.plot(chargecold,label='cold charging',color='#00B8C8')
        # plt.title('Yearly charging of the wells if there are no losses')
        # plt.ylabel('J / day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()                      
                
    if flowtype == 4:
        '''supply total demand in "demand.csv'''
        demand = pd.read_csv(demand,skiprow = [0],delimiter=';')
        for i in well_obj_list:
            if i.type == 'warm':
                i.flow = int(demand.loc[i,0])
                i.charge = int(demand.loc[i,2])
                
            flowwarm = i.flow
            chargewarm = i.charge
                
            if i.type == 'cold':
                i.flow = int(demand.loc[i,1])
                i.charge = int(demand.loc[i,3])
                
            flowcold = i.flow
            chargecold = i.charge
            
        # plt.figure()
        # plt.plot(-flowwarm,label='heating demand')
        # plt.plot(-flowcold,label='cooling demand')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()
        
                        
        # plt.figure()
        # plt.plot(-chargewarm,label='heating charge')
        # plt.plot(-chargecold,label='cooling charge')
        # plt.ylabel('J/day')
        # plt.xlabel('day')
        # plt.legend()
        # plt.show()          
