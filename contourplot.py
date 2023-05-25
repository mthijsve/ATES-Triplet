import numpy as np
import os
import flopy 
import flopy.modflow as mf
import flopy.mt3d as mt3
import flopy.seawat as swt
import flopy.utils.binaryfile as bf

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


'''initialize contourlevels, plot boundaries, wells, layers, etc.'''
SY=0                                                                            # start year of simulation if applicable
Hlevel = np.linspace(-2,2, 100)                                                 # Head levels
Tlevel = np.linspace(Tmin-1, Tmax+1, 100)                                       # Temperature levels
T2level_w= [T_amb-0.5, T_amb+0.5, 25]                                           # levels for temperature contours
colors = ['r','b','k', 'g', 'c', 'y','grey', 'violet', 'dodgerblue','firebrick','coral','yellow', 'lightgreen', 'orange','firebrick','cyan']
X,Z=grid_obj.sumdelr, grid_obj.botm[:,0,0]

Z1 = int(np.around(np.average(well_obj_list[0].L)))
for i in range(len(grid_obj.sumdelr)):
    if grid_obj.sumdelr[i] < 150:                                               # right boundary of contour plots in m from well loaction
        a = i
    ResCol= int(a+1)

days = np.zeros(rl)
for i in range(rl):
    days[i] = perlen*(i+1)

OB = int(dmin_bound/dmin)*1.5                                                   # OutsideBox --> plot window / selection of UCN files
Wellb = np.zeros((len(well_obj_list),4))                                        # Define wells
for x in range(len(well_obj_list)):
    Wellb[x, 0]= well_obj_list[x].C - OB                                        # Left
    Wellb[x, 1]= well_obj_list[x].C + OB                                        # Right
    Wellb[x, 2]= well_obj_list[x].R - OB                                        # Down
    Wellb[x, 3]= well_obj_list[x].R + OB                                        # Up
Wellbmin =  Wellb.min(axis=0)
Wellbmax =  Wellb.max(axis=0)
Left  = int(Wellbmin[0])
Right = int(Wellbmax[1])
Down  = int(Wellbmin[2])
Up    = int(Wellbmax[3])
xw = [0,1]                                                                      # Define filter locations 
welltops = []                                                                   # creating a undefined list
wellbots = []
for i in range(nW):
    welltops.append(well_obj_list[i].ztop * np.ones(len(xw)))                   
    wellbots.append((well_obj_list[i].ztop - well_obj_list[i].ScreenLength) * np.ones(len(xw)))

n_aqt = 0                                                                       # Define aquitards
n_aq = 0

xa = np.ones(len(X[0:ResCol]))
aqttops = []
aqtbots = []
aq_tops = []
aq_bots = []
aqtform_l = []
aqform_l = []

for i in range(nform):                                                          # Assigning which layers are aquitards and which aquifers
    if form_obj_list[i].type == 'aquitard':
        n_aqt += 1
        aqtform_l.append(i)
    else:
        n_aq += 1
        aqform_l.append(i)
for i in range(n_aqt):                                                          # Assigning the thicknesses of the aquifers and between which depth(z) they range
    j = aqtform_l[i]
    aqttops.append(form_obj_list[j].zTop * xa)
    aqtbots.append(form_obj_list[j].zBot * xa)
for i in range(n_aq):
    j = aqform_l[i]
    aq_tops.append(form_obj_list[j].zTop * xa)
    aq_bots.append(form_obj_list[j].zBot * xa)

plots = np.linspace(0,rl-1,21).astype(int)
print(plots)                                         # plot the times indicated in this array
for i in range(len(plots)):

    t=int(plots[i])
    year = SY + int(((t+1) * perlen) / 365)
    month = int((((t+1) * perlen) - (year * 365)) / 30 +1  )
    textstr = 'YY-MM= ' + str(int( year)) +'-'+str(month)

    h_obj = bf.HeadFile(os.path.join(dirs, name+str(t)+'.hds'))
    head = h_obj.get_data(totim=perlen)
    h_obj.close() 
    t_obj = bf.UcnFile(os.path.join(dirs, name+str(t)+'S1'+'.UCN'))
    temp1 = t_obj.get_data(totim=perlen)
    t_obj.close()

    fig, (ax2)  = plt.subplots(nrows=1, sharex=True,figsize=(7,4),layout='tight') 
    im2 = ax2.contour(X[0:ResCol], Z, temp1[:,0,0:ResCol],
                        colors=('black'), linestyles=('dotted'), linewidths=1,levels = T2level_w)
    ax2.clabel(im2, im2.levels,fmt = '%1.1f')
    im2 = ax2.contourf(X[0:ResCol], Z, temp1[:,0,0:ResCol],cmap='seismic', levels = Tlevel) 
    fig.colorbar(im2, ax=ax2, format = '%1.0f')
    for i in range(nW):
        ax2.plot(xw, welltops[i], color = 'k') 
        ax2.plot(xw, wellbots[i], color = 'k')
        ax2.fill_between(xw, welltops[i], wellbots[i], color = 'k')
    for i in range(n_aqt):
        ax2.plot(X[0:ResCol], aqttops[i], color = 'grey', alpha = 0.6)
        ax2.plot(X[0:ResCol], aqtbots[i], color = 'grey', alpha = 0.6)
        ax2.fill_between(X[0:ResCol], aqttops[i], aqtbots[i], color = 'grey', alpha = 0.6)
    ax2.set_title('Temperature [C] of '+well_obj_list[0].type+'-well')
    ax2.set_ylabel('aquifer depth [m]')
    plt.xlabel('distance from well location [m]')
    props = dict(boxstyle='round', facecolor='w', alpha=0.5)
    ax2.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=10, verticalalignment='top', bbox=props)

    direct = os.path.join('Images',str(name)+str(t)+'.png')
    fig.savefig(direct, dpi=300) 