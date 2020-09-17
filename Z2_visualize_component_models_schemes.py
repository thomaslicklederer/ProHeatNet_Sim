import math
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from IPython import get_ipython
import os
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import gurobipy as gp
from gurobipy import GRB

import D_component_models as cm

################ delete old data #########################
get_ipython().magic('reset -f')
if os.name == 'nt':
    get_ipython().magic('cls')
else:
    get_ipython().magic('clear')
plt.close(fig='all')
##########################################################

############################ general parameters ############################

myfigsize = [2*20.0/2.54, 2*10.0/2.54]
myfigsize2= [2*10.0/2.54, 2*10.0/2.54]
myfontsize = 20
mythickness = 5
mylinestyles = [':', '-.', '--', '-']
mylinecolors = ['r','g','b','k']
mylabelpad   = 15
savepath = './figs/'
saveformat = 'png' # jpg, png, eps, svgd
show_pump       = True
show_valve      = True
show_pipe_hy    = True
show_pipe_th    = True
show_hx         = True     
mytolerance     = 10**(-2)
mpl.rcParams['axes.linewidth']=1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.major.width'] = 1.5

############################################################################

############################ define components #############################

myfluid = cm.fluid(rho=1000, cp=4200, mu=(1.0016*(10**(-3))))
mypump  = cm.pump(n_nom=4100, u_ref_1=1, dotV_ref_1=0, Deltap_ref_1=402.21,
    u_ref_2=1, dotV_ref_2=55.33, Deltap_ref_2=0)
myvalve = cm.controlvalve(K_vs=2.5)
mypipe  = cm.pipe(L=100,d_hy=0.022,epsilon=0.011,u_nom=1,zeta_instal=10,
    N_layers = 3, d_layers = [0.022, 0.024, 0.062, 0.262],
    lambda_layers = [395, 0.04, 2], h_ir=6700, h_or=100)
myhx    = cm.heatexchanger(dotV_nom=21.504, Deltap_nom=155, k_nom = 5270, A = 1.13)

############################################################################

############################ define test values ############################

u_test          =   [0.25,0.5,0.75,1]
kappa_test      =   [0.25,0.5,0.75,1]
dotV_vec_test   =   np.linspace(0,30,500)
dotV_vec_pith_test  =   dotV_vec_test[dotV_vec_test>=1]
mu_test         =   -1
dotV_PSM_const_test     = 18 
dotV_HTNW_const_test    = 14
dotV_pi_const_test      = 15
T_in_PSM_const_test     = 45+273.15
T_in_HTNW_const_test    = 65+273.15
T_soil_test     =   12+273.15
T_in_pi_test    =   np.linspace(20,90,300)+273.15
T_in_pi_const_test    =   65+273.15
T_in_hx_test   =   np.linspace(45,65,300)+273.15
dotV_hx_test   =   np.linspace(1,30,300)



############################################################################

############################ initialize result vectors ######################

Deltap_pump     =   np.zeros((len(dotV_vec_test),len(u_test)), dtype=float)
Deltap_valve    =   np.zeros((len(dotV_vec_test),len(kappa_test)), dtype=float)
Deltap_pipe     =   np.zeros(len(dotV_vec_test), dtype=float)
DeltaT_pipe_Tfix     =   np.zeros(len(dotV_vec_pith_test), dtype=float)
DeltaT_pipe_Vfix     =   np.zeros(len(T_in_pi_test), dtype=float)
T_out_pipe_Tfix      =   np.zeros(len(dotV_vec_pith_test), dtype=float)
T_out_pipe_Vfix      =   np.zeros(len(T_in_pi_test), dtype=float)
T_out_hxPSM_Tfix    = np.zeros((len(dotV_hx_test), len(dotV_hx_test)), dtype=float)
T_out_hxHTNW_Tfix   = np.zeros((len(dotV_hx_test), len(dotV_hx_test)), dtype=float)
T_out_hxPSM_Vfix    = np.zeros((len(T_in_hx_test), len(T_in_hx_test)), dtype=float)
T_out_hxHTNW_Vfix   = np.zeros((len(T_in_hx_test), len(T_in_hx_test)), dtype=float)
DeltaT_hxPSM_Tfix   = np.zeros((len(dotV_hx_test), len(dotV_hx_test)), dtype=float)
DeltaT_hxHTNW_Tfix  = np.zeros((len(dotV_hx_test), len(dotV_hx_test)), dtype=float)
DeltaT_hxPSM_Vfix   = np.zeros((len(T_in_hx_test), len(T_in_hx_test)), dtype=float)
DeltaT_hxHTNW_Vfix  = np.zeros((len(T_in_hx_test), len(T_in_hx_test)), dtype=float)

############################################################################

############################ calculate result arrays #######################

if show_pump == True:
    for i in range(len(u_test)):
        for j in range(len(dotV_vec_test)):
            u = u_test[i]
            V = dotV_vec_test[j]
            
            Deltap_pump[j,i] = mypump.calc_Deltap(u, V)
            
if show_valve == True:
    for i in range(len(kappa_test)):
        for j in range(len(dotV_vec_test)):
            kappa = kappa_test[i]
            V = dotV_vec_test[j]
            
            Deltap_valve[j,i] = myvalve.calc_Deltap(myfluid, kappa, V)
            
if show_pipe_hy == True:
    for i in range(len(dotV_vec_test)):
        V = dotV_vec_test[i]
        
        Deltap_pipe[i] = mypipe.calc_Deltap(myfluid, V)
        
if show_pipe_th == True:
    # const inlet temperature, vaiable volume flow
    T_in = T_in_pi_const_test
    for j in range(len(dotV_vec_pith_test)):
        V = dotV_vec_pith_test[j]
        T_out_pipe_Tfix[j], DeltaT_pipe_Tfix[j]  = mypipe.calc_temp(myfluid, V, T_in, T_soil = T_soil_test)
        
    V = dotV_pi_const_test
    for i in range(len(T_in_pi_test)):
        T_in = T_in_pi_test[i]
        T_out_pipe_Vfix[i], DeltaT_pipe_Vfix[i]  = mypipe.calc_temp(myfluid, V, T_in, T_soil = T_soil_test)
            
if show_hx == True:
    # T_in is fixed
    T_in_HTNW   = T_in_HTNW_const_test
    T_in_PSM    = T_in_PSM_const_test
    
    for i in range(len(dotV_hx_test)):
        for j in range(len(dotV_hx_test)):
            V_HTNW  = dotV_hx_test[j]
            V_PSM   = dotV_hx_test[i]
            
            temps_HTNW, temps_PSM = myhx.calc_temp(mu_test, myfluid, V_HTNW, V_PSM, T_in_HTNW, T_in_PSM, tol=mytolerance)
            
            T_out_hxHTNW_Tfix[j,i]   = temps_HTNW[0]
            T_out_hxPSM_Tfix[j,i]    = temps_PSM[0]
            DeltaT_hxHTNW_Tfix[j,i]  = temps_HTNW[1]
            DeltaT_hxPSM_Tfix[j,i]   = temps_PSM[1]
            
    # dotV is fixed
    V_HTNW  = dotV_HTNW_const_test
    V_PSM   = dotV_PSM_const_test
    
    for i in range(len(T_in_hx_test)):
        for j in range(len(T_in_hx_test)):
            T_in_HTNW   =   T_in_hx_test[j]
            T_in_PSM    =   T_in_hx_test[i]
            
            temps_HTNW, temps_PSM = myhx.calc_temp(mu_test, myfluid, V_HTNW, V_PSM, T_in_HTNW, T_in_PSM, tol=mytolerance)
            
            T_out_hxHTNW_Vfix[j,i]   = temps_HTNW[0]
            T_out_hxPSM_Vfix[j,i]    = temps_PSM[0]
            DeltaT_hxHTNW_Vfix[j,i]  = temps_HTNW[1]
            DeltaT_hxPSM_Vfix[j,i]   = temps_PSM[1]

############################################################################

############################ visualize results  ############################

def my_set_ticks(x,y):
    xticks = list(plt.gca().get_xticks())
    yticks = list(plt.gca().get_yticks())
    xmin = np.nanmin(x)
    ymin = np.nanmin(y)
    xmax = np.nanmax(x)
    ymax = np.nanmax(y)
    xticks.insert(0,xmin)
    xticks.append(xmax)
    yticks.insert(0,ymin)
    yticks.append(ymax)
    plt.gca().set_xticks(xticks)
    plt.gca().set_yticks(yticks)
    xticklabels = list(plt.gca().get_xticklabels())
    yticklabels = list(plt.gca().get_yticklabels())
    xticklabels[1:-1] = [str(' ') for i in range(len(xticklabels[1:-1]))]
    yticklabels[1:-1] = [str(' ') for i in range(len(yticklabels[1:-1]))]
    xticklabels[0] = '%3.0f'%xticks[0]
    yticklabels[0] = '%3.0f'%yticks[0]
    xticklabels[-1] = '%3.0f'%xticks[-1]
    yticklabels[-1] = '%3.0f'%yticks[-1]
    plt.gca().set_xticklabels(xticklabels)
    plt.gca().set_yticklabels(yticklabels)
    plt.gca().set_zticklabels([])
    
def my_set_ticks2(x,y):
    xticks = list(plt.gca().get_xticks())
    yticks = list(plt.gca().get_yticks())
    xmin = np.nanmin(x)
    ymin = np.nanmin(y)
    xmax = np.nanmax(x)
    ymax = np.nanmax(y)
    xticks.insert(0,xmin)
    xticks.append(xmax)
    yticks.insert(0,ymin)
    yticks.append(ymax)
    plt.gca().set_xticks(xticks)
    plt.gca().set_yticks(yticks)
    xticklabels = list(plt.gca().get_xticklabels())
    yticklabels = list(plt.gca().get_yticklabels())
    xticklabels[1:-1] = [str(' ') for i in range(len(xticklabels[1:-1]))]
    yticklabels[1:-1] = [str(' ') for i in range(len(yticklabels[1:-1]))]
    xticklabels[0] = '%3.0f'%xticks[0]
    yticklabels[0] = '%3.0f'%yticks[0]
    xticklabels[-1] = '%3.0f'%xticks[-1]
    yticklabels[-1] = '%3.0f'%yticks[-1]
    plt.gca().set_xticklabels(xticklabels)
    plt.gca().set_yticklabels(yticklabels)
    
def my_set_ticks3(x,y):
    xticks = list(plt.gca().get_xticks())
    yticks = list(plt.gca().get_yticks())
    xmin = np.nanmin(x)
    ymin = 0
    xmax = np.nanmax(x)
    ymax = np.nanmax(y)
    xticks.insert(0,xmin)
    xticks.append(xmax)
    yticks.insert(0,ymin)
    yticks.append(ymax)
    plt.gca().set_xticks(xticks)
    plt.gca().set_yticks(yticks)
    xticklabels = list(plt.gca().get_xticklabels())
    yticklabels = list(plt.gca().get_yticklabels())
    xticklabels[1:-1] = [str(' ') for i in range(len(xticklabels[1:-1]))]
    yticklabels[1:-1] = [str(' ') for i in range(len(yticklabels[1:-1]))]
    xticklabels[0] = '%3.0f'%xticks[0]
    yticklabels[0] = '%3.0f'%yticks[0]
    xticklabels[-1] = '%3.0f'%xticks[-1]
    yticklabels[-1] = '%3.0f'%yticks[-1]
    plt.gca().set_xticklabels(xticklabels)
    plt.gca().set_yticklabels(yticklabels)
    
def my_set_ticks4(x,y):
    xticks = list(plt.gca().get_xticks())
    yticks = list(plt.gca().get_yticks())
    xmin = np.nanmin(x)
    ymin = -1000
    xmax = np.nanmax(x)
    ymax = np.nanmax(y)
    xticks.insert(0,xmin)
    xticks.append(xmax)
    yticks.insert(0,ymin)
    yticks.append(ymax)
    plt.gca().set_xticks(xticks)
    plt.gca().set_yticks(yticks)
    xticklabels = list(plt.gca().get_xticklabels())
    yticklabels = list(plt.gca().get_yticklabels())
    xticklabels[1:-1] = [str(' ') for i in range(len(xticklabels[1:-1]))]
    yticklabels[1:-1] = [str(' ') for i in range(len(yticklabels[1:-1]))]
    xticklabels[0] = '%3.0f'%xticks[0]
    yticklabels[0] = '%3.0f'%yticks[0]
    xticklabels[-1] = '%3.0f'%xticks[-1]
    yticklabels[-1] = '%3.0f'%yticks[-1]
    plt.gca().set_xticklabels(xticklabels)
    plt.gca().set_yticklabels(yticklabels)
    
def my_set_ticks5(x,y):
    xticks = list(plt.gca().get_xticks())
    yticks = list(plt.gca().get_yticks())
    xmin = np.nanmin(x)
    ymin = np.nanmin(y)
    xmax = np.nanmax(x)
    ymax = np.nanmax(y)
    xticks.insert(0,xmin)
    xticks.append(xmax)
    yticks.insert(0,ymin)
    yticks.append(ymax)
    plt.gca().set_xticks(xticks)
    plt.gca().set_yticks(yticks)
    xticklabels = list(plt.gca().get_xticklabels())
    yticklabels = list(plt.gca().get_yticklabels())
    xticklabels[1:-1] = [str(' ') for i in range(len(xticklabels[1:-1]))]
    yticklabels[1:-1] = [str(' ') for i in range(len(yticklabels[1:-1]))]
    xticklabels[0] = '%3.0f'%xticks[0]
    yticklabels[0] = '%4.1f'%yticks[0]
    xticklabels[-1] = '%3.0f'%xticks[-1]
    yticklabels[-1] = '%4.1f'%yticks[-1]
    plt.gca().set_xticklabels(xticklabels)
    plt.gca().set_yticklabels(yticklabels)


if show_pump == True:
    plt.figure(num='Pumpe', figsize=myfigsize2)
    
    for i in range(len(u_test)):
        plt.plot(dotV_vec_test, Deltap_pump[:,i], label=''.join(('u=',str(u_test[i]))), \
        linewidth=mythickness, linestyle=mylinestyles[i], color=mylinecolors[i])
        
    #plt.plot(dotV_vec_test, np.zeros(len(dotV_vec_test)), color = 'k', linewidth=mythickness)
    
    plt.title('pump curve', fontsize=myfontsize, fontweight='bold')
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'pressure increase by pump $\Delta p\/\left[hPa\right]$', fontsize=myfontsize)
    plt.gca().set_ylim(ymin=0)
    plt.tick_params(axis='both', labelsize=myfontsize)
    my_set_ticks3(dotV_vec_test, Deltap_pump)


    plt.legend(fontsize = 0.8*myfontsize)
    # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # secax.tick_params(labelsize=myfontsize)
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pre_pump.', saveformat)), format = saveformat)
    
if show_valve == True:
    plt.figure(num='Ventil', figsize=myfigsize2)
    
    for i in range(len(kappa_test)):
        plt.plot(dotV_vec_test, Deltap_valve[:,i], label=''.join((r'$\kappa$=',str(kappa_test[i]))),\
        linewidth=mythickness, linestyle=mylinestyles[i], color=mylinecolors[i])
        
    #plt.plot(dotV_vec_test, np.zeros(len(dotV_vec_test)), color = 'k', linewidth=mythickness)
    my_set_ticks4(dotV_vec_test, Deltap_valve)
    plt.title('control valve curve', fontsize=myfontsize, fontweight='bold')
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'pressure head $\Delta p\/\left[hPa\right]$', fontsize=myfontsize)
    plt.gca().set_ylim(ymax=+10)
    plt.gca().set_ylim(ymin=-1000)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend(fontsize = 0.8*myfontsize)
    # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # secax.tick_params(labelsize=myfontsize)
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pres_valve.', saveformat)), format = saveformat)
    
if show_pipe_hy == True:
    plt.figure(num='Rohr hydraulisch', figsize=myfigsize2)
    
    plt.plot(dotV_vec_test, Deltap_pipe, linewidth=mythickness, color=mylinecolors[-1], linestyle=mylinestyles[-1])

    plt.suptitle('hydraulic pipe curve', fontsize=myfontsize, fontweight='bold')
    plt.title(''.join(('L=',str(mypipe.L),'m,',r'$~~~\zeta^{inst}=$',str(mypipe.zeta_instal))), fontsize=math.floor(myfontsize*0.9))
    my_set_ticks2(dotV_vec_test, Deltap_pipe)
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'pressure head $\Delta p\/\left[hPa\right]$', fontsize=myfontsize)
    plt.gca().set_ylim(ymax=+10)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()

    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pres_pipe_hy.', saveformat)), format = saveformat)
    
if show_pipe_th == True:
    fig=plt.figure(num='Rohr thermisch', figsize=myfigsize)
    #plt.suptitle('thermal pipe curve (insulated)', fontsize=myfontsize, fontweight='bold')
    
    ax1 = fig.add_subplot(1,2,1)
    plt.plot(dotV_vec_pith_test, T_out_pipe_Tfix-273.15, \
    linewidth=mythickness, linestyle=mylinestyles[-1], color=mylinecolors[-1])
    description1 = r'$T_{in} = %4.1f ^\circ C,$ $L= %4.0f m,$ $T_{soil}= %4.1f ^\circ C$'%(T_in_pi_const_test-273.15, mypipe.L, T_soil_test-273.15)
    plt.title('over varying volume flows\n'+description1, fontsize=math.floor(myfontsize*0.9), fontweight = 'bold')
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'outlet temperature $T_{out} ~\left[^\circ C\right]$', fontsize=myfontsize)
    plt.tick_params(axis='both', labelsize=myfontsize)
    my_set_ticks5(dotV_vec_pith_test, T_out_pipe_Tfix-273.15)
    plt.legend()
    
    ax2 = fig.add_subplot(1,2,2)
    plt.plot(T_in_pi_test-273.15, T_out_pipe_Vfix-273.15, \
    linewidth=mythickness, linestyle=mylinestyles[-1], color=mylinecolors[-1])
    description2 = r'$\dot{V} = %4.1f \frac{l}{min},$ $L= %4.0f m,$ $T_{soil}= %4.1f ^\circ C$'%(dotV_pi_const_test, mypipe.L, T_soil_test-273.15)
    plt.title('over varying inlet temperatures\n'+description2, fontsize=math.floor(myfontsize*0.9), fontweight='bold')
    
    plt.gca().set_xlabel(r'inlet temperature $T_{in}\/\left[ ^\circ C\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'outlet temperature $T_{out} ~\left[^\circ C\right]$', fontsize=myfontsize)
    plt.tick_params(axis='both', labelsize=myfontsize)
    my_set_ticks2(T_in_pi_test-273.15, T_out_pipe_Vfix-273.15)
    plt.legend()
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pre_pipe_th.', saveformat)), format = saveformat)
    
if show_hx == True:
    
    #mysuptitle = r'heatexchanger'
    #fig.suptitle(mysuptitle, fontsize=myfontsize, fontweight='bold')
    
    # outlet temperatures over varying inlet temperatures with fixed volume flows
    X,Y = np.meshgrid(T_in_hx_test-273.15, T_in_hx_test-273.15)
    Z1 = T_out_hxPSM_Vfix
    Z2 = T_out_hxPSM_Tfix
    myMin=np.nanmin([np.nanmin(Z1-273.15),np.nanmin(Z2-273.15)])
    myMax=np.nanmax([np.nanmax(Z1-273.15),np.nanmax(Z2-273.15)])
    cmap = mpl.cm.coolwarm
    mynorm = mpl.colors.Normalize(vmin=myMin, vmax=myMax)
    fig = plt.figure(num='Waermetauscher', figsize=myfigsize)

    ax1 = fig.add_subplot(1,2,1, projection='3d')
    surf1 = ax1.plot_surface(X,Y,Z1-273.15,cmap=cmap, norm=mynorm, linewidth=0,antialiased=True)
    ax1.view_init(elev=15, azim=45)
    plt.tick_params(axis='both', labelsize=myfontsize)
    description1 = r'$\dot{V}_{HTNW}= %8.1f \frac{l}{min},\/\dot{V}_{PSM}= %8.1f \frac{l}{min}$'%(dotV_HTNW_const_test, dotV_PSM_const_test)
    plt.title('over vayring inlet temperatures\n'+ description1, fontsize=myfontsize, fontweight='bold')
    plt.gca().set_xlabel(r'$T_{HTNW,in}\/\left[ ^\circ C\right]$', fontsize=myfontsize, labelpad=mylabelpad)
    plt.gca().set_ylabel(r'$T_{PSM,in}\/\left[ ^\circ C\right]$', fontsize=myfontsize, labelpad=mylabelpad)
    plt.gca().set_zlabel(r'$T_{PSM,out}\/\left[^\circ C\right]$', fontsize=myfontsize)
    my_set_ticks(X,Y)
    
    # outlet temperatures over varying volume flows with fixed inlet temperatures
    X,Y = np.meshgrid(dotV_hx_test, dotV_hx_test)
    
    ax2 = fig.add_subplot(1,2,2, projection='3d')
    surf2 = ax2.plot_surface(X,Y,Z2-273.15,cmap=cmap, norm=mynorm, linewidth=0,antialiased=True)
    ax2.view_init(elev=15, azim=45)
    my_set_ticks(X,Y)
    plt.tick_params(axis='both', labelsize=myfontsize)
    description2 = r'$\vartheta_{HTNW,in}= %8.1f ^\circ C,\/\vartheta_{PSM,in}= %8.1f ^\circ C$'%(T_in_HTNW_const_test -273.15, T_in_PSM_const_test -273.15)
    plt.title('over varying volume flows\n' + description2, fontsize=myfontsize, fontweight='bold')
    plt.gca().set_xlabel(r'$\dot{V}_{HTNW}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize, labelpad=mylabelpad)
    plt.gca().set_ylabel(r'$\dot{V}_{PSM}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize, labelpad=mylabelpad)
    plt.gca().set_zlabel(r'$T_{PSM,out}\/\left[^\circ C\right]$', fontsize=myfontsize)
 
    
    plt.legend(fontsize = 0.8*myfontsize)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar=fig.colorbar(surf2, cax=cbar_ax, shrink=0.5, aspect=5)
    cbar.ax.tick_params(labelsize=myfontsize*0.8)
    plt.show(block=False)
    plt.savefig(''.join((savepath, 'pres_hx.', saveformat)), format = saveformat)
    

    #--------------------------------------


