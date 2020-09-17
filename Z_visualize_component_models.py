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

myfontsize = 18
mythickness = 3
mylinestyles = [':', '-.', '--', '-']
mylinecolors = ['r','g','b','k']
savepath = './figs/'
saveformat = 'png' # jpg, png, eps, svgd
show_pump       = True
show_valve      = True
show_pipe_hy    = True
show_pipe_th    = True
show_hx         = True
mytolerance     = 10**(-2)

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
T_in_PSM_const_test     = 45+273.15
T_in_HTNW_const_test    = 65+273.15
T_soil_test     =   12+273.15
T_in_pi_test    =   np.linspace(20,90,4)+273.15
T_in_hx_test   =   np.linspace(20,90,300)+273.15
dotV_hx_test   =   np.linspace(1,30,300)



############################################################################

############################ initialize result vectors ######################

Deltap_pump     =   np.zeros((len(dotV_vec_test),len(u_test)), dtype=float)
Deltap_valve    =   np.zeros((len(dotV_vec_test),len(kappa_test)), dtype=float)
Deltap_pipe     =   np.zeros(len(dotV_vec_test), dtype=float)
DeltaT_pipe     =   np.zeros((len(dotV_vec_pith_test), len(T_in_pi_test)), dtype=float)
T_out_pipe      =   np.zeros((len(dotV_vec_pith_test), len(T_in_pi_test)), dtype=float)
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
    for i in range(len(T_in_pi_test)):
        for j in range(len(dotV_vec_pith_test)):
            T_in = T_in_pi_test[i]
            V = dotV_vec_pith_test[j]
        
            T_out_pipe[j,i], DeltaT_pipe[j,i]  = mypipe.calc_temp(myfluid, V, T_in, T_soil = T_soil_test)
            
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

if show_pump == True:
    plt.figure(num='Pumpe', figsize=[2*21.0/2.54, 2*12.98/2.54])
    
    for i in range(len(u_test)):
        plt.plot(dotV_vec_test, Deltap_pump[:,i], label=''.join(('u=',str(u_test[i]))), \
        linewidth=mythickness, linestyle=mylinestyles[i], color=mylinecolors[i])
        
    #plt.plot(dotV_vec_test, np.zeros(len(dotV_vec_test)), color = 'k', linewidth=mythickness)
    
    plt.title('pump curve', fontsize=myfontsize, fontweight='bold')
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'pressure increase by pump $\Delta p\/\left[hPa\right]$', fontsize=myfontsize)
    plt.gca().set_ylim(ymin=-10)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()
    # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # secax.tick_params(labelsize=myfontsize)
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pump.', saveformat)), format = saveformat)
    
if show_valve == True:
    plt.figure(num='Ventil', figsize=[2*21.0/2.54, 2*12.98/2.54])
    
    for i in range(len(kappa_test)):
        plt.plot(dotV_vec_test, Deltap_valve[:,i], label=''.join((r'$\kappa$=',str(kappa_test[i]))),\
        linewidth=mythickness, linestyle=mylinestyles[i], color=mylinecolors[i])
        
    #plt.plot(dotV_vec_test, np.zeros(len(dotV_vec_test)), color = 'k', linewidth=mythickness)
    
    plt.title('control valve curve', fontsize=myfontsize, fontweight='bold')
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'pressure head $\Delta p\/\left[hPa\right]$', fontsize=myfontsize)
    plt.gca().set_ylim(ymax=+10)
    plt.gca().set_ylim(ymin=-1000)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()
    # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # secax.tick_params(labelsize=myfontsize)
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'valve.', saveformat)), format = saveformat)
    
if show_pipe_hy == True:
    plt.figure(num='Rohr hydraulisch', figsize=[2*21.0/2.54, 2*12.98/2.54])
    
    plt.plot(dotV_vec_test, Deltap_pipe, linewidth=mythickness, color=mylinecolors[-1], linestyle=mylinestyles[-1])
        
    #plt.plot(dotV_vec_test, np.zeros(len(dotV_vec_test)), color = 'k', linewidth=mythickness)
    
    plt.suptitle('hydraulic pipe curve', fontsize=myfontsize, fontweight='bold')
    plt.title(''.join(('L=',str(mypipe.L),'m,',r'$~~~\zeta^{inst}=$',str(mypipe.zeta_instal))), fontsize=math.floor(myfontsize*0.9))
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'pressure head $\Delta p\/\left[hPa\right]$', fontsize=myfontsize)
    plt.gca().set_ylim(ymax=+10)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()
    # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # secax.tick_params(labelsize=myfontsize)
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pipe_hy.', saveformat)), format = saveformat)
    
if show_pipe_th == True:
    fig=plt.figure(num='Rohr thermisch1', figsize=[2*21.0/2.54, 2*12.98/2.54])
    plt.suptitle('thermal pipe curve (insulated)', fontsize=myfontsize, fontweight='bold')
    # ---------------------------------
    #fig.add_subplot(1,2,1)
    
    for i in range(len(T_in_pi_test)):
        T_in_now = '%8.1f' % (T_in_pi_test[i]-273.15)
        plt.plot(dotV_vec_pith_test, T_out_pipe[:,i]-273.15, label=''.join((r'$T_{in}$=',T_in_now)), \
        linewidth=mythickness, linestyle=mylinestyles[i], color=mylinecolors[i])
    
    plt.title(''.join(('outlet temperature, ','L=',str(mypipe.L),'m')), fontsize=math.floor(myfontsize*0.9))
    
    plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'outlet temperature $T_{out} ~\left[^\circ C\right]$', fontsize=myfontsize)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()
    # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # secax.tick_params(labelsize=myfontsize)
    
    plt.show(block=False)
    
    plt.savefig(''.join((savepath, 'pipe_th.', saveformat)), format = saveformat)
    #--------------------------------------
    # fig.add_subplot(1,2,2)
    # for i in range(len(T_in_pi_test)):
        # T_in_now = '%8.1f' % (T_in_pi_test[i]-273.15)
        # plt.plot(dotV_vec_pith_test, DeltaT_pipe[:,i], label=''.join((r'$T_{in}$=',T_in_now)), \
        # linewidth=mythickness, linestyle=mylinestyles[i], color=mylinecolors[i])
    
    # plt.title(''.join(('temperature difference, ','L=',str(mypipe.L),'m')), fontsize=math.floor(myfontsize*0.9))
    
    # plt.gca().set_xlabel(r'volume flow $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # plt.gca().set_ylabel(r'temperature difference $\Delta T~\left[^\circ C\right]$', fontsize=myfontsize)
    # plt.tick_params(axis='both', labelsize=myfontsize)
    # plt.legend()
    # # secax = plt.gca().secondary_xaxis('top', functions=(lambda x: x*60, lambda x: x/60))
    # # secax.set_xlabel(r'volume flux $\dot{V}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    # # secax.tick_params(labelsize=myfontsize)
    
    # plt.show(block=False)
    
    # plt.savefig(''.join((savepath, 'pipe_th.', saveformat)), format = saveformat)
    
if show_hx == True:
    #--------------------------------------
    # outlet temperatures over varying volume flows with fixed inlet temperatures
    X,Y = np.meshgrid(dotV_hx_test, dotV_hx_test)
    Z1 = T_out_hxHTNW_Tfix
    Z2 = T_out_hxPSM_Tfix
    myMin=np.nanmin([np.nanmin(Z1-273.15),np.nanmin(Z2-273.15)])
    myMax=np.nanmax([np.nanmax(Z1-273.15),np.nanmax(Z2-273.15)])
    cmap = mpl.cm.coolwarm
    mynorm = mpl.colors.Normalize(vmin=myMin, vmax=myMax)
    fig = plt.figure(num='Waermetauscher fix T_in', figsize=[2*21.0/2.54, 2*12.98/2.54])
    mysuptitle = r'heatexchanger $\vartheta_{HTNW,in}= %8.1f ^\circ C,\/\vartheta_{PSM,in}= %8.1f ^\circ C $'%(T_in_HTNW_const_test -273.15, T_in_PSM_const_test -273.15)
    fig.suptitle(mysuptitle, fontsize=myfontsize, fontweight='bold')
    ax1 = fig.add_subplot(1,2,1, projection='3d')
    surf1 = ax1.plot_surface(X,Y,Z1-273.15,cmap=cmap, norm=mynorm, linewidth=0,antialiased=True)
    ax1.view_init(elev=15, azim=45)
    plt.tick_params(axis='both', labelsize=myfontsize)
    #fig.colorbar(surf1, shrink=0.5, aspect=5)
    plt.title(r'heat network side (source)', fontsize=myfontsize, fontweight='bold')
    plt.gca().set_xlabel(r'$\dot{V}_{HTNW}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'$\dot{V}_{PSM}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_zlabel(r'$T_{HTNW,out}\/\left[^\circ C\right]$', fontsize=myfontsize)
    ax2 = fig.add_subplot(1,2,2, projection='3d')
    surf2 = ax2.plot_surface(X,Y,Z2-273.15,cmap=cmap, norm=mynorm, linewidth=0,antialiased=True)
    ax2.view_init(elev=15, azim=45)
    #fig.colorbar(surf2, shrink=0.5, aspect=5)
    plt.title(r'prosumer side (sink)', fontsize=myfontsize, fontweight='bold')
    plt.gca().set_xlabel(r'$\dot{V}_{HTNW}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'$\dot{V}_{PSM}\/\left[ \frac{l}{min}\right]$', fontsize=myfontsize)
    plt.gca().set_zlabel(r'$T_{PSM,out}\/\left[^\circ C\right]$', fontsize=myfontsize)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()
    plt.show(block=False)
    plt.savefig(''.join((savepath, 'hx1.', saveformat)), format = saveformat)
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(surf1, cax=cbar_ax, shrink=0.5, aspect=5)
    # --------------------------------------
    
    # outlet temperatures over varying inlet temperatures with fixed volume flows
    X,Y = np.meshgrid(T_in_hx_test-273.15, T_in_hx_test-273.15)
    Z1 = T_out_hxHTNW_Vfix
    Z2 = T_out_hxPSM_Vfix
    myMin=np.nanmin([np.nanmin(Z1-273.15),np.nanmin(Z2-273.15)])
    myMax=np.nanmax([np.nanmax(Z1-273.15),np.nanmax(Z2-273.15)])
    cmap = mpl.cm.coolwarm
    mynorm = mpl.colors.Normalize(vmin=myMin, vmax=myMax)
    fig = plt.figure(num='Waermetauscher fix dotV', figsize=[2*21.0/2.54, 2*12.98/2.54])
    mysuptitle = r'heatexchanger $\dot{V}_{HTNW}= %8.1f \frac{l}{min},\/\dot{V}_{PSM}= %8.1f \frac{l}{min} $'%(dotV_HTNW_const_test, dotV_PSM_const_test)
    fig.suptitle(mysuptitle, fontsize=myfontsize, fontweight='bold')
    ax1 = fig.add_subplot(1,2,1, projection='3d')
    surf1 = ax1.plot_surface(X,Y,Z1-273.15,cmap=cmap, norm=mynorm, linewidth=0,antialiased=True)
    ax1.view_init(elev=15, azim=45)
    plt.tick_params(axis='both', labelsize=myfontsize)
    #fig.colorbar(surf1, shrink=0.5, aspect=5)
    plt.title(r'heat network side (source)', fontsize=myfontsize, fontweight='bold')
    plt.gca().set_xlabel(r'$T_{HTNW,in}\/\left[ ^\circ C\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'$T_{PSM,in}\/\left[ ^\circ C\right]$', fontsize=myfontsize)
    plt.gca().set_zlabel(r'$T_{HTNW,out}\/\left[^\circ C\right]$', fontsize=myfontsize)
    ax2 = fig.add_subplot(1,2,2, projection='3d')
    surf2 = ax2.plot_surface(X,Y,Z2-273.15,cmap=cmap, norm=mynorm, linewidth=0,antialiased=True)
    ax2.view_init(elev=15, azim=45)
    #fig.colorbar(surf2, shrink=0.5, aspect=5)
    plt.title(r'prosumer side (sink)', fontsize=myfontsize, fontweight='bold')
    plt.gca().set_xlabel(r'$T_{HTNW,in}\/\left[ ^\circ C\right]$', fontsize=myfontsize)
    plt.gca().set_ylabel(r'$T_{PSM,in}\/\left[ ^\circ C\right]$', fontsize=myfontsize)
    plt.gca().set_zlabel(r'$T_{PSM,out}\/\left[^\circ C\right]$', fontsize=myfontsize)
    plt.tick_params(axis='both', labelsize=myfontsize)
    plt.legend()
    plt.show(block=False)
    plt.savefig(''.join((savepath, 'hx2.', saveformat)), format = saveformat)
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(surf2, cax=cbar_ax, shrink=0.5, aspect=5)
    #--------------------------------------


