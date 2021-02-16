import math
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pdb
 
class fluid:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all 
    rights reserved
    rho             [kg/(m**3)]     density of fluid
    cp              [J/(kg*K)]      specific isobaric heat capacity
    mu              [Pa*s]          dynamic viscosity'''
    
    # Initializer / Instance Attributes
    def __init__(self, rho=1000, cp=4200, mu=(1.0016*(10**(-3)))):
        self.rho_SI    =   rho     # density of fluid                  # [kg/(m**3)]
        self.cp_SI     =   cp      # specific isobaric heat capacity   # [J/(kg*K)]
        self.mu_SI     =   mu      # dynamic viscosity                 # [Pa*s]
        self.rho    =   rho/1000    # density of fluid                  # [kg/l]
        self.cp     =   cp          # specific isobaric heat capacity   # [J/(kg*K)]
        self.mu     =   mu/(100*60) # dynamic viscosity                 # [hPa*min]
        
    # Define function for showing properties
    def props(self):
        return dict(self.__dict__)
 
class heatexchanger:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all 
    rights reserved
    dotV_nom        [l/min]         nominal volume flow rate
    Deltap_nom      [hPa]           head loss at nominal mass flow rate
    k_nom           [W/((m*2)*K)]   nominal heat transfer coefficient
    A               [m**2]          heat transfer surface'''
    
    # Initializer / Instance Attributes
    def __init__(self, dotV_nom=21.504, Deltap_nom=155, k_nom = 5270, A = 1.13):
        
        # set properties
        self.dotV_nom   =   dotV_nom    # nominal volume flow rate    [l/min]
        self.Deltap_nom =   abs(Deltap_nom) # head loss at 
                                            # nominal mass flow rate   [hPa]
        self.k_nom      =   k_nom       # nominal heat transfer coefficient
                                        # [W/((m**2)*K)]
        self.A          =   A           # heat transfer surface [m**2]
        
    def hy_params(self):
        ## hydraulic paramter
        a_hx_1  =   -(self.Deltap_nom)/((self.dotV_nom)**2)      # [hPa/((l/min)**2)]
        self.a_hx_1 =   a_hx_1
        return a_hx_1
        
    def th_params(self, mode, dotV_HTNW, dotV_PSM, fluid, 
                 tol=10**(-2)):
        ''' mode        [-] -1 = consumer, +1 = producer
            dotV_HTNW   [l/min]     volume flow on the network side
            dotV_PSM    [l/min]     volume flow on the prosumer side
            fluid                   defined fluid'''
        
        # only consider positive volume flows
        dotV_HTNW   =   abs(dotV_HTNW)/(6*(10**4))
        dotV_PSM    =   abs(dotV_PSM)/(6*(10**4))
        
        # allocate variables according to mode (producer/consumer)
        if mode == 1: # producer
            dotV_1 = dotV_PSM
            dotV_2 = dotV_HTNW
        elif mode == -1: # consumer
            dotV_1 = dotV_HTNW
            dotV_2 = dotV_PSM
        else:
            raise ValueError('mode must either be +1 or -1')
        
        # calculate substitution variables
        dotm_1  =   dotV_1 * fluid.rho_SI  # [kg/min]
        dotm_2  =   dotV_2 * fluid.rho_SI
        dotW_1     =   dotm_1 * fluid.cp_SI # [W/K]
        dotW_2     =   dotm_2 * fluid.cp_SI # [W/K]
        R_1     =   dotW_1/dotW_2
        R_2     =   dotW_2/dotW_1
        NTU_1   =   (self.k_nom * self.A) / dotW_1  # [-]
        NTU_2   =   (self.k_nom * self.A) / dotW_2  # [-]
        

        if abs(R_1-1) >= tol:
            try:
                X_1 = (1-math.exp((R_1-1)*NTU_1))/(1-R_1 * math.exp((R_1-1)*NTU_1))
            except OverflowError:
                X_1 = float('inf')
            try:
                X_2 = (1-math.exp((R_2-1)*NTU_2))/(1-R_2 * math.exp((R_2-1)*NTU_2))
            except OverflowError:
                X_2 = float('inf')
        else:
            try:
                X_1 = (NTU_1)/(1+NTU_1)
            except OverflowError:
                X_1 = float('inf')
            try:
                X_2 = (NTU_2)/(1+NTU_2)
            except OverflowError:
                X_2 = float('inf')
            
        b_11 = 1-X_1
        b_12 = -1
        b_13 = X_1
        b_14 = 0
        b_21 = -X_2
        b_22 = 0
        b_23 = X_2-1
        b_24 = 1
        
        # if (abs(b_11/b_13) > 1.5 or abs(b_11/b_13) < 0.5) and (dotV_1 > 1 and dotV_2 > 1):
            # pdb.set_trace()
        # elif (abs(b_21/b_23) > 1.5 or abs(b_21/b_23) < 0.5) and (dotV_1 > 1 and dotV_2 > 1):
            # pdb.set_trace()
        
        params_PSM={}
        params_HTNW={}
        
        # reallocate variables
        if mode == 1: # producer
            params_PSM['T_in_PSM'] = b_11
            params_PSM['T_out_PSM'] = b_12
            params_PSM['T_in_HTNW'] = b_13
            params_PSM['T_out_HTNW'] = b_14
            params_HTNW['T_in_PSM'] = b_21
            params_HTNW['T_out_PSM'] = b_22
            params_HTNW['T_in_HTNW'] = b_23
            params_HTNW['T_out_HTNW'] = b_24
        elif mode == -1: # consumer
            params_HTNW['T_in_HTNW'] = b_11
            params_HTNW['T_out_HTNW'] = b_12
            params_HTNW['T_in_PSM'] = b_13
            params_HTNW['T_out_PSM'] = b_14
            params_PSM['T_in_HTNW'] = b_21
            params_PSM['T_out_HTNW'] = b_22
            params_PSM['T_in_PSM'] = b_23
            params_PSM['T_out_PSM'] = b_24
        else:
            raise ValueError('mode must either be +1 or -1')
            
        self.params_HTNW    =   params_HTNW
        self.params_PSM     =   params_PSM
        
        return params_HTNW, params_PSM       
        
    # Define function for showing properties
    def props(self):
        return dict(self.__dict__)
        
    def calc_Deltap(self, dotV):
        
        a_hx_1 = self.hy_params()
        
        Deltap = a_hx_1 * (dotV**2)
        
        return Deltap
        
    def calc_temp(self, mode, fluid, dotV_HTNW, dotV_PSM, T_in_HTNW, T_in_PSM, tol=10**(-2)):
    
        params_HTNW, params_PSM = self.th_params(mode, dotV_HTNW, dotV_PSM, fluid, 
                 tol=10**(-2))
        
        if mode == 1: # producer
            b_11 = params_PSM['T_in_PSM']
            b_12 = params_PSM['T_out_PSM']
            b_13 = params_PSM['T_in_HTNW']
            b_14 = params_PSM['T_out_HTNW']
            b_21 = params_HTNW['T_in_PSM']
            b_22 = params_HTNW['T_out_PSM']
            b_23 = params_HTNW['T_in_HTNW']
            b_24 = params_HTNW['T_out_HTNW']
        elif mode == -1: # consumer
            b_11 = params_HTNW['T_in_HTNW']
            b_12 = params_HTNW['T_out_HTNW']
            b_13 = params_HTNW['T_in_PSM']
            b_14 = params_HTNW['T_out_PSM']
            b_21 = params_PSM['T_in_HTNW']
            b_22 = params_PSM['T_out_HTNW']
            b_23 = params_PSM['T_in_PSM']
            b_24 = params_PSM['T_out_PSM']
        else:
            raise ValueError('mode must either be +1 or -1')
        
        # only consider positive volume flows
        dotV_HTNW   =   abs(dotV_HTNW)
        dotV_PSM    =   abs(dotV_PSM)
        
        # allocate variables according to mode (producer/consumer)
        if mode == 1: # producer
            dotV_1 = dotV_PSM
            dotV_2 = dotV_HTNW
            T_in_1 = T_in_PSM
            T_in_2 = T_in_HTNW
        elif mode == -1: # consumer
            dotV_1 = dotV_HTNW
            dotV_2 = dotV_PSM
            T_in_1 = T_in_HTNW
            T_in_2 = T_in_PSM
        else:
            raise ValueError('mode must either be +1 or -1')
        
        a = 1/b_12 #1/(b_12-((b_14*b_22)/(b_24)))
        b = -b_11 #((b_14*b_21)/(b_24))-b_11
        c = -b_13 #((b_14*b_23)/(b_24))-b_13
        
        T_out_1 = a*(T_in_1*b+T_in_2*c)
        T_out_2 = (-1/b_24) * (b_21*T_in_1 + b_23*T_in_2)# -((b_21*T_in_1+b_22*T_out_1+b_23*T_in_2)/(b_24))
        
        if mode == 1: # producer
            T_out_PSM = T_out_1
            T_out_HTNW = T_out_2
        elif mode == -1: # consumer
            T_out_PSM = T_out_2
            T_out_HTNW = T_out_1
        else:
            raise ValueError('mode must either be +1 or -1')
        
        DeltaT_HTNW = T_out_HTNW - T_in_HTNW
        DeltaT_PSM  = T_out_PSM - T_in_PSM
        
        temps_HTNW = [T_out_HTNW, DeltaT_HTNW]
        temps_PSM  = [T_out_PSM, T_out_HTNW]
        
        return temps_HTNW, temps_PSM
      
class pipe:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all
    rights reserved
    L               [m]             pipe length
    d_hy            [m]             (hydraulic) inner pipe diameter
    epsilon         [mm]            roughness of inner pipe surface
    u_nom           [m/s]           nominal flow velocity
    zeta_instal     [-]             pressure loss coefficient (Druckverlustbeiwert)
                                    for installations in pipe (such as bends, junctions etc.) 
    N_layers        [-]             number of layers around pipe (for thermal losses),
                                    e.g. pipe wall (1), insulation (2), soil environment (3)
    d_layers        [m]             diameters [m] at the borders of the layers
    lambda_layers   [W/(m*K)]       thermal conductivity of layers
    h_ir            [W/((m**2)*K)]  heat transfer coefficient at inside of pipe
    h_or            [W/((m**2)*K)]  heat transfer coefficient at outside of outest pipe layer'''
    
    # Initializer / Instance Attributes
    def __init__(self, L=15,d_hy=0.022,epsilon=0.011,u_nom=0.5,zeta_instal=10,
                 N_layers = 2, d_layers = [0.022, 0.024, 0.062],
                 lambda_layers = [395, 0.04], h_ir=3500, h_or=float('inf')):
    
        # set properties
        self.L          =   L               # pipe length   [m]
        self.d_hy       =   d_hy            # hydraulic pipe diameter   [m]
        self.epsilon    =   epsilon*(1/1000)# pipe roughness    [m]
        self.u_nom      =   u_nom           # nominal flow speed in pipe [m/s]
        self.zeta_instal=   zeta_instal     # pressure loss coefficient for 
                                            # installations [-]
        self.N_layers   =   N_layers        # number of layers forming the pipe [-]
                                            # e.g. pipe-wall (1), insulation (2)
        self.d_layers   =   d_layers        # diameters [m] at the border of the layers
        self.lambda_layers  =   lambda_layers   # thermal conductivity of layers
                                                # [W/(m*K)]
        self.h_ir       =   h_ir            # [W/((m**2)*K)]
        self.h_or       =   h_or            # [W/((m**2)*K)]
                                                
        if len(self.d_layers)!=(self.N_layers+1):
            raise ValueError('Length of d_layers must be equal to (N_layers+1)!')
        elif len(self.lambda_layers) != self.N_layers:
            raise ValueError('Length of lambda_layers must be equal to N_layers!')
        else:
            pass
            
        if self.d_hy!=self.d_layers[0]:
            raise ValueError('Hydraulic diameter must be the same as innerst layer thermal diameter!')
        else:
            pass
        
        h_layers_rev = 0
        for i in range(self.N_layers):
                h_layers_rev += ((math.log((self.d_layers[i+1])/(self.d_layers[i])))/
                (2*math.pi*self.lambda_layers[i]))
        self.h_layers_rev = h_layers_rev
        
        self.k_th_ir = (1/(self.h_ir*math.pi*self.d_layers[0])+
                        h_layers_rev+
                        +1/(self.h_or*math.pi*self.d_layers[-1]))**(-1) # [W/((m)*K)]
                        
        self.R_th_ir = 1/(self.k_th_ir) # [(m*K)/(W)]
    
    def hy_params(self, fluid):
        '''fluid                   defined fluid'''
    
        A_hy    =   (1.0/4.0) * math.pi * (self.d_hy**2) # [m^2]
        
        Re_nom  =   (fluid.rho_SI * self.u_nom * self.d_hy) / (fluid.mu_SI)   #   [-]
        
        if Re_nom <= 2000: #Hagen-Poiseuille flow
            f_D = 64 / Re_nom
        elif Re_nom >= 4000:# Swamee and Jain 1976;
            f_D = 0.25/((math.log10((self.epsilon/(3.7*self.d_hy))+(5.74/(Re_nom**0.9))))**2)
        elif ( Re_nom >2000 ) and ( Re_nom < 4000 ): # transition region, linear interpolation between laminar and turbulent
            x1 = 2000
            x2 = 4000
            y1 = 64 / x1
            y2 = 0.25/((math.log10((self.epsilon/(3.7*self.d_hy))+(5.74/(x2**0.9))))**2)
            m = (y2-y1)/(x2-x1)
            t = y1 - m * x1
            f_D = m * Re_nom + t
        
        # moody equation - old, artifact
        # f_D   =   0.0055 * (1+(2*(10**4)*(self.epsilon/self.d_hy)+(10**6/
        # Re_nom))**(1.0/3.0))     # [-]
        
        # hydralic resistance of pipe
        a_pi_1    =   -(8/(math.pi**2))*fluid.rho_SI*(1/(self.d_hy**4))*\
            ((self.L/self.d_hy)*f_D+self.zeta_instal)# [Pa*(s**2)/(m**3**2)]
        a_pi_1    =   a_pi_1 * 1/(100*(60*1000)**2) #[Pa*(min**2)/(l**2)]
        
        return  a_pi_1
        
    def th_params(self, fluid, dotV, T_soil = 12+273.15):
        ''' fluid       defined fluid
            dotV        volume flow [l/min]
            T_soil      constant temperature of surrounding soil [K]'''
        
        dotV2 = abs(dotV)*(1/60000)
        
        A_circ = math.pi * self.d_hy * self.L
        
        try:
            
            s_pi = -(self.L/(fluid.rho_SI * fluid.cp_SI * self.R_th_ir * dotV2))
            b_pi_1 = math.exp(s_pi)
            b_pi_3 = (math.exp(s_pi)-1)*T_soil
            b_pi_2 = -1
            
        except:
            s_pi = float('inf')
            b_pi_1 = 0
            b_pi_3 = 0
            b_pi_2 = 0
                      
        return [b_pi_1, b_pi_2, b_pi_3]
                                                    
    # Define function for showing properties
    def props(self):
        return dict(self.__dict__)
        
    def calc_Deltap(self,fluid,dotV):
    
        g = 9.81 # [m/(s**2)]
        
        a_pi_1 = self.hy_params(fluid)
        
        Deltap = a_pi_1 * (dotV**2)
        
        return Deltap
        
    def calc_temp(self, fluid, dotV, T_in, T_soil = 12+273.15):
        
        T_soil_loc = T_soil
    
        [b_pi_1, b_pi_2, b_pi_3] = self.th_params(fluid, dotV, T_soil = T_soil_loc)
        
        T_out = (b_pi_3-b_pi_1*T_in)/(b_pi_2)
        
        DeltaT = T_out - T_in
        
        return T_out, DeltaT
        
class pump:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all
    rights
    reserved
    n_nom           [RPM = 1/min]   nominal pump speed
    u_ref_1         [-]             controll variable ref. oper. state 1
    dotV_ref_1      [l/min]         volume flow ref. oper. state 1
    Deltap_ref_1    [hPa]           pressure difference ref. oper. state 1
    u_ref_2         [-]             controll variable ref. oper. state 2
    dotV_ref_2      [l/min]         volume flow ref. oper. state 2
    Deltap_ref_2    [hPa]           pressure difference ref. oper. state 2
    check_valve     [-]             indicator whether pump has a check valve in serial
    K_vs_cv         [m³/h]          flow coefficient of check valve in flow direction'''
    
    # Initializer / Instance Attributes
    def __init__(self, n_nom=4100, u_ref_1=1, dotV_ref_1=0, Deltap_ref_1=402.21, 
                    u_ref_2=1, dotV_ref_2=55.33, Deltap_ref_2=0, check_valve = True, K_vs_cv = 6):
            
        self.n_nom  =   n_nom
        self.u_ref_1    =   u_ref_1
        self.dotV_ref_1 =   dotV_ref_1
        self.Deltap_ref_1   =   Deltap_ref_1
        self.u_ref_2    =   u_ref_2
        self.dotV_ref_2 =   dotV_ref_2
        self.Deltap_ref_2   =   Deltap_ref_2
        self.K_vs_cv = K_vs_cv
        self.check_valve = check_valve
    
    def hy_params(self, fluid):
        
        a_pu_1  =   (self.Deltap_ref_1-((self.u_ref_1/self.u_ref_2)**2)*self.Deltap_ref_2)/\
            ((self.dotV_ref_1**2)-((self.u_ref_1/self.u_ref_2)**2)*(self.dotV_ref_2**2))
                    
        a_pu_2  =   (1/(self.u_ref_2**2))*(self.Deltap_ref_2-(self.dotV_ref_2**2)*a_pu_1)
        
        if self.check_valve==True:
            a_cv = - (1/(self.K_vs_cv*1000/60)**2) * (fluid.rho/1000) * (10**3)
            a_pu_1 = a_pu_1 + a_cv
        
        self.a_pu_1 = a_pu_1
        self.a_pu_2 = a_pu_2
        
        return a_pu_1, a_pu_2
    
    # Define function for showing properties
    def props(self):
        return dict(self.__dict__)
        
    def calc_Deltap(self, u, dotV):
        
        a_pu_1, a_pu_2 = self.hy_params()
        
        Deltap = a_pu_1 * (dotV**2) + a_pu_2 * (u**2)
        
        return Deltap
        
class controlvalve:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all
    rights reserved
    K_vs            [(m**3)/h]      Flow coefficient at full opening'''
    
    def __init__(self, K_vs=2.5):
        # set properties 
        self.K_vs   =   K_vs*(1000/60)
        
    def hy_params(self, fluid):
    
        a_va_1  = - (10**3) * (fluid.rho_SI/1000) * (1/(self.K_vs**2))
        
        self.a_va_1 = a_va_1
        
        return a_va_1
    
    def calc_Deltap(self,fluid,kappa,dotV):
        
        a_va_1 = self.hy_params(fluid)
        
        if kappa > 0:
            Deltap = a_va_1*(kappa**(-2))*(dotV**2)
        else:
            Deltap = float('-inf')
        
        return Deltap
    
        
        
                    
                    