# ################################################################################
# create parameterized components

## fluid
###         rho             [kg/(m**3)]     density of fluid
###         cp              [J/(kg*K)]      specific isobaric heat capacity
###         mu              [Pa*s]          dynamic viscosity
myfluid         =   cm.fluid()

## pumps
###         n_nom           [RPM = 1/min]   nominal pump speed
###         u_ref_1         [-]             controll variable ref. oper. state 1
###         dotV_ref_1      [l/min]         volume flow ref. oper. state 1
###         Deltap_ref_1    [hPa]           pressure difference ref. oper. state 1
###         u_ref_2         [-]             controll variable ref. oper. state 2
###         dotV_ref_2      [l/min]         volume flow ref. oper. state 2
###         Deltap_ref_2    [hPa]           pressure difference ref. oper. state 2
mypump1     =      cm.pump()
                    
## control valves
###         K_vs            [(m**3)/h]      Flow coefficient at full opening
myvalve1    =      cm.controlvalve()

## heat exchangers
###         dotV_nom        [l/min]         nominal volume flow rate
###         Deltap_nom      [hPa]           head loss at nominal mass flow rate
###         k_nom           [W/((m*2)*K)]   nominal heat transfer coefficient
###         A               [m**2]          heat transfer surface
myhx1       =      cm.heatexchanger()

## pipes
###         L               [m]             pipe length
###         d_hy            [m]             (hydraulic) inner pipe diameter
###         epsilon         [mm]            roughness of inner pipe surface
###         u_nom           [m/s]           nominal flow velocity
###         zeta_instal     [-]             pressure loss coefficient (Druckverlustbeiwert)
###                                         for installations in pipe (such as bends, junctions etc.) 
###         N_layers        [-]             number of layers around pipe (for thermal losses),
###                                         e.g. pipe wall (1), insulation (2), soil environment (3)
###         d_layers        [m]             diameters [m] at the borders of the layers
###         lambda_layers   [W/(m*K)]       thermal conductivity of layers
###         h_ir            [W/((m**2)*K)]  heat transfer coefficient at inside of pipe
###         h_or            [W/((m**2)*K)]  heat transfer coefficient at outside of outest pipe layer
mypipe1    =      cm.pipe()

# ################################################################################

# define components for each edge automatically
self.components =      {}
for edge in e_vec:
    self.components[edge] = {}
    self.components[edge]['pipe'] = mypipe1
    # filter for connections that connect warm and cold subnetwork and add 
    # devices
    if (('h' in edge[0] and 'c' in edge[1]) or ('h' in edge[1] and 'c' in 
                                                edge[0])):
        self.components[edge]['heatexchanger'] = myhx1
        self.components[edge]['pump'] = mypump1
        self.components[edge]['controlvalve'] = myvalve1
        self.components[edge]['pipe'] = mypipe1
