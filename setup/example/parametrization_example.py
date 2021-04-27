# ################################################################################
# create parameterized components

## fluid
###         rho             [kg/(m**3)]     density of fluid
###         cp              [J/(kg*K)]      specific isobaric heat capacity
###         mu              [Pa*s]          dynamic viscosity
myfluid         =   cm.fluid(rho=1000, cp=4200, mu=(1.0016*(10**(-3))))

## pumps
###         n_nom           [RPM = 1/min]   nominal pump speed
###         u_ref_1         [-]             controll variable ref. oper. state 1
###         dotV_ref_1      [l/min]         volume flow ref. oper. state 1
###         Deltap_ref_1    [hPa]           pressure difference ref. oper. state 1
###         u_ref_2         [-]             controll variable ref. oper. state 2
###         dotV_ref_2      [l/min]         volume flow ref. oper. state 2
###         Deltap_ref_2    [hPa]           pressure difference ref. oper. state 2
mypump1     =      cm.pump(n_nom=4100, u_ref_1=1, dotV_ref_1=0, Deltap_ref_1=402.21, 
                    u_ref_2=1, dotV_ref_2=55.33, Deltap_ref_2=0, check_valve = True, K_vs_cv = 6)
                    
## control valves
###         K_vs            [(m**3)/h]      Flow coefficient at full opening
myvalve1    =      cm.controlvalve(K_vs=2.5)

## heat exchangers
###         dotV_nom        [l/min]         nominal volume flow rate
###         Deltap_nom      [hPa]           head loss at nominal mass flow rate
###         k_nom           [W/((m*2)*K)]   nominal heat transfer coefficient
###         A               [m**2]          heat transfer surface
myhx1       =      cm.heatexchanger(dotV_nom=21.504, Deltap_nom=155, k_nom = 5270, A = 1.13)

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
mypipe1    =     cm.pipe(L=50,d_hy=0.022,epsilon=0.011,u_nom=0.22,zeta_instal=10,
                 N_layers = 2, d_layers = [0.022, 0.024, 0.062],
                 lambda_layers = [395, 0.04], h_ir=3500, h_or=float('inf'))
mypipe2    =      cm.pipe(L=10,d_hy=0.022,epsilon=0.011,u_nom=0.22,zeta_instal=3.5,
                 N_layers = 2, d_layers = [0.022, 0.024, 0.062],
                 lambda_layers = [395, 0.04], h_ir=3500, h_or=float('inf'))

# ################################################################################

# define fluid
self.myfluid    =   myfluid

# define components
self.components =      {}

## pipes in the subnetworks
### ('1h', '2h')
edge = ('1h', '2h')
self.components[edge] = {}
self.components[edge]['pipe'] = mypipe1

### ('2h', '3h')
edge = ('2h', '3h')
self.components[edge] = {}
self.components[edge]['pipe'] = mypipe1

### ('1c', '2c')
edge = ('1c', '2c')
self.components[edge] = {}
self.components[edge]['pipe'] = mypipe1

### ('2c', '3c')
edge = ('2c', '3c')
self.components[edge] = {}
self.components[edge]['pipe'] = mypipe1

# ### ('1h', '3h')
# edge = ('1h', '3h')
# self.components[edge] = {}
# self.components[edge]['pipe'] = mypipe1

# ### ('1c', '3c')
# edge = ('1c', '3c')
# self.components[edge] = {}
# self.components[edge]['pipe'] = mypipe1

## components in the edges with the substations (prosumers)
### ('1h', '1c')
edge = ('1h', '1c')
self.components[edge] = {}
self.components[edge]['heatexchanger'] = myhx1
self.components[edge]['pump'] = mypump1
self.components[edge]['controlvalve'] = myvalve1
self.components[edge]['pipe'] = mypipe2

### ('2h', '2c')
edge = ('2h', '2c')
self.components[edge] = {}
self.components[edge]['heatexchanger'] = myhx1
self.components[edge]['pump'] = mypump1
self.components[edge]['controlvalve'] = myvalve1
self.components[edge]['pipe'] = mypipe2

### ('3h', '3c')
edge = ('3h', '3c')
self.components[edge] = {}
self.components[edge]['heatexchanger'] = myhx1
self.components[edge]['pump'] = mypump1
self.components[edge]['controlvalve'] = myvalve1
self.components[edge]['pipe'] = mypipe2
