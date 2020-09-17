import math
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from IPython import get_ipython
import os
import sympy as sym
import gurobipy as gp
from gurobipy import GRB
from pprint import pprint

import B_setup as su
import C_network as net
import E_hydraulic_model as hm
import F_thermal_model as tm

################ delete old data #########################
get_ipython().magic('reset -f')
if os.name == 'nt':
    get_ipython().magic('cls')
else:
    get_ipython().magic('clear')
plt.close(fig='all')
##########################################################

############################ set-up ########################################
path_config     =   "set-up/BA_configuration.py"
path_param      =   "set-up/BB_parametrization.py" 
path_scenario   =   "set-up/BC_scenario.py"
dicenbr         =   7

mysetup = su.setup(path_config, path_param, path_scenario, dicenbr)

############################################################################

############################ prepare solutions container ###################

mysolutions = mysetup.setup_solutions()

############################################################################

############################ network representation ########################

mygraph = net.networkgraph(mysetup)

############################################################################

############################ hydraulic model ###############################

hy_model = hm.hy_prob(mysetup, mygraph, mysolutions)
hy_model.problem_formulation()
hy_model.solve_group()
mysolutions = hy_model.save_solution()

############################################################################

############################ thermal model #################################

th_model = tm.th_prob(mysetup, mygraph, mysolutions)
th_model.problem_formulation()
th_model.solve()
th_model.do_some_calculations()
mysolutions = th_model.save_solution()

############################################################################
pprint(vars(mysetup))
pprint(mysolutions)

sum = 0

for key in mysolutions['Q_trnsf'].keys():
    sum += mysolutions['Q_trnsf'][key]
    
print('sum: ', sum)





