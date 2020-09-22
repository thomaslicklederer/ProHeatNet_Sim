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
import G_visualization as vis

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
try:
    th_model = tm.th_prob(mysetup, mygraph, mysolutions)
    th_model.problem_formulation()
    th_model.solve()
    th_model.do_some_calculations()
    mysolutions = th_model.save_solution()
except:
    pass

############################################################################

############################ print solution in console #####################

pprint(vars(mysetup))
pprint(mysolutions)


try:
    sum = 0
    sum2 = 0
    sum3 = 0

    for key in mysolutions['Q_trnsf'].keys():
        sum += mysolutions['Q_trnsf'][key]
        
    for key in mysolutions['Q_trnsf2'].keys():
        sum2 += mysolutions['Q_trnsf2'][key]

    for key in mysolutions['Q_loss'].keys():
        sum3 += abs(mysolutions['Q_loss'][key])
        
    print('sum_trnsf1: ', sum)
    print('sum_trnsf2: ', sum2)
    print('sum_loss: ', sum3)
except:
    pass

############################################################################

############################ visualization #################################

# graph
vis.visualize_graph(mysetup, mygraph, mysolutions)

# scenario
vis.visualize_scenario(mysetup, mygraph, mysolutions)

# hydraulic solution
vis.visualize_hy_solution(mysetup, mygraph, mysolutions, withsymbols = False, withnumbers = True)

# thermal solution
vis.visualize_th_solution(mysetup, mygraph, mysolutions, th_model, withsymbols = False, withnumbers = True)

############################################################################



