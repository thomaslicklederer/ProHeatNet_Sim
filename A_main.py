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
import pprint
import pickle

import B_setup as su
import C_network as net
import E_hydraulic_model as hm
import F_thermal_model as tm
import G_visualization as vis
import H_toExcel as ex

################ delete old data #########################
get_ipython().magic('reset -f')
if os.name == 'nt':
    get_ipython().magic('cls')
else:
    get_ipython().magic('clear')
plt.close(fig='all')
##########################################################

# setup_name = "example"

# ############################ set-up ########################################
# path_config     =   "setup/example/configuration_example.py"
# path_param      =   "setup/example/parametrization_example.py" 
# path_scenario   =   "setup/example/scenario_example.py"

mysetup = su.setup(path_config, path_param, path_scenario)

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

############################ print solution in console #####################

pprint.pprint(vars(mysetup))
pprint.pprint(mysolutions)

############################################################################

############################ export results to logfiles #####################
# textfile
with open("results/results.txt", "w") as fout:
    fout.write("########## set-up ##########\n\n")
    
    fout.write(pprint.pformat(vars(mysetup)))
    
    fout.write("\n\n########## solutions ##########\n\n")

    fout.write(pprint.pformat(mysolutions))
    
    fout.write("\n\n##############################")
    
#pickle
myfile = open('results/results.pkl','wb')
myresults = {'PSM': mysetup.PSM, 'v_vec': mysetup.v_vec,
            'e_vec': mysetup.e_vec, 'topology': mysetup.topology,
            'scenario': mysetup.scenario, 'solutions': mysolutions}
pickle.dump(myresults, myfile)
myfile.close()   

#excel
ex.excelexport(setup_name)

############################################################################

############################ visualization #################################

# # graph
# vis.visualize_graph(mysetup, mygraph, mysolutions)

# # scenario
# vis.visualize_scenario(mysetup, mygraph, mysolutions, withnumbers=True)

# # hydraulic solution
# vis.visualize_hy_solution(mysetup, mygraph, mysolutions, withsymbols = False, withnumbers = True, data = 'both')

# # thermal solution
# vis.visualize_th_solution(mysetup, mygraph, mysolutions, th_model, withsymbols = False, withnumbers = True)

# # prosumer solutions
# vis.visualize_prosumer_results(mysetup, mygraph, mysolutions, th_model, withnumbers=True)

############################################################################

############################ plausability check ############################

sum_Q_trnsf = 0
sum_Q_trnsf2 = 0
sum_Q_loss = 0

for key in mysolutions['Q_trnsf'].keys():
    sum_Q_trnsf += mysolutions['Q_trnsf'][key]
for key in mysolutions['Q_trnsf2'].keys():
    sum_Q_trnsf2 += mysolutions['Q_trnsf2'][key]
for key in mysolutions['Q_loss'].keys():
    sum_Q_loss += mysolutions['Q_loss'][key]
    
print('sum_Q_trnsf:\t', float(sum_Q_trnsf))
print('sum_Q_trnsf2:\t', float(sum_Q_trnsf2))
print('sum_Q_loss:\t', float(sum_Q_loss))

############################################################################



