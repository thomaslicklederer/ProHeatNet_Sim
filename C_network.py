import math
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpltlib
import networkx as nx
import pandas as pd
from IPython import get_ipython
import os
import sympy as sym
import matplotlib.image as mpimg
from pprint import pprint
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

import D_component_models as cm

class networkgraph:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all 
    rights reserved
    DESCRIPTION TO BE ADDED'''
    
    # Initializer / Instance Attributes
    def __init__(self,setup):
                
        # create directed graph
        G = nx.DiGraph()
        G.add_nodes_from(setup.v_vec)
        G.add_edges_from(setup.e_vec)
        G_prim = nx.Graph(G)
        self.G_prim = G_prim
        self.G = G
    
    # Define function for showing properties
    def props(self):
        print('\n%%%%%%%%%%%%%%%%% NETWORKGRAPH-PROPERTIES %%%%%%%%%%%%%%%%%')
        pprint(vars(self))
        print('%%%%%%%%%%%%%%%%% NETWORKGRAPH-PROPERTIES %%%%%%%%%%%%%%%%%\n')
    