class setup:
    """Code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, 2020,
    all rights reserved
    DESCRIPTION TO BE ADDED"""
    
    def __init__(self, path_config, path_param, path_scenario, dicenbr):
        import D_component_models as cm
        import math
        import random
        import numpy as np
        from pprint import pprint
            
        exec(open(path_config).read(), locals())
        
        #### --------------------------------------------------------------------
        # sort topology
        self.topology    =   sorted(self.topology, key=lambda tup: (tup[0],tup[1]) )
        
        # create prosumer list
        PSM = list(range(1,self.M+1,1))
        self.PSM = PSM
        
        # create sorted nodes vector
        temp1 = ('h','c')
        temp2 = [str(elem) + temp1[0] for elem in self.PSM]
        temp3 = [str(elem) + temp1[1] for elem in self.PSM]
        v_vec = temp2 + temp3
        del temp1, temp2, temp3
        self.v_vec = v_vec

        # create sorted edges vector
        temp1 = [(str(elem)+'h', str(elem)+'c') for elem in self.PSM]
        temp2 = [(str(elem[0])+'h', str(elem[1])+'h') for elem in self.topology]
        temp3 = [(str(elem[0])+'c', str(elem[1])+'c') for elem in self.topology]
        e_vec = temp2 + temp1 + temp3
        del temp1, temp2, temp3
        self.e_vec = e_vec
        
        # edges within hot subnetwork
        e_h_vec = []
        # edges within cold subnetwork
        e_c_vec = []
        # edges interconnecting the subnetworks
        e_hc_vec = []
        
        for edge in self.e_vec:
            if (('h' in edge[0] and 'c' in edge[1]) or ('h' in edge[1] and 'c' in 
                                                        edge[0])):
                e_hc_vec.append(edge)
            if 'h' in edge[0] and 'h' in edge[1]:
                e_h_vec.append(edge)
            if 'c' in edge[0] and 'c' in edge[1]:
                e_c_vec.append(edge)
                
        self.e_h_vec    =   e_h_vec
        self.e_hc_vec    =   e_hc_vec
        self.e_c_vec    =   e_c_vec
        
        # edge_prosumer_vector
        e_PSM_dict = {}
        for i in range(len(self.e_hc_vec)):
            key = self.e_hc_vec[i]
            v = key[0]
            if 'h' in v:
                temp = 'h'
            if 'c' in v:
                temp = 'c'
            P = int(str(v.replace(temp, '')))
            value = P
            e_PSM_dict[key] = value
        self.e_PSM_dict = e_PSM_dict
        #### --------------------------------------------------------------------
                
        exec(open(path_param).read(), locals())
        exec(open(path_scenario).read(), locals())
               
        print('\n%%%%%%%%%%%%%%%%% SETUP-BEGIN %%%%%%%%%%%%%%%%%')
        pprint(vars(self))
        print('%%%%%%%%%%%%%%%%% SETUP-END %%%%%%%%%%%%%%%%%\n')
        
    def setup_solutions(self):
        
        # volume flows        
        dotV = {}
        
        # pressure differences
        Deltap = {}
        
        # temperatures
        T = {}
        
        # temperature differences
        DeltaT = {}
        
        # transferred heat for prosumers
        Q_trnsf = {}
        
        # heat losses in pipes
        Q_loss = {}
        
        solutions = {}
        solutions['dotV']   = dotV
        solutions['Deltap'] = Deltap
        solutions['T'] = T
        solutions['DeltaT'] = DeltaT
        solutions['Q_trnsf'] = Q_trnsf
        solutions['Q_loss'] = Q_loss
        
        for i in self.PSM:
            solutions['dotV'][i] = self.scenario[i]['dotV_sec_in']
        
        return solutions
