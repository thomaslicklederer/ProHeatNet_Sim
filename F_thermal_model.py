import numpy as np
import scipy as sp
from pprint import pprint
import gurobipy as gp
from gurobipy import GRB

class th_prob:
    '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all 
    rights reserved
    DESCRIPTION TO BE ADDED'''
    
    # Initializer / Instance Attributes
    def __init__(self, setup, graph, solution_dict):
    
        self.setup  =   setup
        self.graph  =   graph
        self.solution_dict = solution_dict
        
    def problem_formulation(self):
        
        self.preparation()
        
        # create a new gurobi model
        self.gumo = gp.Model("thermal_model")
        
        # create variables
        T = self.gumo.addMVar(shape= len(self.general_stuff['T_vec']), obj=1, vtype=GRB.CONTINUOUS,
                                name="T", lb = 275, ub = 370)
                                
        # DeltaT = gumo.addMVar(shape= len(self.general_stuff.DeltaT_vec), obj=1, vtype=GRB.CONTINUOUS,
                                # name="DeltaT", lb = -100, ub = 100)
                                
        # create objective function: minimize pressure losses
        self.gumo.setObjective(1)                        # just solve
        
        # add constraints
        self.gumo.addConstr(self.M_mat @ T == 0)
        self.gumo.addConstr(self.L_mat @ T == 0)
        self.gumo.addConstr(self.E_mat @ T == self.R_vec)
        self.gumo.addConstr(self.S_T_in_mat @ T == self.T_sec_in_vec)
        
    def solve(self):
        
        self.gumo.optimize()
        
        self.T_solvec = self.gumo.X
                
    def do_some_calculations(self):
        
        self.DeltaT_solvec = self.C_mat * self.T_solvec

        self.Q_trnsf = {}
        self.Q_loss = {}
        self.Q_trnsf2 = {}
        
        index = []
        for edge in self.setup.e_hc_vec:
            V = self.solution_dict['dotV'][edge]
            index = []
            for i in range(len(self.general_stuff['dotV_vec'])):
                if self.general_stuff['dotV_vec'][i] == edge:
                    index.append(i)
            if len(index) > 1:
                raise ValueError('too many indices')
            dT = self.DeltaT_solvec[index]
            self.Q_trnsf[edge] = self.setup.myfluid.rho * self.setup.myfluid.cp * 1/60 * V * dT
        
        index = []
        for edge in self.setup.PSM:
            V = self.solution_dict['dotV'][edge]
            index = []
            for i in range(len(self.general_stuff['dotV_vec'])):
                if self.general_stuff['dotV_vec'][i] == edge:
                    index.append(i)
            if len(index) > 1:
                raise ValueError('too many indices')
            dT = self.DeltaT_solvec[index]
            self.Q_trnsf2[edge] = self.setup.myfluid.rho * self.setup.myfluid.cp * 1/60 * V * dT
            
        index = []
        for edge in self.setup.e_h_vec:
            V = self.solution_dict['dotV'][edge]
            index = []
            for i in range(len(self.general_stuff['dotV_vec'])):
                if self.general_stuff['dotV_vec'][i] == edge:
                    index.append(i)
            if len(index) > 1:
                raise ValueError('too many indices')
            dT = self.DeltaT_solvec[index]
            self.Q_loss[edge] = self.setup.myfluid.rho * self.setup.myfluid.cp * 1/60 * V * dT
            
        index = []
        for edge in self.setup.e_c_vec:
            V = self.solution_dict['dotV'][edge]
            index = []
            for i in range(len(self.general_stuff['dotV_vec'])):
                if self.general_stuff['dotV_vec'][i] == edge:
                    index.append(i)
            if len(index) > 1:
                raise ValueError('too many indices')
            dT = self.DeltaT_solvec[index]
            self.Q_loss[edge] = self.setup.myfluid.rho * self.setup.myfluid.cp * 1/60 * V * dT

        
            
    def save_solution(self):
        # Edges_Temps_dict = self.general_stuff['Edges_Temps_dict']
        
        self.solution_dict['Q_trnsf'] = self.Q_trnsf
        self.solution_dict['Q_trnsf2'] = self.Q_trnsf2
        self.solution_dict['Q_loss'] = self.Q_loss

        for i in range(len(self.general_stuff['T_vec'])):
            name = self.general_stuff['T_vec'][i]
            value = self.T_solvec[i]-273.15
            self.solution_dict['T'][name] = value
            
        for i in range(len(self.general_stuff['dotV_vec'])):
            name = self.general_stuff['dotV_vec'][i]
            value = self.DeltaT_solvec[i]
            self.solution_dict['DeltaT'][name] = value
            
        return self.solution_dict
        

    
    def preparation(self):
    
        self.general_stuff()
        
        self.M_mat = self.capacity_flow_mat()
                
        self.L_mat = self.same_T_out()
        
        self.E_mat, self.R_vec = self.th_comp_mods()
        
        self.C_mat = self.temps_to_diffs()
    
        self.S_T_in_mat, self.T_sec_in_vec = self.set_PSM_inlet_temps()
        
        self.dotV_vec = self.set_PSM_flows()
    
        pprint(self.general_stuff)
        
    def general_stuff(self):
        
        self.general_stuff = {}
        dotV_dict = self.solution_dict['dotV']
        
        import networkx as nx
        
        # incidence_matrix
        B_mat = nx.linalg.graphmatrix.incidence_matrix(self.graph.G, nodelist=self.setup.v_vec,
                                                       edgelist = self.setup.e_vec, 
                                                       oriented=True)
        self.B_mat = B_mat
        
        # create vector with entries for each end of an edge (edge-ende-vector) 
        # at the same time create vector that determines the signum for temperatures
        # based on if they are at the beginning (-1) or end (1) of an edge
        # at the same time create vectors that represent the association between the
        # temperatures and the edges / nodes
        ee_vec = []
        T_sign = []
        T_e_assoc = []
        T_v_assoc = []
        T_PSM_assoc = []
        for node_row in range(B_mat.shape[0]):
            for edge_col in range(B_mat.shape[1]):
                if B_mat[node_row, edge_col] == -1:
                    ee_vec.append(self.setup.e_vec[edge_col])
                    T_sign.append(-1)
                    T_e_assoc.append(edge_col)
                    T_v_assoc.append(node_row)
                elif B_mat[node_row, edge_col] == 1:
                    ee_vec.append((self.setup.e_vec[edge_col][1],self.setup.e_vec[edge_col][0]))
                    T_sign.append(1)
                    T_e_assoc.append(edge_col)
                    T_v_assoc.append(node_row)
                elif B_mat[node_row, edge_col] == 0:
                    pass
                else:
                    raise ValueError('Entries of B_mat must either be 0 or 1 or -1!')
        
        for node_nbr in T_v_assoc:
            T_PSM_assoc.append(int(self.setup.v_vec[node_nbr][:-1]))
        

        # add temperatures for secondary side
        T_sec_vec = []
        T_vec = [i for i in ee_vec]
        for i in self.setup.PSM:
            T_sec_vec.append(''.join(('PSM', str(i), 'h')))
            T_vec.append(''.join(('PSM', str(i), 'h')))
            T_v_assoc.append(None)
            T_sign.append(1)
            T_PSM_assoc.append(i)
            T_e_assoc.append(''.join(('PSM', str(i))))
            T_sec_vec.append(''.join(('PSM', str(i), 'c')))
            T_vec.append(''.join(('PSM', str(i), 'c')))
            T_v_assoc.append(None)
            T_sign.append(-1)
            T_PSM_assoc.append(i)
            T_e_assoc.append(''.join(('PSM', str(i))))
                        
        
        temp = []
        for i in range(len(T_vec)):
            if T_v_assoc[i] != None:
                temp.append([self.setup.v_vec[T_v_assoc[i]], self.setup.e_vec[T_e_assoc[i]], T_sign[i], T_PSM_assoc[i]])
            else:
                temp.append([[],T_e_assoc[i],T_sign[i],T_PSM_assoc[i]])
                

        T_dict = dict(zip(T_vec, temp))
               
        
        # create DeltaT_vec
        DeltaT_vec = [entry for entry in self.setup.e_vec]
        for i in self.setup.PSM:
            DeltaT_vec.append(''.join(('PSM', str(i))))
            
        # mu_dict
        mu_dict = {}
        for i in self.setup.PSM:
            PSM = i
            mu = self.setup.scenario[i]['mu']
            mu_dict[PSM] = mu
            
        # dotV_sec_dict
        dotV_sec_dict = {}
        for i in self.setup.PSM:
            PSM = i
            V = self.setup.scenario[i]['dotV_sec_in']
            dotV_sec_dict[PSM] = V
            
        # Edges_Temps_dict
        Edges_Temps_dict={}
        for prosumer in list(mu_dict.keys()):
            Edges_Temps_dict[prosumer]={}
            for i in range(len(T_vec)):
                temp = T_vec[i]
                if (str(prosumer) in temp[0]) and (str(prosumer) in temp[1]):
                    if 'h' in temp[0]:
                        if mu_dict[prosumer] == -1:
                            Edges_Temps_dict[prosumer]['T_in_HTNW'] = temp
                        if mu_dict[prosumer] == 1:
                            Edges_Temps_dict[prosumer]['T_out_HTNW'] = temp
                    if 'c' in temp[0]:
                        if mu_dict[prosumer] == 1:
                            Edges_Temps_dict[prosumer]['T_in_HTNW'] = temp
                        if mu_dict[prosumer] == -1:
                            Edges_Temps_dict[prosumer]['T_out_HTNW'] = temp
                    try:
                        Edges_Temps_dict[prosumer]['dotV_HTNW'] = dotV_dict[temp]
                        Edges_Temps_dict[prosumer]['edge'] = temp
                    except KeyError:
                        pass
                                
                if (str(prosumer) in temp):                    
                    if 'h' in temp:
                        if mu_dict[prosumer] == 1:
                            Edges_Temps_dict[prosumer]['T_in_PSM'] = temp
                        if mu_dict[prosumer] == -1:
                            Edges_Temps_dict[prosumer]['T_out_PSM'] = temp
                    if 'c' in temp:                                
                        if mu_dict[prosumer] == -1:
                            Edges_Temps_dict[prosumer]['T_in_PSM'] = temp
                        if mu_dict[prosumer] == 1:
                            Edges_Temps_dict[prosumer]['T_out_PSM'] = temp
                    Edges_Temps_dict[prosumer]['dotV_PSM'] = dotV_sec_dict[prosumer]

        for temp_diff in DeltaT_vec:
            if ('h' in temp_diff[0] and 'h' in temp_diff[1]) or ('c' in temp_diff[0] and 'c' in temp_diff [1]):
                Edges_Temps_dict[temp_diff]={}
                signum = np.sign(dotV_dict[temp_diff])
                if signum == 1:
                    Edges_Temps_dict[temp_diff]['T_in'] = temp_diff
                    Edges_Temps_dict[temp_diff]['T_out'] = (temp_diff[1], temp_diff[0])
                elif signum == -1:
                    Edges_Temps_dict[temp_diff]['T_in'] = (temp_diff[1], temp_diff[0])
                    Edges_Temps_dict[temp_diff]['T_out'] = temp_diff
                Edges_Temps_dict[temp_diff]['dotV'] = dotV_dict[temp_diff]
                Edges_Temps_dict[temp_diff]['edge'] = temp_diff

        
        # save stuff
        self.general_stuff['B_mat'] = B_mat
        self.general_stuff['ee_vec'] = ee_vec
        self.general_stuff['T_sign'] = T_sign
        self.general_stuff['T_e_assoc'] = T_e_assoc
        self.general_stuff['T_v_assoc'] = T_v_assoc
        self.general_stuff['T_PSM_assoc'] = T_PSM_assoc
        self.general_stuff['T_sec_vec'] = T_sec_vec
        self.general_stuff['T_vec'] = T_vec
        self.general_stuff['T_dict'] = T_dict
        self.general_stuff['DeltaT_vec'] = DeltaT_vec
        self.general_stuff['mu_dict'] = mu_dict
        self.general_stuff['dotV_sec_dict'] = dotV_sec_dict
        self.general_stuff['Edges_Temps_dict'] = Edges_Temps_dict
                        
    
    def capacity_flow_mat(self):
        '''create M_mat: matrix for balance of capacity flows in/out nodes'''
                
        ee_vec = self.general_stuff['ee_vec']
        T_dict = self.general_stuff['T_dict']
        T_sec_vec = self.general_stuff['T_sec_vec']
        rho = self.setup.myfluid.rho
        dotV_dict = self.solution_dict['dotV']
        
        # Mixing in Nodes: Calculate Matrix for Energy Balance in nodes
        M_mat = np.zeros((len(self.setup.v_vec),len(ee_vec)))
        for j in range(len(ee_vec)):
            for i in range(len(self.setup.v_vec)):
                if T_dict[ee_vec[j]][0]==self.setup.v_vec[i]:
                    temp_edge = T_dict[ee_vec[j]][1]
                    M_mat[i,j]=dotV_dict[temp_edge]*T_dict[ee_vec[j]][2]*rho
        # concatenate submatrix with zeros for temperatures on secondary prosumer side
        M_mat = np.concatenate((M_mat, np.zeros((len(self.setup.v_vec),len(T_sec_vec)))),axis=1)
        
        M_mat = sp.sparse.csc_matrix(M_mat)
        
        return M_mat
        
    def same_T_out(self):
        '''create L_mat: matrix for all outflowing temperatures of node have to be the same'''
        M_mat = self.M_mat.todense()
        
        M_mat_sign = np.where(M_mat<0, -1, 0)
        outgoing_fluxes = np.count_nonzero(M_mat_sign < 0, axis=1)
        relev_nodes = list(np.argwhere(outgoing_fluxes>1))
        relev_nodes = [int(i) for i in relev_nodes]
        M_mat_sign_tilde = M_mat_sign[relev_nodes,:]
        relev_Ts = []
        for row in range(M_mat_sign_tilde.shape[0]):
            temp=list(np.argwhere(M_mat_sign_tilde[row]==-1))
            temp = [int(i) for i in temp]
            relev_Ts.append(temp)
        L_mat = []
        for i in range(len(relev_Ts)):
            node_Ts = relev_Ts[i]
            for j in range(1, len(node_Ts)):
                L_mat.append(list(np.zeros(M_mat_sign_tilde.shape[1])))
                L_mat[i+j-1][node_Ts[0]]=-1
                L_mat[i+j-1][node_Ts[j]]=1
        
        if L_mat == []:
            L_mat = np.zeros((1,np.shape(M_mat)[1]))
        
        L_mat = sp.sparse.csc_matrix(L_mat)
        
        return L_mat
        
    def th_comp_mods(self):
        '''create E_mat and R_vec: matrix and right hand side for summary of thermal component models'''
        
        e_vec = self.setup.e_vec
        mu_dict = self.general_stuff['mu_dict']
        T_vec = self.general_stuff['T_vec']
        dotV_dict = self.solution_dict['dotV']
        dotV_sec_dict = self.general_stuff['dotV_sec_dict']
        DeltaT_vec = self.general_stuff['DeltaT_vec']
        T_dict = self.general_stuff['T_dict']
        Edges_Temps_dict = self.general_stuff['Edges_Temps_dict']
        
        R_vec = np.zeros(len(e_vec)+len(list(mu_dict.keys())))

        ## E1_mat for the interconnecting edges
        E1_mat = np.zeros((2*len(list(mu_dict.keys())), len(list(T_dict.keys()))))
        mycounter = -2
        for prosumer in list(mu_dict.keys()):
            mycounter += 2
            mu = mu_dict[prosumer]
            dotV_HTNW = Edges_Temps_dict[prosumer]['dotV_HTNW']
            dotV_PSM = Edges_Temps_dict[prosumer]['dotV_PSM']
            T_in_HTNW = Edges_Temps_dict[prosumer]['T_in_HTNW']
            T_in_PSM = Edges_Temps_dict[prosumer]['T_in_PSM']
            T_out_HTNW = Edges_Temps_dict[prosumer]['T_out_HTNW']
            T_out_PSM = Edges_Temps_dict[prosumer]['T_out_PSM']
            params_HTNW, params_PSM = self.setup.components[Edges_Temps_dict[prosumer]['edge']]\
                ['heatexchanger'].th_params(mu, dotV_HTNW, dotV_PSM, self.setup.myfluid)

 
            for i in range(len(T_vec)):
                curTemp = T_vec[i]
                if curTemp == T_in_HTNW:
                    E1_mat[mycounter, i] = params_HTNW['T_in_HTNW']
                    E1_mat[mycounter+1, i] = params_PSM['T_in_HTNW']
                elif curTemp == T_out_HTNW:
                    E1_mat[mycounter, i] = params_HTNW['T_out_HTNW']
                    E1_mat[mycounter+1, i] = params_PSM['T_out_HTNW']
                elif curTemp == T_in_PSM:
                    E1_mat[mycounter, i] = params_HTNW['T_in_PSM']
                    E1_mat[mycounter+1, i] = params_PSM['T_in_PSM']
                elif curTemp == T_out_PSM:
                    E1_mat[mycounter, i] = params_HTNW['T_out_PSM']
                    E1_mat[mycounter+1, i] = params_PSM['T_out_PSM']
            R_vec[mycounter] = 0
            R_vec[mycounter+1] = 0
            

        ## E2_mat for the edges within the subnetworks
        E2_mat = np.zeros((len(e_vec)-len(list(mu_dict.keys())), len(list(T_dict.keys()))))
        mycounter = -1

        for j in range(len(DeltaT_vec)):
            temp_diff = DeltaT_vec[j]
            if ('h' in temp_diff[0] and 'h' in temp_diff[1]) or ('c' in temp_diff[0] and 'c' in temp_diff [1]):
                mycounter += 1                
                try:
                    T_in = Edges_Temps_dict[temp_diff]['T_in']
                    T_out = Edges_Temps_dict[temp_diff]['T_out']
                except KeyError:
                    print('fail')
                
                dotV = Edges_Temps_dict[temp_diff]['dotV']
                [b_pi_1, b_pi_2, b_pi_3] = self.setup.components[Edges_Temps_dict[temp_diff]['edge']]\
                    ['pipe'].th_params(self.setup.myfluid, dotV, T_soil = self.setup.scenario['T_soil'])
                # b_pi1, b_pi2, b_pi3 = self.networkgraph.G.edges[Edges_Temps_dict[temp_diff]['edge']]\
                    # ['pipe'].th_paras(self.fluid, dotV)
                for i in range(len(T_vec)):
                    curTemp = T_vec[i]
                    if curTemp == T_in:
                        E2_mat[mycounter, i] = b_pi_1
                    if curTemp == T_out:
                        E2_mat[mycounter, i] = b_pi_2
                R_vec[np.shape(E1_mat)[0]+mycounter] = b_pi_3
   
        E_mat = np.concatenate((E1_mat, E2_mat), axis=0)
        E_mat = sp.sparse.csc_matrix(E_mat)
        E_mat = E_mat
        R_vec = R_vec
        
        return E_mat, R_vec
        
    def temps_to_diffs(self):
        '''create C_mat: relation between temperature differences and temperatures'''
        
        T_vec = self.general_stuff['T_vec']
        DeltaT_vec = self.general_stuff['DeltaT_vec']
        T_dict = self.general_stuff['T_dict']
        
        C_mat = np.zeros((len(DeltaT_vec), len(T_vec)))
        for diff_row in range(len(DeltaT_vec)):
            for temp_column in range(len(T_vec)):
                if T_dict[T_vec[temp_column]][1] == DeltaT_vec[diff_row]:
                    C_mat[diff_row,temp_column] = T_dict[T_vec[temp_column]][2]
        C_mat = sp.sparse.csc_matrix(C_mat)
        
        return C_mat
        
    def set_PSM_inlet_temps(self):
        '''create S_T_in_mat and T_sec_in_vec: matrix and vector to determine the inlet temperatures at the secondary side'''
    
        T_vec = self.general_stuff['T_vec']
        T_sec_vec = self.general_stuff['T_sec_vec']
        Edges_Temps_dict = self.general_stuff['Edges_Temps_dict']
        T_sec_in = []
        index = []
        
        for P in self.setup.PSM:
            name_in = Edges_Temps_dict[P]['T_in_PSM']
            value_in = self.setup.scenario[P]['T_sec_in']
            T_sec_in.append(value_in)
            for i in range(len(T_vec)):
                temp = T_vec[i]
                if temp == name_in:
                    index.append(i)
        
        S_T_in_mat = np.zeros((len(T_sec_in), len(T_vec)))
        
        for j in range(len(index)):
            k = index[j]
            S_T_in_mat[j,k]=1
        
        T_sec_in_vec = np.array(T_sec_in)
        
        return S_T_in_mat, T_sec_in_vec

    def set_PSM_flows(self):
        '''create S_V_mat and dotV_sec_vec: matrix and vector to determine the flows at the secondary side'''
        
        self.general_stuff['dotV_vec'] = self.setup.e_vec + self.setup.PSM
        name_vec = self.general_stuff['dotV_vec']
        
        dotV_vec = []
        
        for i in range(len(name_vec)):
            name = name_vec[i]
            dotV_vec.append(self.solution_dict['dotV'][name])
        
        return dotV_vec