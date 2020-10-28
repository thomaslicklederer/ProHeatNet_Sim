import numpy as np
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import itertools
import networkx as nx
import scipy as sp


class hy_prob:
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
        
        # create optimization problems
        self.problem_group = []
        
        for cntr in range(np.shape(self.sigma_mat)[1]):
            self.problem_group.append(self.single_prob(cntr))
            
    def solve_group(self):
       
        solutions={}
        fails = 0
        objVal_vec = []
        varSol_vec = []
        cases_vec = []
        iterations_vec = []
        
        for cntr in range(np.shape(self.sigma_mat)[1]):
            
            m = self.problem_group[cntr]
            
            m.optimize()
            
            case="".join(("m",str(cntr+1)))
            solutions[case]={}
            if m.Status == 2:
                solutions[case]["objVal"]=m.objVal
                for var in m.getVars():
                    solutions[case][var.VarName]=var.X
                objVal_vec.append(m.objVal)
                varSol_vec.append(m.X)
                iterations_vec.append(cntr)
                cases_vec.append(case)

            else:
                fails += 1
                solutions[case]["objVal"]=-9999
                for var in m.getVars():
                    solutions[case][var.VarName]=-9999
            
        
        # find best solution
        bestobj = min(objVal_vec)
        
        bestcase = cases_vec[objVal_vec.index(bestobj)]
        bestind = iterations_vec[objVal_vec.index(bestobj)]
        
        # recreate and recalculate the model with the best solution
        case="best"
        solutions[case]={}
            
        m = self.single_prob(bestind)
        m.optimize()
                
        if m.Status == 2:
                solutions[case]["objVal"]=m.objVal
                for var in m.getVars():
                    solutions[case][var.VarName]=var.X
                objVal_vec.append(m.objVal)
                varSol_vec.append(m.X)
                iterations_vec.append(cntr)
                cases_vec.append(case)
            
        self.dotV_solvec = m.X[:int(len(m.X)/2)]
        self.Deltap_solvec = m.X[int(len(m.X)/2):]
    
    def save_solution(self):
    
        for i in range(len(self.setup.e_vec)):
            edge = self.setup.e_vec[i]
            dotV = self.dotV_solvec[i]
            Deltap = self.Deltap_solvec[i]
            
            self.solution_dict['dotV'][edge]=dotV
            self.solution_dict['Deltap'][edge]=Deltap

        
        return self.solution_dict

    def preparation(self):
        
        # calculate matrices
        self.calc_matrices()
        
        # calculate parameter vectors
        self.create_param_vecs()
        
        
    def single_prob(self,cntr):

        
        # create a new gurobi model
        gum = gp.Model("hydraulic_model")
        gum.params.LogToConsole = 0
        
        # create variables
        # dotV
        dotV = gum.addMVar(shape= len(self.setup.e_vec), obj=1, vtype=GRB.CONTINUOUS,
                                name="dotV", lb = -50, ub = 50)
        
        # Deltap
        Deltap = gum.addMVar(shape= len(self.setup.e_vec), obj=1, vtype=GRB.CONTINUOUS,
                                name="dotV", lb = -5000, ub = 5000)
        
        # create objective function: minimize pressure losses
        gum.setObjective(Deltap @ self.Q @ Deltap, GRB.MINIMIZE)
        
        # add constraints
        ## 1. Kirchhoff law
        gum.addConstr(self.B_mat @ dotV == 0, "K1")
        
        ## 2. Kirchhoff law
        gum.addConstr(self.O_mat @ Deltap == 0, "K2")
        
        ## edge constraints
        ### edges in hot and cold subnetworks
        for i in range(len(self.setup.e_vec)):
            edge = self.setup.e_vec[i]
            if (edge in self.setup.e_h_vec) or (edge in self.setup.e_c_vec):
                a_pi_1 = self.setup.components[edge]['pipe'].hy_params(self.setup.myfluid)
                r_e = a_pi_1
                y_e = 0
                gum.addConstr(Deltap[i] == self.sigma_mat[i,cntr]*(r_e * (dotV[i]@dotV[i]) + y_e))
                gum.addConstr(-self.sigma_mat[i,cntr]*dotV[i] <= 0)
        
        ### edges that interconnect the subnetworks
        for i in range(len(self.setup.e_vec)):
            edge = self.setup.e_vec[i]
            if (edge in self.setup.e_hc_vec):
                P = self.setup.e_PSM_dict[edge]
                Pind = self.setup.PSM.index(P)
                if self.sigma_mat[i,cntr]*self.mu_vec[Pind]>0:
                    raise ValueError('sigma_mat[i,cntr]*mu_vec[Pind]>0')
                try:
                    a_pi_1 = self.setup.components[edge]['pipe'].hy_params(self.setup.myfluid)
                except:
                    a_pi_1 = 0
                a_hx_1 = self.setup.components[edge]['heatexchanger'].hy_params()
                a_pu_1, a_pu_2 = self.setup.components[edge]['pump'].hy_params()
                a_va_1 = self.setup.components[edge]['controlvalve'].hy_params(self.setup.myfluid)
                gamma_minus = 0.5*(1-self.mu_vec[Pind])
                gamma_plus = 0.5*(1+self.mu_vec[Pind])
                if (self.kappa_vec[Pind]<=(10**(-2))) and (gamma_minus != 0):
                    r_e = 99999999
                elif gamma_minus == 0:
                    r_e = a_pi_1 + a_hx_1 + gamma_plus*(a_pu_1) + gamma_minus * (a_va_1)
                else:
                    r_e = a_pi_1 + a_hx_1 + gamma_plus*(a_pu_1) + gamma_minus * (a_va_1) * (self.kappa_vec[Pind]**(-2))
                y_e = gamma_plus*(a_pu_2)*(self.u_vec[Pind]**2)
                gum.addConstr(Deltap[i] == self.sigma_mat[i,cntr]*(r_e * (dotV[i]@dotV[i]) + y_e))
                gum.addConstr(-self.sigma_mat[i,cntr]*dotV[i] <= 0)

        # Adjust parameters
        gum.params.NonConvex = 2
        
        
        return gum
       
    def create_param_vecs(self):

        
        # pi_P
        self.pi_vec = [self.setup.scenario[P]['pi'] for P in self.setup.PSM]
            
        # mu_P
        self.mu_vec = [self.setup.scenario[P]['mu'] for P in self.setup.PSM]
        
        # kappa_P
        self.kappa_vec = [self.setup.scenario[P]['kappa'] for P in self.setup.PSM]
        
        # u_P
        self.u_vec = [self.setup.scenario[P]['u'] for P in self.setup.PSM]
        
        # sigma_e
        ## number of edges with variable volume flow directions
        var_edges       = self.setup.e_h_vec + self.setup.e_c_vec
        nbr_var_edges   = len(var_edges)
        
        ## generate different combinations of flow directions
        mycombs = [p for p in itertools.product((-1,+1), repeat=int(nbr_var_edges/2))]
        
        ## initiate sigma_mat
        sigma_mat = np.zeros((len(self.setup.e_vec),len(mycombs)),dtype=int)
        
        ## find pairs of connections in the subnetworks
        pairs = []
        for i in range(len(self.setup.e_vec)):
            edge = self.setup.e_vec[i]
            if ('h' in edge[0] and 'h' in edge[1]): 
                a = i
                temp1 = edge[0].replace('h','c')
                temp2 = edge[1].replace('h','c')
                b = self.setup.e_vec.index((temp1,temp2))
                pairs.append((a,b))
                
        ## set sigma_mat values for interconnecting edges
        for j in range(np.shape(sigma_mat)[1]):
            for i in range(np.shape(sigma_mat)[0]):
                if not (self.setup.e_vec[i] in var_edges):
                    e = self.setup.e_vec[i]
                    P = self.setup.e_PSM_dict[e]
                    for k in range(len(self.setup.PSM)):
                        if P == self.setup.PSM[k]:
                            mu = self.mu_vec[k]
                    if mu == -1: # consumer
                        sigma_mat[i,j] = 1
                    elif mu == 1: # producer
                        sigma_mat[i,j] = -1
        
        ## set sigma_mat values for edges with variable flow directions
        for j in range(np.shape(sigma_mat)[1]):
            combination = mycombs[j]
            for k in range(len(combination)):
                pair = pairs[k]
                
                sigma_mat[pair[0],j]=combination[k]
                sigma_mat[pair[1],j]=-combination[k]
                
        self.sigma_mat = sigma_mat
        
        ## Q matrix for objective function
        self.Q = np.diag((10**(-4))*np.ones(len(self.setup.e_vec)))
        
        
    def calc_matrices(self):
        '''Calculate matrices for hydraulic problem formulation.'''
        
       
        
        # incidence_matrix
        B_mat = nx.linalg.graphmatrix.incidence_matrix(self.graph.G, nodelist=self.setup.v_vec,
                                                       edgelist = self.setup.e_vec, 
                                                       oriented=True)
        self.B_mat = B_mat
        
        #cycle basis matrix
        self.O_mat = self.cycle_basis_mat()
        
        
        
    def cycle_basis_mat(self):
 
        
        # cycle basis matrix
        loops = nx.cycle_basis(self.graph.G_prim,root=None)

        num_loops=len(loops)
        if num_loops != len(self.setup.e_vec)-len(self.setup.v_vec)+1:
            raise Exception('number of calculated basic circles in the network is '\
                'wrong')
        
        ## generate list with loops in form of adjacent nodes
        loops_adj=[0]*num_loops
        for i in range(len(loops_adj)):
            loops_adj[i]=[]
        del i
        for i in range(num_loops):
            for j in range(len(loops[i])):
                if j<(len(loops[i])-1):
                    loops_adj[i].append((loops[i][j], loops[i][j+1]))
                elif j==(len(loops[i])-1):
                    loops_adj[i].append((loops[i][j], loops[i][0]))
                else:
                    print('fail')
        del i,j

        ## convert edges to {-1,0,1} according to loops and create O_mat
        O_mat_list = [[0 for i in range(len(self.setup.e_vec))] for j in 
                           range(len(loops_adj))]
        for i in range(len(loops_adj)):
            for j in range(len(loops_adj[i])):
                for k in range(len(self.setup.e_vec)):
                    if loops_adj[i][j] == self.setup.e_vec[k]:
                        O_mat_list[i][k] = 1
                    elif (loops_adj[i][j][1], loops_adj[i][j][0]) == self.setup.e_vec[k]:
                        O_mat_list[i][k] = -1
        ### if number of -1s is larger then the number of 1s in a row, than mulitply by 
        ### -1
        for i in range(len(O_mat_list)):
            if O_mat_list[i].count(-1) > O_mat_list[i].count(1):
                O_mat_list[i] = [x*-1 for x in O_mat_list[i]]
        ## convert to sparse matrix
        O_mat = sp.sparse.csr_matrix(O_mat_list)
        #self.O_mat = O_mat 
        
        return O_mat