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
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

def visualize_graph(setup, graph, solutions):
        
        # initializations
        G = graph.G
        
        # parameters
        myfontsize = 20
        myfigsize=[2*21.0/2.54, 2*12.98/2.54]
        myarrowstyle = mpltlib.patches.ArrowStyle.CurveFilledAB(head_length=0.6, head_width=0.4) #CurveB, CurveFilledB
        savepath = './figs/'
        saveformat = 'png' # jpg, png, eps, svgd
        
        # initialize figure
        fig1 = plt.figure(num='graph', figsize=myfigsize)
        plt.title('graph', {'fontsize': 30,
            'fontweight' : 'bold'})
        
        # set-up node color
        for node in G.nodes():
            if 'h' in node:
                G.nodes[node]['color']='r'
            if 'c' in node:
                G.nodes[node]['color']='b'
        
        nodes_colors = {node: G.nodes[node]['color'] for node in 
                        G.nodes()}

        colorslist = list(nodes_colors.values())
        
        # draw
            
        nx.drawing.nx_pylab.draw(G, 
                                 pos=setup.coordinates,
                                 fontsize = 20,
                                 with_labels=True, node_color = colorslist,
                                 arrowsize=15, arrowstyle=myarrowstyle,
                                 node_size=700, width=3,
                                 font_weight = 'bold', alpha=1)
        
        xcoords= [setup.coordinates[node][0]
                  for node in setup.coordinates]
        ycoords= [setup.coordinates[node][1]
                  for node in setup.coordinates]
        
        plt.gca().set_xlim(min(xcoords)-30, max(xcoords)+30)
        plt.gca().set_ylim(min(ycoords)-30, max(ycoords)+30)
        
        plt.axis('off')
        plt.show(block=False)
        
        plt.savefig(''.join((savepath, 'vis_graph.', saveformat)), format = saveformat)
    
def visualize_scenario(setup, graph, solutions):
        # initializations
        G = graph.G
        img_valve = mpimg.imread('figs/valve_img.png')
        img_hx = mpimg.imread('figs/hx_img.png')
        img_pump = mpimg.imread('figs/pump_img.png')
        
        mu_dict = {P: setup.scenario[P]['mu'] for P in setup.PSM}
        u_dict = {P: setup.scenario[P]['u'] for P in setup.PSM}
        kappa_dict = {P: setup.scenario[P]['kappa'] for P in setup.PSM}
        T_sec_dict = {P: setup.scenario[P]['T_sec_in_degree'] for P in setup.PSM}
        dotV_sec_dict = {P: setup.scenario[P]['dotV_sec_in'] for P in setup.PSM}
        
        # parameters
        myfontsize = 20
        myfigsize=[2*21.0/2.54, 2*12.98/2.54]
        myarrowstyle = mpltlib.patches.ArrowStyle.CurveFilledAB(head_length=0.6, head_width=0.4) #CurveB, CurveFilledB
        savepath = './figs/'
        saveformat = 'png' # jpg, png, eps, svgd
        
        # initialize figure
        fig = plt.figure(num='graph', figsize=myfigsize)
        plt.title('graph', {'fontsize': 30,
            'fontweight' : 'bold'})
        
        # set-up node color
        for node in G.nodes():
            if 'h' in node:
                G.nodes[node]['color']='r'
            if 'c' in node:
                G.nodes[node]['color']='b'
        
        nodes_colors = {node: G.nodes[node]['color'] for node in 
                        G.nodes()}

        colorslist = list(nodes_colors.values())
        
        # draw
            
        nx.drawing.nx_pylab.draw(G, 
                                 pos=setup.coordinates,
                                 fontsize = 20,
                                 with_labels=True, node_color = colorslist,
                                 arrowsize=15, arrowstyle=myarrowstyle,
                                 node_size=700, width=3,
                                 font_weight = 'bold', alpha=1)
        
        xcoords= [setup.coordinates[node][0]
                  for node in setup.coordinates]
        ycoords= [setup.coordinates[node][1]
                  for node in setup.coordinates]
        
        plt.gca().set_xlim(min(xcoords)-30, max(xcoords)+30)
        plt.gca().set_ylim(min(ycoords)-30, max(ycoords)+30)
        
        plt.axis('off')
        plt.show(block=False)
        

                
        for PSM in mu_dict:
            for i in range(len(setup.e_vec)):
                edge = setup.e_vec[i]
                if (str(PSM) in edge[0]) and (str(PSM) in edge[1]):
                    coord_hx = []
                    coord_hx.append(setup.coordinates[edge[0]][0]
                                     +(setup.coordinates[edge[1]][0]-
                                     setup.coordinates[edge[0]][0])*(1/2))
                    coord_hx.append(setup.coordinates[edge[0]][1]
                                     +(setup.coordinates[edge[1]][1]-
                                       setup.coordinates[edge[0]][1])*(2/5))
                    coord_act = [coord_hx[0], 
                                 setup.coordinates[edge[0]][1]
                                     +(setup.coordinates[edge[1]][1]-
                                       setup.coordinates[edge[0]][1])*(3/5)]
                    
                    coord_act_disp = plt.gca().transData.transform(coord_act)
                    coord_hx_disp = plt.gca().transData.transform(coord_hx)
                    
                    imagebox_hx = OffsetImage(img_hx, zoom=0.1)
                    ab_hx = AnnotationBbox(imagebox_hx, (coord_hx[0],coord_hx[1]), 
                                           frameon=False)
                    plt.gca().add_artist(ab_hx)
                    
                    fig.canvas.draw()
                    renderer = fig.canvas.renderer
                    inv1=plt.gca().transData.inverted()
                    inv2=plt.gca().transAxes.inverted()
                    
                    U_pump_str=r'$u$ = %3.1f' % (u_dict[PSM])
                    U_valve_str=r'$\kappa = %3.1f$' % (kappa_dict[PSM])
                    dotm_sec_str = r'$\dot{V}_{sec} = %4.2f \frac{l}{min}$' % (dotV_sec_dict[PSM])
                    T_sec_str = r'$T_{sec,in} = %5.1f ^\circ C$' % (T_sec_dict[PSM])
                    
                    myExtend0=imagebox_hx.get_extent(renderer)
                    xpos_text_disp0a = coord_hx_disp[0]+0.55*myExtend0[0]
                    xpos_text_disp0b = coord_hx_disp[0]+0.55*myExtend0[0]
                    ypos_text_disp0a = coord_hx_disp[1]+0.4*myExtend0[0]
                    ypos_text_disp0b = coord_hx_disp[1]-0.4*myExtend0[0]
                    coord_text_data0a=inv1.transform((
                            xpos_text_disp0a,ypos_text_disp0a))
                    coord_text_data0b=inv1.transform((
                            xpos_text_disp0b,ypos_text_disp0b))
                    plt.text(coord_text_data0a[0], coord_text_data0a[1], dotm_sec_str, 
                             horizontalalignment='left',
                             verticalalignment='center',
                             fontsize = myfontsize)
                    plt.text(coord_text_data0b[0], coord_text_data0b[1], T_sec_str, 
                             horizontalalignment='left',
                             verticalalignment='center',
                             fontsize = myfontsize)
                    
                    if mu_dict[PSM]==1:
                        
                        imagebox_pump = OffsetImage(img_pump, zoom=0.19)
                        myExtend1=imagebox_pump.get_extent(renderer)
                        xpos_text_disp1 = coord_act_disp[0]+0.55*myExtend1[0]
                        coord_text_data1=inv1.transform((
                            xpos_text_disp1,coord_act_disp[1]))
                        ab_pump = AnnotationBbox(imagebox_pump, 
                                                 (coord_act[0],coord_act[1]),
                                                 box_alignment=(0.5, 0.5),
                                                 pad=0,
                                                 frameon=False)
                        plt.gca().add_artist(ab_pump)
                        plt.text(coord_text_data1[0], coord_text_data1[1], U_pump_str, 
                             horizontalalignment='left',
                             verticalalignment='center',
                             fontsize = myfontsize)
                    elif mu_dict[PSM]== -1:
                        # position valve image
                        imagebox_valve = OffsetImage(img_valve, zoom=0.4)
                        myExtend2=imagebox_valve.get_extent(renderer)
                        
                        xpos_valve_disp = coord_act_disp[0]-0.31*myExtend2[0]
                        
                        coord_valve_data=inv1.transform((xpos_valve_disp,
                                                         coord_act_disp[1]))
                        #coord_valve_axes=inv2.transform((xpos_valve_disp,coord_act_disp[1]))
                        
                        ab_valve = AnnotationBbox(imagebox_valve, 
                                                  coord_valve_data,
                                                  box_alignment=(0.5, 0.5),
                                                  pad=0,
                                                  frameon=False)
                                                    #xycoords='axes fraction'
                        plt.gca().add_artist(ab_valve)
                        
                        # position text for control value U
                        xpos_text_disp2 = xpos_valve_disp+0.55*myExtend2[0]
                        coord_text_data2=inv1.transform((
                            xpos_text_disp2,coord_act_disp[1]))
                        plt.text(coord_text_data2[0], coord_text_data2[1], U_valve_str, 
                             horizontalalignment='left',
                             verticalalignment='center',
                             fontsize = myfontsize)
        
        plt.axis('off')
        plt.show(block=False)
        
        plt.savefig(''.join((savepath, 'vis_scenario.', saveformat)), format = saveformat)
    
def visualize_hy_solution(setup, graph, solutions):
    def visualize_sol_graph(self, withsymbols=False, withnumbers=True):
        '''code by M.Sc. Thomas Licklederer, Technical University of Munich, MSE, all
        rights reserved DESCRIPTION TO BE ADDED'''
        
        # combine the solution values with the correct edges
        dmdt_soldict = dict(zip(self.networkgraph.e_vec, self.dmdt_solvec))
        Deltap_soldict = dict(zip(self.networkgraph.e_vec, self.Deltap_solvec))
        
        # assign the solution values to the Graph
        for edge in self.networkgraph.e_vec:
            self.networkgraph.G.edges[edge]['dmdt_sol']=dmdt_soldict[edge]
            self.networkgraph.G.edges[edge]['Deltap_sol']=Deltap_soldict[edge]
        
        # initialize new Graph that represents the solution
        self.networkgraph.Gsol_hy = nx.DiGraph()
        self.networkgraph.Gsol_hy.add_nodes_from(self.networkgraph.v_vec)
        myedgelabels = {}
        for i in range(len(self.networkgraph.e_vec)):
            edge = self.networkgraph.e_vec[i]
            if np.sign(dmdt_soldict[edge])==-1:
                self.networkgraph.Gsol_hy.add_edge(edge[1],edge[0])
                self.networkgraph.Gsol_hy.edges[(edge[1],edge[0])]['dmdt']=-dmdt_soldict[(edge[0],edge[1])]
                self.networkgraph.Gsol_hy.edges[(edge[1],edge[0])]['Deltap']=-Deltap_soldict[(edge[0],edge[1])]
                myedgelabels[(edge[1],edge[0])]=[]
                myedgelabels[(edge[1],edge[0])].append('{0:6.2f} kg/min'.format(-dmdt_soldict[(edge[0],edge[1])]*60))
                myedgelabels[(edge[1],edge[0])].append('{0:6.1f} mbar'.format(-Deltap_soldict[(edge[0],edge[1])]/100))
            elif np.sign(dmdt_soldict[edge])==1:
                self.networkgraph.Gsol_hy.add_edge(edge[0],edge[1])
                self.networkgraph.Gsol_hy.edges[(edge[0],edge[1])]['dmdt']=dmdt_soldict[(edge[0],edge[1])]
                self.networkgraph.Gsol_hy.edges[(edge[0],edge[1])]['Deltap']=Deltap_soldict[(edge[0],edge[1])]
                myedgelabels[(edge[0],edge[1])]=[]
                myedgelabels[(edge[0],edge[1])].append('{0:6.2f} kg/min'.format(-dmdt_soldict[(edge[0],edge[1])]*60))
                myedgelabels[(edge[0],edge[1])].append('{0:6.1f} mbar'.format(-Deltap_soldict[(edge[0],edge[1])]/100))
            elif np.sign(dmdt_soldict[edge])==0:
                self.networkgraph.Gsol_hy.add_edge(edge[0],edge[1])
                self.networkgraph.Gsol_hy.edges[(edge[0],edge[1])]['dmdt']=dmdt_soldict[(edge[0],edge[1])]
                self.networkgraph.Gsol_hy.edges[(edge[0],edge[1])]['Deltap']=Deltap_soldict[(edge[0],edge[1])]
                myedgelabels[(edge[0],edge[1])]=[]
                myedgelabels[(edge[0],edge[1])].append('{0:6.2f} kg/min'.format(-dmdt_soldict[(edge[0],edge[1])]*60))
                myedgelabels[(edge[0],edge[1])].append('{0:6.1f} mbar'.format(-Deltap_soldict[(edge[0],edge[1])]/100))
            else:
                raise ValueError('')
        
        # plot graph
        
        myfontsize = 16
        myfigsize=[2*21.0/2.54, 2*12.98/2.54]
        myarrowstyle = mpltlib.patches.ArrowStyle.CurveFilledB(head_length=0.6, head_width=0.4) #CurveB
        
        fig1 = plt.figure(num='hydraulic solution', figsize=myfigsize)
        plt.title('hydraulic solution', {'fontsize': 30,
            'fontweight' : 'bold'})
        
        for node in self.networkgraph.Gsol_hy.nodes():
            if 'w' in node:
                self.networkgraph.Gsol_hy.nodes[node]['color']='r'
            if 'c' in node:
                self.networkgraph.Gsol_hy.nodes[node]['color']='b'
        
        nodes_colors = {node: self.networkgraph.Gsol_hy.nodes[node]['color'] for node in 
                        self.networkgraph.Gsol_hy.nodes()}

        colorslist = list(nodes_colors.values())
            
        nx.drawing.nx_pylab.draw(self.networkgraph.Gsol_hy, 
                                 pos=self.networkgraph.coordinates,
                                 fontsize = 'xx-large',
                                 with_labels=True, node_color = colorslist,
                                 arrowsize=15, arrowstyle=myarrowstyle,
                                 node_size=700, width=3,
                                 font_weight = 'bold', alpha=1)
        
        xcoords= [self.networkgraph.coordinates[node][0]
                  for node in self.networkgraph.coordinates]
        ycoords= [self.networkgraph.coordinates[node][1]
                  for node in self.networkgraph.coordinates]
        
        plt.gca().set_xlim(min(xcoords)-30, max(xcoords)+30)
        plt.gca().set_ylim(min(ycoords)-30, max(ycoords)+30)
                
        img_valve = mpimg.imread('pictures/valve_img.png')
        img_hx = mpimg.imread('pictures/hx_img.png')
        img_pump = mpimg.imread('pictures/pump_img.png')
        
        if withsymbols == True:
            for PSM in self.mu_dict:
                for i in range(len(self.networkgraph.e_vec)):
                    edge = self.networkgraph.e_vec[i]
                    if (str(PSM) in edge[0]) and (str(PSM) in edge[1]):
                        coord_hx = []
                        coord_hx.append(self.networkgraph.coordinates[edge[0]][0]
                                        +(self.networkgraph.coordinates[edge[1]][0]-
                                        self.networkgraph.coordinates[edge[0]][0])*(1/2))
                        coord_hx.append(self.networkgraph.coordinates[edge[0]][1]
                                        +(self.networkgraph.coordinates[edge[1]][1]-
                                        self.networkgraph.coordinates[edge[0]][1])*(2/5))
                        coord_act = [coord_hx[0], 
                                    self.networkgraph.coordinates[edge[0]][1]
                                        +(self.networkgraph.coordinates[edge[1]][1]-
                                        self.networkgraph.coordinates[edge[0]][1])*(3/5)]
                        coord_act_disp = plt.gca().transData.transform(coord_act)
                        
                        imagebox_hx = OffsetImage(img_hx, zoom=0.1)
                        ab_hx = AnnotationBbox(imagebox_hx, (coord_hx[0],coord_hx[1]), 
                                            frameon=False)
                        plt.gca().add_artist(ab_hx)
                        
                        fig1.canvas.draw()
                        renderer = fig1.canvas.renderer
                        inv1=plt.gca().transData.inverted()
                        inv2=plt.gca().transAxes.inverted()
                        
                        U_pump_str="U = %s" % (self.U_pumps_dict[PSM])
                        U_valve_str="U = %s" % (self.U_valves_dict[PSM])
                        
                        if self.mu_dict[PSM]==1:
                            imagebox_pump = OffsetImage(img_pump, zoom=0.19)
                            myExtend1=imagebox_pump.get_extent(renderer)
                            xpos_text_disp1 = coord_act_disp[0]+0.55*myExtend1[0]
                            coord_text_data1=inv1.transform((
                                xpos_text_disp1,coord_act_disp[1]))
                            ab_pump = AnnotationBbox(imagebox_pump, 
                                                    (coord_act[0],coord_act[1]),
                                                    box_alignment=(0.5, 0.5),
                                                    pad=0,
                                                    frameon=False)
                            plt.gca().add_artist(ab_pump)
                            plt.text(coord_text_data1[0], coord_text_data1[1], U_pump_str, 
                                horizontalalignment='left',
                                verticalalignment='center',
                                fontsize = myfontsize)
                        elif self.mu_dict[PSM]==-1:
                            # position valve image
                            imagebox_valve = OffsetImage(img_valve, zoom=0.4)
                            myExtend2=imagebox_valve.get_extent(renderer)
                            
                            xpos_valve_disp = coord_act_disp[0]-0.31*myExtend2[0]
                            
                            coord_valve_data=inv1.transform((xpos_valve_disp,
                                                            coord_act_disp[1]))
                            coord_valve_axes=inv2.transform((xpos_valve_disp,coord_act_disp[1]))
                            
                            ab_valve = AnnotationBbox(imagebox_valve, 
                                                    coord_valve_data,
                                                    box_alignment=(0.5, 0.5),
                                                    pad=0,
                                                    frameon=False)
                                                        #xycoords='axes fraction'
                            plt.gca().add_artist(ab_valve)
                            
                            # position text for control value U
                            xpos_text_disp2 = xpos_valve_disp+0.55*myExtend2[0]
                            coord_text_data2=inv1.transform((
                                xpos_text_disp2,coord_act_disp[1]))
                            plt.text(coord_text_data2[0], coord_text_data2[1], U_valve_str, 
                                horizontalalignment='left',
                                verticalalignment='center',
                                fontsize = myfontsize)
        if withnumbers == True:
            nx.drawing.nx_pylab.draw_networkx_edge_labels(self.networkgraph.Gsol_hy, 
                                 pos=self.networkgraph.coordinates,
                                 fontsize = 20,
                                 edge_labels = myedgelabels, node_color = colorslist,
                                 arrowsize=15, arrowstyle=myarrowstyle,
                                 node_size=700, width=3,
                                 font_weight = 'normal', alpha=1)
        
        plt.axis('off')
        plt.show(block=False)
    
def visualize_th_solution(setup, graph, solutions):
    pass
    