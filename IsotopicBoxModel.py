#!/usr/bin/env python
'''
A box model written in Python to simulate isotopic evolution in a human body

See README for details

Laurent Pouilloux and Klervia Jaouen, 
Ecole Normale Superieure de Lyon

This tools released under the GNU Public
License, version 3 or later.
'''
from pprint import pformat
from os import path, mkdir
from random import gauss
from execo_engine import Engine, ParamSweeper, sweep, slugify, logger
from execo.log import set_style
from numpy import linspace, array, zeros, absolute
from scipy.integrate import odeint
import matplotlib.pyplot as plt 
from pydot import Dot, Node, Edge


class IsotopicBoxModel(Engine):
    """ 
    This is the main engine that is used to created custom isotopic models
    """ 
    def __init__(self):
        """ Add options for the number of measures, migration bandwidth, number of nodes
        walltime, env_file or env_name, stress, and clusters and initialize the engine """
        super(IsotopicBoxModel, self).__init__()
        self.init_plots()
        logger.setLevel('INFO')
        logger.info(set_style('\n\n                 Welcome to the human isotopic Box Model\n', 'log_header'))
        logger.debug(pformat(self.__dict__))
            
    def initial_state(self, outdir = None):
        """ Convert the dict given from parameters to Numpy array """
        logger.info(set_style('Initial boxes configuration\n', 'log_header')+
                    ''.ljust(8)+''.join( [ set_style(box.rjust(10), 'emph') for box in self.Boxes.iterkeys() ])+
                    set_style('\n'+'Delta'.ljust(8), 'object_repr')+''.join( [ str(box['Delta']).rjust(10) for box in self.Boxes.itervalues()])+
                    set_style('\n'+'Mass'.ljust(8), 'object_repr')+''.join( [ str(box['Mass']).rjust(10) for box in self.Boxes.itervalues() ])
                    ) 
        if outdir is None:
            outdir = self.result_dir+'/'
        self.plot_state(self.Boxes.keys(), array( [ box['Delta'] for box in self.Boxes.itervalues()]), 
                        name = '_initial', outdir = outdir)
        
        self._Mass = array( [ box['Mass'] for box in self.Boxes.itervalues() ] )
        self._Flux = array( [ box.values() for box in self.Flux.values() ])
        self._Partcoeff = array( [ box.values() for box in self.Partcoeff.values() ])
        
        f = open(outdir+'/Delta.initial', 'w')
        for box, value in self.Boxes.iteritems():
            f.write(box+' '+str(value['Delta'])+'\n')
        f.close()
        return [ box['Delta'] for box in self.Boxes.itervalues() ]
        
    
    def compute_evolution(self, Delta, func = None, outdir = None):
        """ """
        logger.info(set_style('Computing evolution', 'log_header'))
        if func is None:
            func = self.evol_ratio
        if outdir is None:
            outdir = self.result_dir
        Ratio = [ ( delta/1e3+1e0)*self.standard for delta in Delta ] 
        Ratio = odeint(func, Ratio, self.time)
        Delta = ((Ratio/self.standard)-1.0)*1000;        
        self.plot_evolution(Delta, outdir = outdir)
        return Delta
    
    
        
    def evol_ratio(self, ratio, t):
        """ The evolution function that can be used for isotopic ratio evolution"""
        rationew = zeros(ratio.size)
        for ii in range(ratio.size):
            outflux = 0;
            influx = 0;
            for jj in range(ratio.size):
                outflux = outflux + self._Flux[ii][jj]/self._Mass[ii]*self._Partcoeff[ii][jj]*ratio[ii]
                influx = influx + self._Flux[jj][ii]/self._Mass[ii]*self._Partcoeff[jj][ii]*ratio[jj]
            rationew[ii]= influx - outflux
        return rationew;
    
    def final_state(self, Delta_final, outdir = None):
        """ """
        if outdir is None:
            outdir = self.result_dir
        logger.info(set_style('Final boxes state\n', 'log_header')+
                    ''.ljust(8)+''.join( [ set_style(box.rjust(10), 'emph') for box in self.Boxes.iterkeys() ])+
                    set_style('\n'+'Delta'.ljust(8), 'object_repr')+''.join( [ str(round(delta, 7)).rjust(10) for delta in Delta_final if absolute(delta) < 1000]))
        self.plot_state(self.Boxes.keys(), Delta_final, name = '_final', outdir = outdir)
        f = open(outdir+'/Delta.final', 'w')
        for box in self.Boxes.iterkeys():
            idx = self.Boxes.keys().index(box)
            f.write(box+' '+str(Delta_final[idx])+'\n')
        f.close()
        
        
    def init_plots(self):
        """ Define the colors and shape of the model boxes"""
        self.plots_conf = {
            "diet":   {'color': "#BBFFB5", 'shape': "rectangle"},
            "plasma": {'color': "#D6DE42", 'shape': "ellipse"},
            "RBC":    {'color': "#FF3A25", 'shape': "ellipse"},
            "liver":  {'color': "#93BB8F", 'shape': "ellipse"}, 
            "urine":  {'color': "#FBFF93", 'shape': "rectangle"},
            "feces":  {'color': "#BF9285", 'shape': "rectangle"},
            "menses": {'color': "#FF0600", 'shape': "rectangle"},
            "kidney": {'color': "#620A00", 'shape': "ellipse"},
            "muscle": {'color': "#FF2933", 'shape': "ellipse"},
            "skin":   {'color': "#FFE486", 'shape': "ellipse"} ,
            "bone":   {'color': "#C9C985", 'shape': "ellipse"}  
            }
        self.color_chars = '0123456789ABCDEF'
        
    
    def plot_state(self, boxes, deltas, name = '', outdir = None):
        """ Make a graph of a given state """
        graph = Dot(graph_type='digraph', fontname="Verdana", size="10, 5", fixedsize= True)
        i_box = 0
        for box in boxes:            
            textcolor = 'white' if sum( [ self.color_chars.index(col) for col in self.plots_conf[box]['color'].split('#')[1] ] ) < 35 else 'black' 
            node_box = Node(box, style="filled", label = '<<font POINT-SIZE="10" color="'+textcolor+'">'+box+'<br/> '+
                            "%.7f" % round(deltas[i_box], 7)+'</font>>',
                fillcolor = self.plots_conf[box]['color'], shape = self.plots_conf[box]['shape'])
            i_box += 1
            graph.add_node(node_box)
            
        for box_from, boxes_to in self.Flux.iteritems():
            for box_to, flux in boxes_to.iteritems():
                if flux !=0:
                    if flux > 0:
                        edge = Edge(box_from, box_to,  label = '<<font POINT-SIZE="10">'+str(flux)+'</font>>')
                    elif flux < 0:
                        edge = Edge(box_to, box_from,  label = '<<font POINT-SIZE="10">'+str(flux)+'</font>>')                
                    graph.add_edge(edge)
        
        if outdir is None:
            outdir = self.result_dir
            
        outfile = outdir+'/state'+name+'.png'
        graph.write_png(outfile)
        logger.info('State has been saved to '+set_style(outfile, 'emph'))
        
    def plot_evolution(self, Delta, outdir = None):
        """ Draw a graph of the boxes evolution through years"""
        plt.figure()
        i_box = 0
        for box in self.Boxes:
            # remove deriving boxes 
            if absolute(Delta[0,i_box]-Delta[-1,i_box]) < 1000:
                plt.plot(self.time/365., Delta[:,i_box], label=box,
                         color = self.plots_conf[box]['color'])
            i_box += 1
        plt.legend()
        plt.xlabel(r"Years")
        plt.ylabel(self.delta_name)
        if outdir is None:
            outdir = self.result_dir
        outfile = outdir+'/evolution.png'
        plt.savefig(outfile)
        logger.info('Evolution has been saved to '+set_style(outfile, 'emph'))    
        
        
        
        