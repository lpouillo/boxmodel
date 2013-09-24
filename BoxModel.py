#!/usr/bin/env python
'''
A box model written in Python

See README for details

Laurent Pouilloux and Klervia Jaouen, 
Ecole Normale Superieure de Lyon

This tools released under the GNU Public
License, version 3 or later.
'''
from pprint import pformat
from random import gauss
from execo_engine import Engine, ParamSweeper, sweep, slugify, logger
from numpy import linspace, array, zeros
from scipy.integrate import odeint
import matplotlib.pyplot as plt 
import pydot as P


class BoxModel(Engine):
    """ 
    This is the main engine that is used to created custom isotopic models
    """ 
    def __init__(self):
        """ Add options for the number of measures, migration bandwidth, number of nodes
        walltime, env_file or env_name, stress, and clusters and initialize the engine """
        super(BoxModel, self).__init__()
    
    def initial_state(self):
        """ Convert the dict given from parameters to Numpy array """
        logger.debug('MASS\n'+pformat([box['Mass'] for box in self.Boxes.itervalues() ])) 
        self._Mass = array( [box['Mass'] for box in self.Boxes.itervalues() ] )
        self._Flux = array( [ box.values() for box in self.Flux.values() ])
        self._Partcoeff = array( [ box.values() for box in self.Partcoeff.values() ])
        Ratio = array( [ (box['Delta']/1e3+1e0)*self.standard for box in self.Boxes.itervalues() ] )    
        return Ratio
    
    def evol_ratio(self, ratio, t):
        """ The evolution function that can be used in  """
        rationew = zeros(ratio.size)
        for ii in range(ratio.size):
            outflux = 0;
            influx = 0;
            for jj in range(ratio.size):
                outflux = outflux + self._Flux[ii][jj]/self._Mass[ii]*self._Partcoeff[ii][jj]*ratio[ii]
                influx = influx + self._Flux[jj][ii]/self._Mass[ii]*self._Partcoeff[jj][ii]*ratio[jj]
            rationew[ii]= influx - outflux
        return rationew;
    
    def plot_state(self, boxes, deltas, name = ''):
        colors = {"diet": "#BBFFB5", "plasma": "#FFC66D", "RBC": "#BA0400", 
                "liver": "#93BB8F", "urine": "#FBFF93", "feces": "#BF9285", "menses": "#FF0600"}
        shapes = {"diet": "rectangle", "plasma": "ellipse", "RBC": "ellipse", 
                "liver": "ellipse", "urine": "rectangle", "feces": "rectangle", "menses": "rectangle"}
        graph = P.Dot(graph_type='digraph', fontname="Verdana", ratio = "1")
        i_box = 0
        for box in boxes:
            node_box = P.Node(box, style="filled", label = box+'\n '+str(round(deltas[i_box], 7)),
                fillcolor=colors[box], shape = shapes[box])
            i_box += 1
            graph.add_node(node_box)
            
        for box_from, boxes_to in self.Flux.iteritems():
            for box_to, flux in boxes_to.iteritems():
                if flux > 0:
                    edge = P.Edge(box_from, box_to,  label = flux)
                    graph.add_edge(edge)
        graph.write_png(self.result_dir+'/state'+name+'.png')

        
    def plot_evolution(self, Delta):
        i_box = 0
        for box in self.Boxes:
            if box not in [ 'feces', 'urine' ]:
                plt.plot(self.time/365., Delta[:,i_box], label=box)
            i_box += 1
        plt.legend()
        plt.xlabel(r"Years")
        plt.ylabel(r"$\delta^{66}Zn$$(permil)$")
        plt.savefig(self.result_dir+'/evolution.png')         
        
        
        
        