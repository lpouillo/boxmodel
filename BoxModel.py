#!/usr/bin/env python
'''
A box model written in Python

This tool aims to calculate the evolution of a system composed of boxes
that exchange flux.


It requires SciPy for computation, Matplotlib for model output
and execo_engine for parameter range exploration and engine characteristics.

use with 
execo-run BoxModel -ML


More documentation can be found in the README file.

Laurent Pouilloux and Klervia Jaouen, 
Ecole Normale Superieure de Lyon

'''
from random import gauss
from scipy.integrate import odeint
from execo_engine import Engine, ParamSweeper, sweep, slugify, logger
import numpy as np
import matplotlib.pyplot as plt 
import pydot as P


class BoxModel(Engine):
    """ 
    This is the main engine that allow to define the simulation parameter, the model boxes and
    range of scenarios for your experiments
    """ 
    def __init__(self):
        """ Add options for the number of measures, migration bandwidth, number of nodes
        walltime, env_file or env_name, stress, and clusters and initialize the engine """
        super(BoxModel, self).__init__()
        self.default_parameters()

    def run(self):
        """ Execute the engine and compute the results """
        Ratio = self.initial_state()
        
        Ratio = odeint(self.evol_ratio, Ratio, self.time)
        Delta_final = ((Ratio/self.standard_IRMM)-1.0)*1000;
        self.plot_evolution(Delta_final)
        self.plot_state(self.Boxes.keys(), Delta_final[:,0], name = '_initial')
        self.plot_state(self.Boxes.keys(), Delta_final[:,-1], name = '_final')
            
    def evol_ratio(self, ratio, t):
        """ The evolution function coming from """
        rationew = np.zeros(ratio.size)
        for ii in range(ratio.size):
            outflux=0;
            influx=0;
            for jj in range(ratio.size):
                outflux = outflux + self._Flux[ii][jj]/self._Mass[ii]*self._Partcoeff[ii][jj]*ratio[ii]
                influx = influx + self._Flux[jj][ii]/self._Mass[ii]*self._Partcoeff[jj][ii]*ratio[jj]
            rationew[ii]= influx - outflux
        return rationew;
        
    def default_parameters(self):
        """ Define the boxes, the flux and the partition coefficients """
        n_timestep = 10000
        self.time = np.linspace(0, 13870.0, n_timestep)  # temps
         
        self.standard_IRMM = 0.0637 
        self.Boxes = { 
            "diet":     {'Delta':  1.0e0, 'Mass':  1e12}, 
            "plasma":   {'Delta': 1.51e0, 'Mass':   3e0}, 
            "RBC":      {'Delta': 2.74e0, 'Mass': 2.5e3}, 
            "liver":    {'Delta': 1.35e0, 'Mass':  10e2}, 
            "urine":    {'Delta':  1.0e0, 'Mass': 1e-10}, 
            "feces":    {'Delta':  0.1e0, 'Mass':  1e-0}, 
            "menses":   {'Delta':  2.5e0, 'Mass':  1e-2} 
        }
        self.Flux = {
            "diet":     {"diet": 0.0, "plasma":  1.3, "RBC":  0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 },
            "plasma":   {"diet": 0.0, "plasma":  0.0, "RBC": 24.4, "liver": 5.0, "urine": 0.1, "feces": 0.0, "menses": 0.0 },
            "RBC":      {"diet": 0.0, "plasma": 23.9, "RBC":  0.0, "liver": 0.0, "urine": 0.0, "feces": 0.5, "menses": 0.0 },
            "liver":    {"diet": 0.0, "plasma":  4.3, "RBC":  0.0, "liver": 0.0, "urine": 0.0, "feces": 0.7, "menses": 0.0 },
            "urine":    {"diet": 0.0, "plasma":  0.0, "RBC":  0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 },
            "feces":    {"diet": 0.0, "plasma":  0.0, "RBC":  0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 },
            "menses":   {"diet": 0.0, "plasma":  0.0, "RBC":  0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 }
            }
        self.Partcoeff = {
            "diet":     {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 },
            "plasma":   {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 },
            "RBC":      {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 },
            "liver":    {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 },
            "urine":    {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 },
            "feces":    {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 },
            "menses":   {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses": 1e0 }
            }
    
    def initial_state(self):
        """ Convert the dict given from parameters to Numpy array """ 
        self._Mass = np.array( [box['Mass'] for box in self.Boxes.itervalues() ] )
        self._Flux = np.array( [ box.values() for box in self.Flux.values() ])
        self._Partcoeff = np.array( [ box.values() for box in self.Partcoeff.values() ])
        Ratio = np.array( [ (box['Delta']/1e3+1e0)*self.standard_IRMM for box in self.Boxes.itervalues() ] )    
        return Ratio
    
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