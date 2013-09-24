#!/usr/bin/env python
'''
A box model written in Python

This tool aims to calculate the evolution of a system composed of boxes
that exchange flux.


It requires SciPy for computation, Matplotlib for model output
and execo_engine for parameter range exploration.


More documentation can be found in the README file.

Laurent Pouilloux and Klervia Jaouen, 
Ecole Normale Superieure de Lyon

'''
from random import gauss
from scipy.integrate import odeint
from execo_engine import Engine, ParamSweeper, sweep, slugify, logger
import numpy as np
import matplotlib.pyplot as plt 


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
        
        Ratio = odeint(self.evol_ratio, self.initial_state(), self.time)
        Delta_final = ((Ratio/self.standard_IRMM)-1.0)*1000;
        print Delta_final[:,]
#        plt.xkcd()
#        p2 = plt.plot(self.time, Delta_final[:,1], label="liver")
        
        
        i_box = 0
        for box, value in self.Delta.iteritems():
            print box, value
            if box not in [ 'feces', 'urine' ]:
                plt.plot(self.time/365., Delta_final[:,i_box], label=box)
            i_box += 1
        
        
        
        
        plt.legend()
        
                
        plt.xlabel(r"Years")
        plt.ylabel(r"$\delta^{66}Zn$$(permil)$")

        plt.savefig('test.png')        
        

#            Ratio.append((Delta[ii]/1e3+1e0)*Ratio_standard_IRMM)
        
        
    def initial_state(self):
        self.Mass = np.array( self.Mass.values() )
        self.Flux = np.array( [ box.values() for box in self.Flux.values() ])
        self.Partcoeff = np.array( [ box.values() for box in self.Partcoeff.values() ])
        Ratio = np.array( [ (delta/1e3+1e0)*self.standard_IRMM for delta in self.Delta.values() ] )
        return Ratio
        
        

#        plot_state(Boxes, Ratio, data_flux, i_mass, i_flux)
        
                
    def default_parameters(self):
        """ Define the boxes, the flux and the partition coefficients """
        n_timestep = 10000
        self.time = np.linspace(0, 13870.0, n_timestep)  # temps
         
        self.standard_IRMM = 0.0637 
        self.Delta = {
             "diet":     1.0e0, 
             "plasma":  1.51e0, 
             "RBC":     2.74e0, 
             "liver":   1.35e0, 
             "urine":    1.0e0, 
             "feces":    0.1e0, 
             "menses":   2.5e0
             }
        self.Mass = { 
            "diet":     1e12, 
            "plasma":   3e0, 
            "RBC":      2.5e3, 
            "liver":    10e2, 
            "urine":    1e-10, 
            "feces":    1e-0, 
            "menses":   1e-2
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
       
#    def init_ratio(self):
#        """ Calculate """
        
    def evol_ratio(self, ratio, t):
        rationew = np.zeros(ratio.size)
        for ii in range(ratio.size):
            outflux=0;
            influx=0;
            for jj in range(ratio.size):
                # le nouveau ratio est calcule a partir des flux entrants et sortants
                outflux = outflux + self.Flux[ii][jj]/self.Mass[ii]*self.Partcoeff[ii][jj]*ratio[ii]
                influx = influx + self.Flux[jj][ii]/self.Mass[ii]*self.Partcoeff[jj][ii]*ratio[jj]
            rationew[ii]= influx - outflux
        return rationew;