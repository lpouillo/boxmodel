#!/usr/bin/env python
from BoxModel import *

class FeModel(BoxModel):
    """ 
    This is the main engine that allow to define the simulation parameter, the model boxes and
    range of scenarios for your experiments
    """ 
    def __init__(self):
        """ Add options for the number of measures, migration bandwidth, number of nodes
        walltime, env_file or env_name, stress, and clusters and initialize the engine """
        super(FeModel, self).__init__()
        self.parameters()
    
    def run(self):
        """ Execute the engine and compute the results """
        Ratio = self.initial_state()
        Ratio = odeint(self.evol_ratio, Ratio, self.time)
        Delta_final = ((Ratio/self.standard_IRMM)-1.0)*1000;
        self.plot_evolution(Delta_final)
        self.plot_state(self.Boxes.keys(), Delta_final[:,0], name = '_initial')
        self.plot_state(self.Boxes.keys(), Delta_final[:,-1], name = '_final')
 
    
    def parameters(self):
        """ Define the boxes, the flux and the partition coefficients """
        n_timestep = 10000
        self.time = linspace(0, 13870.0, n_timestep)  # temps
         
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