#!/usr/bin/env python
from IsotopicBoxModel import *

class ZnCalibration(IsotopicBoxModel):
    """ 
    A simple engine that perform the computation 
    """ 
    def __init__(self):
        """ Add options for the number of measures, migration bandwidth, number of nodes
        walltime, env_file or env_name, stress, and clusters and initialize the engine """
        super(ZnCalibration, self).__init__()
        self.parameters()
        self.delta_name = r"$\delta^{66}Zn$"
    
    def run(self):
        """ Execute the engine and compute the results """        
        Delta = self.initial_state()
        
        Delta = self.compute_evolution(Delta)
        
        self.final_state(Delta[-1,:])
 
    
    def parameters(self):
        """ Define the time parameters, the isotopic standard and  the boxes, flux and partition coefficients """
        n_timestep = 1000000
        self.time = linspace(0, 13870.0, n_timestep)  # temps
         
        # JMC standard
        self.standard = 0.565203
         
        self.Boxes = { 
            "diet":     {'Delta':    0e0, 'Mass':  1e12}, 
            "plasma":   {'Delta':    0e0, 'Mass':   3e0}, 
            "RBC":      {'Delta':    0e0, 'Mass': 2.5e1}, 
            "liver":    {'Delta':    0e0, 'Mass': 1.3e2}, 
            "urine":    {'Delta':    0e0, 'Mass':   1e1}, 
            "feces":    {'Delta':    0e0, 'Mass':   1e1}, 
            "muscle":   {'Delta':    0e0, 'Mass': 1.5e3},
            "bone":     {'Delta': 0.48e0, 'Mass': 7.7e2},
            "skin":     {'Delta':    0e0, 'Mass': 1.6e2},
            "kidney":   {'Delta':    0e0, 'Mass':   2e1}
        }
        self.Flux = {
            "diet":   {"diet": 0.0, "plasma": 4e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 8e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "plasma": {"diet": 0.0, "plasma": 0e0, "RBC": 0.18e0, "liver": 2.64e0, "urine": 0e0, "feces": 3e0, "muscle": 0.0035e0, "bone": 0.01, "skin": 0.5e0, "kidney":2.5e0},
            "RBC":    {"diet": 0.0, "plasma": 0.18e0, "RBC": 0e0, "liver": 0.0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "liver":  {"diet": 0.0, "plasma": 2.64e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "urine":  {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "feces":  {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "muscle": {"diet": 0.0, "plasma": 0.0035e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "bone":   {"diet": 0.0, "plasma": 0.01, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "skin":   {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "kidney": {"diet": 0.0, "plasma": 2.0e0, "RBC": 0e0, "liver": 0e0, "urine": 0.5e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0}
            }
        # Coefficients de partage
     #Coeff Monte Carlo
       # coeff_KU = 1 / 0.9998e0
        #coeff_PRBC = 1.00025e0
        #coeff_PS = 1.000275
        #coeff_PM = 0.99986
        #coeff_PB = 1.0003
        #coeff_PL = 0.99939
        #coeff_PD = 1.000
        #Coeff souris
        #coeff_KU = 1e0
        #coeff_PRBC = 1.00061e0
        #coeff_PS = 1.00027
        #coeff_PM = 0.99987
        #coeff_PB = 1.00053
        #coeff_PL = 0.99965
        #coeff_PD = 1.000
        
        #Coeff humain
        coeff_KU = 1e0
        coeff_PRBC = 1.00025e0
        coeff_PS = 0.99965
        coeff_PM = 0.99986
        coeff_PB = 1.0003
        coeff_PL = 0.99939
        coeff_PD = 1.000
        

        self.Partcoeff = {
            "diet":   {"diet": 1.0, "plasma": 1.00018e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone": 1e0, "skin": 1e0, "kidney":1e0},
            "plasma": {"diet": 1.0, "plasma": 1e0, "RBC": coeff_PRBC, "liver": coeff_PL, "urine": 1e0, "feces": coeff_PD, "muscle": coeff_PM, "bone":coeff_PB, "skin": coeff_PS, "kidney":1/coeff_KU},
            "RBC":    {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1/coeff_PRBC, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "liver":  {"diet": 1.0, "plasma": 1/coeff_PL, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "urine":  {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "feces":  {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "muscle": {"diet": 1.0, "plasma": 1/coeff_PM, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "bone":   {"diet": 1.0, "plasma": 1/coeff_PB, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "skin":   {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "kidney": {"diet": 1.0, "plasma": coeff_KU, "RBC": 1e0, "liver": 1e0, "urine": coeff_KU, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0}
            }
            
