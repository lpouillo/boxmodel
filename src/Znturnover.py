#!/usr/bin/env python
from IsotopicBoxModel import *

#from guppy import hpy


class Znturnover(IsotopicBoxModel):
    """
A simple engine that perform the computation
"""
    def __init__(self):
        """ Add options for the number of measures, migration bandwidth, number of nodes
walltime, env_file or env_name, stress, and clusters and initialize the engine """
        super(Znturnover, self).__init__()
        self.parameters()
        self.delta_name = r"$\delta^{66}Zn$"
    
    def run(self):
        """ Execute the engine and compute the results """
        self.parameters = {'flux_diet': range(7, 18), 'flux_bone': [0.0029 * i for i in range(1, 28, 1) ]}
        sweeps = sweep(self.parameters)
        sweeper = ParamSweeper( path.join(self.result_dir, "sweeps"), sweeps)
#        h = hpy()
#        print h.heap()
        while len(sweeper.get_remaining()) >0:
            comb = sweeper.get_next()
            comb_dir = self.result_dir +'/'+ slugify(comb)
            try:
                mkdir(comb_dir)
            except:
                pass
            self.set_flux(comb['flux_diet'], comb['flux_bone'])
            Delta = []
            
            Delta = self.initial_state(outdir = comb_dir)
            Delta = self.compute_evolution(Delta, outdir = comb_dir)
            self.final_state(Delta[-1,:], outdir = comb_dir)
            sweeper.done(comb)
#            print h.heap()
            gc.collect()
            logger.info('Combination done\n')
 
        
        logger.info('All combinations have been done, result can be found in '+self.result_dir)
        self.peace_flag(['RBC'], 'flux_bone', 'flux_diet')
        
        
        
    def set_flux(self, flux_diet, flux_bone):
        k = 0.33
        t = flux_diet/12
        logger.info('Using flux_diet = '+str(flux_diet)+' and flux_bone = '+str(flux_bone))
        self.Flux = {
            "diet": {"diet": 0.0, "plasma": k*flux_diet, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": (1-k)*flux_diet, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "plasma": {"diet": 0.0, "plasma": 0e0, "RBC": t*0.18e0, "liver": t*2.64e0, "urine": 0e0, "feces": 0.75*k*flux_diet, "muscle": t*0.0035e0, "bone": flux_bone, "skin": 0.125*k*flux_diet, "kidney":0.625*k*flux_diet},
            "RBC": {"diet": 0.0,"plasma": t*0.18e0, "RBC": 0e0, "liver": 0.0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "liver": {"diet": 0.0,"plasma": t*2.64e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "urine": {"diet": 0.0,"plasma": 0e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "feces": {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "muscle": {"diet": 0.0,"plasma": t*0.0035e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "bone": {"diet": 0.0, "plasma": flux_bone, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "skin": {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0},
            "kidney": {"diet": 0.0,"plasma": 0.5*k*flux_diet, "RBC": 0e0, "liver": 0e0, "urine": 0.125*k*flux_diet, "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney":0e0}}
        logger.debug(pformat(self.Flux))
    
    def parameters(self):
        """ Define the time parameters, the isotopic standard and the boxes, flux and partition coefficients """
        n_timestep = 100
        self.time = linspace(0, 18250.0, n_timestep) # temps
         
        # JMC standard
        self.standard = 0.565203
         
        self.Boxes = {
            "diet": {'Delta': 0e0, 'Mass': 1e12},
            "plasma": {'Delta': 0e0, 'Mass': 3e0},
            "RBC": {'Delta': 0e0, 'Mass': 2.5e1},
            "liver": {'Delta': 0e0, 'Mass': 1.3e2},
            "urine": {'Delta': 0e0, 'Mass': 1e1},
            "feces": {'Delta': 0e0, 'Mass': 1e1},
            "muscle": {'Delta': 0e0, 'Mass': 1.5e3},
            "bone": {'Delta': 0.48e0, 'Mass': 7.7e2},
            "skin": {'Delta': 0e0, 'Mass': 1.6e2},
            "kidney": {'Delta': 0e0, 'Mass': 2e1}
        }
        
        #coeff_KU=1/0.9993e0;
        coeff_KU=1/0.9998e0;
        coeff_PRBC=1.00025e0
        coeff_PS=1.000275;
        coeff_PM=0.99986;
        coeff_PB=1.0003;
        coeff_PL=0.99939;
        #coeff_PD=1.00025;
        coeff_PD=1.000;

        self.Partcoeff = {
            "diet": {"diet": 1.0, "plasma": 1.00018e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone": 1e0, "skin": 1e0, "kidney":1e0},
            "plasma": {"diet": 1.0, "plasma": 1e0, "RBC": coeff_PRBC, "liver": coeff_PL, "urine": 1e0, "feces": coeff_PD, "muscle": coeff_PM, "bone":coeff_PB, "skin": coeff_PS, "kidney":1/coeff_KU},
            "RBC": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1/coeff_PRBC, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "liver": {"diet": 1.0, "plasma": 1/coeff_PL, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "urine": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "feces": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "muscle": {"diet": 1.0, "plasma": 1/coeff_PM, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "bone": {"diet": 1.0, "plasma": 1/coeff_PB, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "skin": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0},
            "kidney": {"diet": 1.0, "plasma": coeff_KU, "RBC": 1e0, "liver": 1e0, "urine": coeff_KU, "feces": 1e0, "muscle": 1e0, "bone":1e0, "skin": 1e0, "kidney":1e0}
            }
            

    def peace_flag(self, boxes, param1, param2 ):
        """ Create a contour plot of the boxes final value as a function of two parameters"""
        logger.info('Drawing contour plot of '+','.join( [ set_style( box, 'log_header') for box in boxes] )+\
                    ' final Delta for '+set_style(param1, 'emph')+' and '+set_style(param2, 'emph')+' ...')
        i_box = [ self.Boxes.keys().index(box) for box in boxes]
        
        x = self.parameters[param1]
        y = self.parameters[param2]
        
        Delta = []
        for yy in y:
            Delta.append([])
            for xx in x:    
                infile = self.result_dir +slugify([param1, xx, param2, yy])+'/Delta.final'
                
                f = open(infile)
                for i, line in enumerate(f):
                    if i in i_box:
                        vals = line.split(' ')
                        Delta[y.index(yy)].append( vals[1].rstrip() )
                f.close()
                
        plt.contourf(x, y, Delta)
        plt.xlabel(param1)
        plt.ylabel(param2)

        plt.colorbar(fraction=0.1)

#        plt.show()
        plt.savefig('test.png')

        
        
        