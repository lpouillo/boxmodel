#!/usr/bin/env python
from pprint import pformat
from IsotopicBoxModel import IsotopicBoxModel, sweep, ParamSweeper, slugify, \
    logger, path, mkdir, linspace, style
from numpy import arange


class ZnCompoIsoDiet(IsotopicBoxModel):
    """
    An engine that perform the computation of the evolution of Zn isotopic
    ratio as a function of Diet value
    """
    def __init__(self):
        """Set delta_name"""
        super(ZnCompoIsoDiet, self).__init__()
        self.delta_name = r"$\delta^{66}Zn$"
        self.time = linspace(0, 18250.0, 100000)
        # JMC standard
        self.standard = 0.565203

    def run(self):
        """ Execute the engine and compute the results """
        parameters = {'delta_diet': arange(-0.1, 1.5, 0.5),
                      'coeff_DP': arange(0.9995, 1.0005, 0.0003),
                      'flux_DP': range(8, 17),
                      'flux_PB': arange(0.01, 0.1, 0.03)}
        sweeps = sweep(parameters)
        sweeper = ParamSweeper(path.join(self.result_dir, "sweeps"), sweeps)
        logger.info('Engine will treat %s models',
                    style.emph(len(sweeper.get_remaining())))

        while len(sweeper.get_remaining()) > 0:
            comb = sweeper.get_next()
            logger.info(style.comb(' Performing new combination ') + '\n%s',
                        pformat(comb))
            comb_dir = self.result_dir + '/' + slugify(comb)
            try:
                mkdir(comb_dir)
            except:
                pass
            self.set_boxes(comb['delta_diet'])
            self.set_flux(comb['flux_DP'], comb['flux_PB'])
            self.set_partcoeff(comb['coeff_DP'])
            Delta = []
            Delta = self.initial_state(outdir=comb_dir)
            Delta = self.compute_evolution(Delta, outdir=comb_dir)
            self.final_state(Delta[-1, :], outdir=comb_dir)
            sweeper.done(comb)
            logger.info('Combination done\n')

        logger.info('All combinations have been done, result can be found in '
                    + self.result_dir)

    def set_flux(self, flux_DP, flux_PB):
        """Change value"""
        self.Flux = {
            "diet": {"diet": 0.0, "plasma": flux_DP, "RBC": 0e0,
                     "liver": 0e0, "urine": 0e0, "feces": flux_DP,
                     "muscle": 0e0, "bone": 0e0, "skin": 0e0, "kidney": 0e0},
            "plasma": {"diet": 0.0, "plasma": 0e0, "RBC": 0.18e0,
                       "liver": 2.64e0, "urine": 0e0,
                       "feces": 0.75 * flux_DP, "muscle": 0.0035e0,
                       "bone": flux_PB, "skin": 0.125 * flux_DP,
                       "kidney": 0.625 * flux_DP},
            "RBC": {"diet": 0.0, "plasma": 0.18e0, "RBC": 0e0,
                    "liver": 0.0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0,
                    "bone": 0e0, "skin": 0e0, "kidney": 0e0},
            "liver": {"diet": 0.0, "plasma": 2.64e0, "RBC": 0e0,
                      "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0,
                      "bone": 0e0, "skin": 0e0, "kidney": 0e0},
            "urine": {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0,
                      "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0,
                      "skin": 0e0, "kidney": 0e0},
            "feces": {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0,
                      "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0,
                      "skin": 0e0, "kidney": 0e0},
            "muscle": {"diet": 0.0, "plasma": 0.0035e0, "RBC": 0e0,
                       "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0,
                       "bone": 0e0, "skin": 0e0, "kidney": 0e0},
            "bone": {"diet": 0.0, "plasma": flux_PB, "RBC": 0e0,
                     "liver": 0e0, "urine": 0e0, "feces": 0e0, "muscle": 0e0,
                     "bone": 0e0, "skin": 0e0, "kidney": 0e0},
            "skin": {"diet": 0.0, "plasma": 0e0, "RBC": 0e0, "liver": 0e0,
                     "urine": 0e0, "feces": 0e0, "muscle": 0e0, "bone": 0e0,
                     "skin": 0e0, "kidney": 0e0},
            "kidney": {"diet": 0.0, "plasma": 0.5 * flux_DP, "RBC": 0e0,
                       "liver": 0e0, "urine": 0.125 * flux_DP,
                       "feces": 0e0, "muscle": 0e0, "bone": 0e0, "skin": 0e0,
                       "kidney": 0e0}}

    def set_boxes(self, delta_diet):
        """ Define the time parameters, the isotopic standard and the boxes,
        flux and partition coefficients """
        self.Boxes = {
            "diet":     {'Delta': delta_diet,   'Mass': 1e12},
            "plasma":   {'Delta': 0e0,          'Mass': 3e0},
            "RBC":      {'Delta': 0e0,          'Mass': 2.5e1},
            "liver":    {'Delta': 0e0,          'Mass': 1.3e2},
            "urine":    {'Delta': 0e0,          'Mass': 1e1},
            "feces":    {'Delta': 0e0,          'Mass': 1e1},
            "muscle":   {'Delta': 0e0,          'Mass': 1.5e3},
            "bone":     {'Delta': 0.48e0,       'Mass': 7.7e2},
            "skin":     {'Delta': 0e0,          'Mass': 1.6e2},
            "kidney":   {'Delta': 0e0,          'Mass': 2e1}
        }

    def set_partcoeff(self, coeff_DP):
        """Change the value of the Partition Coefficient between Diet
        and Plasma"""
        coeff_KU = 1 / 0.9998e0
        coeff_PRBC = 1.00025e0
        coeff_PS = 1.000275
        coeff_PM = 0.99986
        coeff_PB = 1.0003
        coeff_PL = 0.99939
        coeff_PD = 1.000

        self.Partcoeff = {
            "diet": {"diet": 1.0, "plasma": coeff_DP, "RBC": 1e0, "liver": 1e0,
                     "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone": 1e0,
                     "skin": 1e0, "kidney": 1e0},
            "plasma": {"diet": 1.0, "plasma": 1e0, "RBC": coeff_PRBC,
                       "liver": coeff_PL, "urine": 1e0, "feces": coeff_PD,
                       "muscle": coeff_PM, "bone": coeff_PB, "skin": coeff_PS,
                       "kidney": 1 / coeff_KU},
            "RBC": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0,
                    "liver": 1 / coeff_PRBC, "urine": 1e0, "feces": 1e0,
                    "muscle": 1e0, "bone": 1e0, "skin": 1e0, "kidney": 1e0},
            "liver": {"diet": 1.0, "plasma": 1 / coeff_PL, "RBC": 1e0,
                      "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0,
                      "bone": 1e0, "skin": 1e0, "kidney": 1e0},
            "urine": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0,
                      "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone": 1e0,
                      "skin": 1e0, "kidney": 1e0},
            "feces": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0,
                      "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone": 1e0,
                      "skin": 1e0, "kidney": 1e0},
            "muscle": {"diet": 1.0, "plasma": 1 / coeff_PM, "RBC": 1e0,
                       "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0,
                       "bone": 1e0, "skin": 1e0, "kidney": 1e0},
            "bone": {"diet": 1.0, "plasma": 1 / coeff_PB, "RBC": 1e0,
                     "liver": 1e0, "urine": 1e0, "feces": 1e0, "muscle": 1e0,
                     "bone": 1e0, "skin": 1e0, "kidney": 1e0},
            "skin": {"diet": 1.0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0,
                     "urine": 1e0, "feces": 1e0, "muscle": 1e0, "bone": 1e0,
                     "skin": 1e0, "kidney": 1e0},
            "kidney": {"diet": 1.0, "plasma": coeff_KU, "RBC": 1e0,
                       "liver": 1e0, "urine": coeff_KU, "feces": 1e0,
                       "muscle": 1e0, "bone": 1e0, "skin": 1e0, "kidney": 1e0}
            }


if __name__ == "__main__":
    engine = ZnCompoIsoDiet()
    engine.start()
