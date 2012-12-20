#!/usr/bin/python
''' 
Programme permettant de calculer les delta finaux pour les RBC en fonction des variations
du flux diet-plasma et le fractionnement a l'alimentation
'''

from scipy.integrate import odeint
from random import gauss
from pylab import *

# Definition de la fonction d'evolution
def evol_ratio(ratio,t):
	rationew=np.zeros(ratio.size)
	for ii in range(ratio.size):
		outflux=0;
		influx=0;
		for jj in range(ratio.size):
			# le nouveau ratio est calcule a partir des flux entrants et sortants
			outflux = outflux + flux[ii][jj]/Mass[ii]*coeff[ii][jj]*ratio[ii]
			influx = influx + flux[jj][ii]/Mass[ii]*coeff[jj][ii]*ratio[jj]
		rationew[ii]= influx - outflux
	return rationew;

''' 
Parametres communs a chaque simulation 
'''

Ratio_standard_IRMM=0.0637 # rapport isotopique du standard IRMM
n_timestep=100000
time = np.linspace(0, 13870.0, n_timestep) # temps

# Conditions initiales des boites
Delta={"diet": 1.0e0, "plasma": 1.51e0, "RBC": 2.74e0, "liver": 1.35e0, "urine": 1.0e0, "feces": 0.1e0, "menses": 2.5e0}
Boxes=Delta.keys()
Delta=array(Delta.values())
Ratio=[];
for ii in range(Delta.size):
	Ratio.append((Delta[ii]/1e3+1e0)*Ratio_standard_IRMM)
Ratio=array(Ratio);




# Coefficients de partage
coeff_P_RBC=1.0006202e0;
coeff_P_L=0.9999143e0;
coeff_D_P=0.9990997e0;

data_coeff={
	"diet": {"diet": 1e0, "plasma": 1/coeff_D_P, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses":1e0 },
	"plasma": {"diet": 1e0, "plasma": 1e0, "RBC": coeff_P_RBC, "liver": coeff_P_L, "urine": 1e0, "feces": 1e0, "menses":1e0 },
	"RBC": {"diet": 1e0, "plasma": 1/coeff_P_RBC, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses":1e0 },
	"liver": {"diet": 1e0, "plasma": 1/coeff_P_L, "RBC":1e0, "liver": 1e0, "urine": 1e0, "feces": 1.0e0, "menses":1e0 },
	"urine": {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses":1e0 },
	"feces": {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine": 1e0, "feces": 1e0, "menses":1e0 },
	"menses": {"diet": 1e0, "plasma": 1e0, "RBC": 1e0, "liver": 1e0, "urine":1e0, "feces": 1e0, "menses":1e0 }}
# tranformation en list puis en array des coefficients de partage
list_coeff=[]
for box in data_coeff.values():
	list_coeff.append(box.values()) 
coeff=array(list_coeff);

''' 
Parametres specificiques pour chaque simulation
- variation de la masse du foie
- variation du flux de la diette vers le plasma
''' 

# definition de la gamme de masse
mass_min=0.1
mass_max=1e3
n_mass=20
mass_L=np.linspace(mass_min, mass_max, n_mass)

# definition de la gamme de flux
flux_min=0.0
flux_max=0.7
n_flux=20
flux_DP=np.linspace(flux_min, flux_max, n_flux)

# Creation du tableau pour recevoir les donnees absicce premier, ordonnee dans le range()
delta_RBC_final= [[None] * n_flux for i in range(n_mass)]

## Debut des simulations 
for i_mass in range(mass_L.size):
	# boucle sur le premier parametre
	print mass_L[i_mass]/mass_L[mass_L.size-1]*100
	Mass={"diet": 1e12, "plasma": 3e0, "RBC": 2.5e3, "liver": mass_L[i_mass], "urine": 1e-10, "feces": 1e-0, "menses": 1e-2}
	Mass=array(Mass.values())

	for i_flux in range(flux_DP.size):
		# boucle sur le second parametre
		
		## Ceci sert a reinitialiser les ratios pour ne pas avoir un tableau ratio qui contient tous les resultats des simulations
		Ratio=[];
		for ii in range(Delta.size):
			Ratio.append((Delta[ii]/1e3+1e0)*Ratio_standard_IRMM)
		Ratio=array(Ratio);
		
		# Flux entre les boites
		data_flux={
			"diet": {"diet": 0.0, "plasma":1.3 , "RBC": 0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 },
			"plasma": {"diet": 0.0, "plasma": 0.0, "RBC": 24.4+flux_DP[i_flux], "liver":5.0 , "urine": 0.1, "feces": 0.0, "menses": 0.0 },
			"RBC": {"diet": 0.0, "plasma": 23.9, "RBC": 0.0, "liver": 0.0, "urine": 0.0, "feces": 0.5, "menses": flux_DP[i_flux] },
			"liver": {"diet": 0.0, "plasma": 4.3+flux_DP[i_flux], "RBC": 0.0, "liver": 0.0, "urine": 0.0, "feces": 0.7-flux_DP[i_flux], "menses": 0.0 },
			"urine": {"diet": 0.0, "plasma": 0.0, "RBC": 0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 },
			"feces": {"diet": 0.0, "plasma": 0.0, "RBC": 0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0 },
			"menses": {"diet": 0.0, "plasma": 0.0, "RBC": 0.0, "liver": 0.0, "urine": 0.0, "feces": 0.0, "menses": 0.0  }}
	
		# tranformation en list puis en array
		list_flux=[]
		for box in data_flux.values():
			list_flux.append(box.values())
		flux=array(list_flux);
		
		Ratio=odeint(evol_ratio,Ratio,time)
		
		delta_RBC_final[i_mass][i_flux]=((Ratio[n_timestep-1][6]/Ratio_standard_IRMM)-1.0)*1000;


# creation de la figure
contourf(flux_DP,mass_L,delta_RBC_final)
xlabel(r"Menses $(mg.day^{-1})$")
ylabel(r"$M{}_{liver}$ $(mg)$")

colorbar(fraction=0.1)
show()

