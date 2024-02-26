import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import newton

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from src.thermodynamicData import ThermodynamicData
from src.atmoSpec import Atmosphere

# data = ThermodynamicData()

atmosphere = Atmosphere()
atmosphere.problem = 21
atmosphere.algorithm = 1
T = 2000
P = 100
hc = 0
ho = 4
hn = 0

         # Calcul d'un point de spéciation pour une valeur H/O donnée

if atmosphere.problem == 20:
    if atmosphere.algorithm == 2:
         atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
         atmosphere.calculateReactionConstants()
         initials = [0.501, 1e-10]
         atmosphere.calculateSpeciation(initials)
    elif atmosphere.algorithm == 1:     
         atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
         atmosphere.calculateReactionConstants()
         atmosphere.calculateSpeciation(None)
         
         # Calcul de la courbe de speciation
elif atmosphere.problem == 21:
    if atmosphere.algorithm == 2:
         atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
         atmosphere.calculateReactionConstants()
         initials = [0.1, 0.8]
         atmosphere.calculateSpeciationCurve(1e-9, 40, initials)
    elif atmosphere.algorithm == 1:     
         atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
         atmosphere.calculateReactionConstants()
         atmosphere.calculateSpeciation(None)
         atmosphere.calculateSpeciationCurve(1e-9, 40, None)
    

