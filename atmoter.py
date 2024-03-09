import numpy as np

import argparse


from src.thermodynamicData import ThermodynamicData
from src.atmoSpec import Atmosphere


if __name__ == "__main__":

     parser = argparse.ArgumentParser(prog="atmoter", description="Spéciation atmosphère")
     parser.add_argument('cmd', help='Commande')
     parser.add_argument('--ho'   , type=float, help='Valeur du rapport H/O')
     parser.add_argument('--hc'   , type=float, help='Valeur du rapport H/C')
     parser.add_argument('--hn'   , type=float, help='Valeur du rapport H/N')
     parser.add_argument('--t'    , type=float, help='Valeur de la température T en Kelvin')
     parser.add_argument('--p'    , type=float, help='Valeur de la pression P en bar')
     parser.add_argument('--algo' , type=int, choices=[1, 2], help='Algorithme (1 ou 2)')
     parser.add_argument('--ch4'  , type=float, help='Coefficient stoechiométrique du CH4 dans la réaction')
     parser.add_argument('--co2'  , type=float, help='Coefficient stoechiométrique du CO2 dans la réaction')
     parser.add_argument('--co'   , type=float, help='Coefficient stoechiométrique du CO dans la réaction')
     parser.add_argument('--h2o'  , type=float, help='Coefficient stoechiométrique du H2O dans la réaction')
     parser.add_argument('--h2'   , type=float, help='Coefficient stoechiométrique du H2 dans la réaction')
     parser.add_argument('--o2'   , type=float, help='Coefficient stoechiométrique du O2 dans la réaction')
     parser.add_argument('--nh3'  , type=float, help='Coefficient stoechiométrique du NH3 dans la réaction')
     parser.add_argument('--n2'   , type=float, help='Coefficient stoechiométrique du N2 dans la réaction')
     parser.add_argument('--o'    , type=float, help='Coefficient stoechiométrique du O dans la réaction')
     
     args = parser.parse_args()
     
     print ("Atmoter")

     atmosphere = Atmosphere()
     atmosphere.problem = 21

     data = ThermodynamicData()

     T  = args.t  if args.t  is not None else 2000
     P  = args.p  if args.p  is not None else 100
     ho = args.ho if args.ho is not None else 4
     hc = args.ho if args.hc is not None else 0
     hn = args.ho if args.hn is not None else 0

     atmosphere.algorithm = args.algo if args.algo is not None else 1 
     
         #---------------------------------------------------------------
         # Affichage des propriétés JANAF pour la substance sélectionnée
         #---------------------------------------------------------------

     if args.cmd == 'janaf':
        data.printJANAFData (T)

         #---------------------------------------------------------------
         # Calcul des paramètres thermodynamiques d'une réaction
         #---------------------------------------------------------------

     elif args.cmd == 'reaction':
         
         species = {k: getattr(args, k) for k in ['ch4', 'co2', 'co', 'h2o', 'h2', 'o2', 'nh3', 'n2', 'o'] if getattr(args, k) is not None}

         print ("Reaction")
         data.printReaction(species)
         data.calculateParametersForReaction(species,T)

         #---------------------------------------------------------------
         # Calcul de la fugacité d'O2 (tampon IW)
         #---------------------------------------------------------------

     elif args.cmd == 'fO2':
          atmosphere.problem = 30

          atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
          atmosphere.calculateReactionConstants()
          print (f"T              (K)  : {atmosphere.T:10.4g}")
          print (f"Fugacité d'O2  (bar): {atmosphere.fO2:10.4g}")


         #---------------------------------------------------------------
         # Calcul d'un point de spéciation pour une valeur H/O donnée
         #---------------------------------------------------------------

     elif args.cmd == 'speciation':
          atmosphere.problem = 20

          if atmosphere.algorithm == 2:
              atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
              atmosphere.calculateReactionConstants()
              initials = [0.501, 1e-10]
              atmosphere.calculateSpeciation(initials)
          elif atmosphere.algorithm == 1:     
              atmosphere.calculateMolesNumber(T, P, hc, ho, hn)
              atmosphere.calculateReactionConstants()
              atmosphere.calculateSpeciation(None)

         #---------------------------------------------------------------
         # Calcul de la courbe de speciation
         #---------------------------------------------------------------


     elif args.cmd == 'curve':
          atmosphere.problem = 21     
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

     print()