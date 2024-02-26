import numpy as np
from src.thermodynamicData import ThermodynamicData
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from scipy import optimize
import matplotlib.pyplot as plt

#===========================================================================
class Atmosphere:
#===========================================================================
    """ La classe pour faire les calculs de spéciation dans l'atmosphère
    """
    
    HC = 0
    HN = 0
    HO = 0
    P = 0
    T = 0
    algorithm = 2
    
#--------------------------------------------------------
    def __init__(self):
#--------------------------------------------------------
        self.data = ThermodynamicData()
    
#--------------------------------------------------------
    def calculateMolesNumber (self,t,p,hc,ho,hn):
#--------------------------------------------------------

        self.T = t
        self.P = p
        self.V = 1
        self.n0_tot = self.P * 1e5 * self.V / ( self.data.R * self.T )

        if self.problem == 10 :
            self.HC = hc
            self.HO = ho
            self.HN = hn
            ntot = 1000
            self.nH = ntot / (1/self.HC + 1 + 1/self.HO + 1/self.HN )
            self.nO = self.nH / self.HO
            self.nC = self.nH / self.HC
            self.nN = self.nH / self.HN
            self.xH = self.nH / ntot
            self.xC = self.nC / ntot
            self.xO = self.nO / ntot
            self.xN = self.nN / ntot
        elif self.problem == 1 :    
            self.HN = hn
            self.nH = self.n / ( 1 +  1/self.HN )
            self.nN = self.nH / self.HN
            self.xH = self.nH / self.n
            self.xN = self.nN / self.n
        elif self.problem == 20 or self.problem == 21 :        # system O H
            ntot = 1000     # number of atoms
            self.HO = ho
            self.HC = 0
            self.HN = 0
            self.nH = ntot / ( 1 + 1 / self.HO)   
            self.nO = ntot - self.nH  
            self.xH = self.nH / ntot
            self.xO = self.nO / ntot

#--------------------------------------------------------
    def calculateSpeciation(self,initials):
#--------------------------------------------------------

        if self.problem == 20 or self.problem == 21:
            
            if self.algorithm == 2:
                 solution = optimize.root (self.equations2, initials, jac=self.dequations2, method='hybr')
                 print(solution)
                 return solution.success
            elif self.algorithm == 1: 
                 solution = root_scalar ( self.equation2, method='bisect', bracket=[0,1] )
                 print(solution)
                 return solution.converged

#--------------------------------------------------------
    def dequations2( self, variables ):
#--------------------------------------------------------
        xH2O, xO2 = variables
        
        xH2 = 1 - xH2O - xO2

        nH = 2 * ( 1 - xO2 )
        nO =     ( 2 * xO2 + xH2O )
        
        fH2O = self.P * xH2O
        fO2  = self.P * xO2
        fH2  = self.P * ( 1 - xH2O - xO2 )

        eq1 = self.HO * ( 2 * xO2 + xH2O ) - 2 * ( 1 - xO2 )
        eq2 = self.K4 * self.P**0.5 * xO2**0.5 * (1 - xH2O - xO2 ) - xH2O
        
        deq1dxH2O = self.HO 
        deq1dxO2  = self.HO * 2 + 2

        deq2xH2O = - self.K4 * self.P**0.5 * xO2**0.5 - 1
        deq2dxO2  = 0.5 * xO2**(-0.5) * self.K4 * self.P**0.5 * (1 - xH2O - xO2 ) - self.K4 * self.P**0.5 * xO2**0.5 
        
        return np.array( [ [deq1dxH2O, deq1dxO2 ] , [deq2xH2O, deq2dxO2]  ])



#--------------------------------------------------------
    def equations2( self, variables ):
#--------------------------------------------------------
        xH2O, xO2 = variables
        
        xH2 = 1 - xH2O - xO2

        if xH2O <= 0 or xH2 <= 0 or xO2 <= 0 :
            return 1e10,1e10

        nH = 2 * ( 1 - xO2 )
        nO = 2 * xO2         + xH2O 
        
        fH2O = self.P * xH2O
        fO2  = self.P * xO2
        fH2  = self.P * ( 1 - xH2O - xO2 )

        eq1 = self.HO * ( 2 * xO2 + xH2O ) - 2 * ( 1 - xO2 )
        eq2 = self.K4 * self.P**0.5 * xO2**0.5 * (1 - xH2O - xO2 ) - xH2O
        
        self.xH2O = xH2O
        self.xO2  = xO2
        self.xH2  = xH2
        self.nH = nH
        self.nO = nO
        
        self.fH2O = fH2O
        self.fO2  = fO2
        self.fH2  = fH2
        
        print(f"H/O {self.HO:10.3f} H2O {xH2O:10.3f} H2  {xH2:10.3f} O2  {xO2:10.3f} eq1  {eq1:10.3f} eq2  {eq2:10.3f} ")   
        
        return eq1, eq2


#--------------------------------------------------------
    def equation2 ( self, variable ):
#--------------------------------------------------------

        dzeta = float(variable)
        
        if dzeta < 0 or dzeta > 1:
            return 1e10
        
        if self.HO <= 2:
            nh = 100
            no = nh / self.HO
            
            nh2  = nh/2 * (1-dzeta)
            no2  = no/2 - nh/4*dzeta
            nh2o = nh/2 * (1+dzeta)
            
            nt = nh2 + no2 + nh2o
            
            xh2  = nh2 / nt
            xo2  = no2 / nt
            xh2o = nh2o / nt
        else:
            no = 100
            nh = self.HO * no
            
            nh2  = nh/2 - no * dzeta
            no2  = no/2 * (1-dzeta)
            nh2o = no * dzeta   

            nt = nh2 + no2 + nh2o
            
            xh2  = nh2 / nt
            xo2  = no2 / nt
            xh2o = nh2o / nt
        
        
        fH2O = self.P * xh2o
        fO2  = self.P * xo2
        fH2  = self.P * xh2

        eq0 = self.K4 * self.P**0.5 * xo2**0.5 * (1 - xh2o - xo2 ) - xh2o
        
        self.xH2O = xh2o
        self.xO2  = xo2
        self.xH2  = xh2
        self.nH = nh
        self.nO = no
        
        self.fH2O = fH2O
        self.fO2  = fO2
        self.fH2  = fH2
        
        print(f"H/O {self.HO:10.3f} dzeta {dzeta:10.3f}  H2O {xh2o:10.3f} H2  {xh2:10.3f} O2  {xo2:10.3f} eq0  {eq0:10.3f}  ")   
        
        return eq0



#-------------------------------------------------------------------
    def calculateSpeciationCurve(self, homin, homax, initials):
#-------------------------------------------------------------------

         increment = 0.01
         hoArray = []
         xH2OArray = []
         xH2Array = []
         xO2Array = []

         with open('data.txt', 'w') as file:
              file.write("IGOR\n")
              file.write("WAVES\t\t\tho\t\t\t xH2O\t\t\t xO2\t\t\t xH2\n")
              file.write("BEGIN\n")

            # Boucle sur les valeurs de ho
              ho = homin
              while ho <= homax:
                     self.calculateMolesNumber(self.T, self.P, self.HC, ho, self.HN)
                     success = self.calculateSpeciation(initials)
                     
                     if success:            
                         initials = [self.xH2O, self.xO2 ]

                         file.write(f"{ho:10.8g}\t\t\t {self.xH2O:10.8g}\t\t\t {self.xO2:10.8g}\t\t\t {self.xH2:10.8g}\n")

                         hoArray.append (ho)
                         xH2OArray.append (self.xH2O)
                         xO2Array.append (self.xO2)
                         xH2Array.append (self.xH2)
                     else:
                         initials = [0.999, 1e-4]   

             # Écrire les résultats dans le fichier
                     
             # Incrémenter ho
                     
                     ho += increment

              file.write("END")

         plt.figure(figsize=(10, 6))

         # Tracé des courbes
         plt.plot(hoArray, xH2OArray, label='xH2O')
         plt.plot(hoArray, xO2Array, label='xO2')
         plt.plot(hoArray, xH2Array, label='xH2')

         # Configurations supplémentaires
         plt.xlabel('HO')
         plt.ylabel('x')
         # plt.yscale('log')  # Échelle logarithmique sur l'axe des ordonnées
         plt.title('Spéciation en fonction de H/O')
         plt.legend()
         plt.grid(True)

        # Affichage du graphique
         plt.show() 

#--------------------------------------------------------
    def calculateReactionConstants (self):
#--------------------------------------------------------

        self.data.calculate(self.T)
        self.K0 = self.data.K0
        self.K2 = self.data.K2
        self.K3 = self.data.K3
        self.K4 = self.data.K4
        self.K5 = self.data.K5
        self.K6 = self.data.K6
        self.K7 = self.data.K7
        self.K8 = self.data.K8
        self.K9 = self.data.K9

        self.fO2 = self.K0**2
        self.fO = (self.K6*self.fO2)**(1/2)
        self.xO2 = self.fO2 / self.P 
        self.xOx = self.fO  / self.P 

#--------------------------------------------------------
    def calculateFugacities (self, xH2O, xCO, xN2):
#--------------------------------------------------------

     self.xH2O = xH2O
     self.xCO  = xCO
    
     self.fH2O = self.xH2O * self.P 
     self.fCO  = self.xCO  * self.P 
    
     self.fCO2 = self.fCO  * self.fO2**0.5 / self.K2
     self.fCH4 = self.fCO  * self.fH2O**2  / self.fO2**1.5 / self.K3
     self.fH2  = self.fH2O                 / self.fO2**0.5 / self.K4

     self.xCO2 = self.fCO2 / self.P 
     self.xCH4 = self.fCH4 / self.P 
     self.xH2  = self.fH2  / self.P 

     self.xN2 = xN2
     self.fN2 = self.xN2 * self.P 
     self.fNH3 = np.sqrt( self.fN2 / self.K5 * self.fH2O**3 / self.fO2**1.5 )
     self.xNH3 = self.fNH3 / self.P 

#--------------------------------------------------------
    def print(self):
#--------------------------------------------------------
     
     if self.problem == 10:
        print (f"T    {self.T :.2f} K" )
        print (f"P    {self.P:.2f} bar" )
        print (f"HC   {self.HC :.2e}" )
        print (f"HO   {self.HO :.2e}" )
        print (f"HN   {self.HN :.2e}" )
        print (f"n    {self.n  :.2f}" )
        print (f"nH   {self.nH :.2f}" )
        print (f"nO   {self.nO :.2f}" )
        print (f"nC   {self.nC :.2f}" )
        print (f"nN   {self.nN :.2f}" )
        print (f"xH   {self.xH :.2f}" )
        print (f"xO   {self.xO :.2f}" )
        print (f"xC   {self.xC :.2f}" )
        print (f"xN   {self.xN :.2f}" )
        print (f"xO2  {self.xO2  :.2f}" )
        print (f"xOx  {self.xOx  :.2f}" )
        print (f"xCO  {self.xCO  :.2f}" )
        print (f"xCO2 {self.xCO2 :.2f}" )
        print (f"xCH4 {self.xCH4 :.2f}" )
        print (f"xH2O {self.xH2O :.2f}" )
        print (f"xNH3 {self.xNH3 :.2f}" )
        print (f"xN2  {self.xN2  :.2f}" )
        print (f"xH2  {self.xH2  :.2f}" )
        self.xsum = self.xO2 + self.xOx + self.xCO2 + self.xCO + self.xCH4 + self.xH2O + self.xNH3 + self.xN2 + self.xH2
        print (f"xsum  {self.xsum  :.2f}" )
     elif self.problem == 1:
        print (f"T    {self.T :.2f} K" )
        print (f"P    {self.P:.2f} bar" )
        print (f"HN   {self.HN :.2e}" )
        print (f"n    {self.n  :.2f}" )
        print (f"nH   {self.nH :.2f}" )
        print (f"nN   {self.nN :.2f}" )
        print (f"xH   {self.xH :.2f}" )
        print (f"xN   {self.xN :.2f}" )
        print (f"xNH3 {self.xNH3 :.2f}" )
        print (f"xN2  {self.xN2  :.2f}" )
        print (f"xH2  {self.xH2  :.2f}" )
        self.xsum =  self.xNH3 + self.xN2 + self.xH2
        print (f"xsum  {self.xsum  :.2f}" )
        
        
    # #--------------------------------------------------------
    # def equations1( self, variables ):
    # #--------------------------------------------------------
    #      yCH4, yH2, yO2, yN2, yCO = variables

     
    #      nc = yCH4 + yCO 
    #      nh = 4 * yCH4 + 2 * yH2 
    #      no = 2 * yO2 + yCO 
    #      nn = 2 * yN2
    #      ntot = nc + nh + no + nn
    #      xc = nc / ntot
    #      xh = nh / ntot
    #      xo = no / ntot
    #      xn = nn / ntot

     
    #      eq1 = self.xC - xc
    #      eq2 = self.xH - xh
    #      eq3 = self.xO - xo
    #      eq4 = self.xN - xn
    #      eq5 = yCH4 + yH2 + yO2 + yN2 + yCO - 1
     
     
    #      print (f" CH4 {yCH4:10.3g} CO  {yCO:10.3g}  O2  {yO2:10.3g} H2  {yH2:10.3g} N2  {yN2:10.3g} eq1 {eq1:10.3g} eq2 {eq2:10.3g} eq3 {eq3:10.3g} eq4 {eq4:10.3g} eq5 {eq5:10.3g}")
    #      return eq1, eq2, eq3, eq4, eq5
        




# #--------------------------------------------------------
#     def equations2( self, variable ):
# #--------------------------------------------------------
    
#      alpha = variable
    
#      n0_CH4 = self.n0_CH4
#      n0_H2  = self.n0_H2
#      n0_O2  = self.n0_O2
#      n0_CO  = self.n0_CO
#      n0_N2  = self.n0_N2
    
#      n0_CO2 = self.n0_CO2
    
#      n0_tot = self.n0_tot

#      n_CO   = n0_CO * (1-alpha)
#      n_O2   = n0_O2 * (1-alpha/2)
    
#      n_CO2  = n0_CO2 + n0_CO * alpha
    
#      ntot = n0_CH4 + n0_H2 + n_O2 + n_CO + n_CO2 + n0_N2
    
#      xCO  = n_CO / ntot
#      xO2  = n_O2 / ntot
#      xCO2 = n_CO2 / ntot
    
#      fCO  = xCO  * self.P
#      fO2  = xO2  * self.P
#      fCO2 = xCO2 * self.P
    
    
#      self.ntot = ntot
#      self.n_CO2 = n_CO2
#      self.n_CO  = n_CO
#      self.n_O2  = n_O2
#      self.xCO   = xCO
#      self.xCO2  = xCO2 
#      self.xO2   = xO2
#      eq0 = self.K2 * fCO * fO2**0.5 - fCO2
    
#      print (f"alpha: {alpha:10.4f} eq0 {eq0:10.4f}")
    
#      return eq0


        
# #--------------------------------------------------------
#     def equations( self, variables ):
# #--------------------------------------------------------

#      if self.problem == 10:
#         xH2O, xCO, xN2 = variables

#         self.calculateFugacities(xH2O, xCO, xN2)


#         eq7 = self.HC * (    self.xCO  + self.xCO2 + self.xCH4                          ) - (2 * self.xH2 + 2 * self.xH2O + 4 * self.xCH4 + 3 * self.xNH3) 
#         eq8 = self.HN * (2 * self.xN2  + self.xNH3                                      ) - (2 * self.xH2 + 2 * self.xH2O + 4 * self.xCH4 + 3 * self.xNH3) 
#         eq9 = self.HO * (2 * self.xCO2 + self.xCO  + self.xH2O + self.xOx + 2 * self.xO2) - (2 * self.xH2 + 2 * self.xH2O + 4 * self.xCH4 + 3 * self.xNH3) 

#         print (f"eq7  {eq7  :.2g}" )
#         print (f"eq8  {eq8  :.2g}" )
#         print (f"eq9  {eq9  :.2g}" )
    
#         return eq7, eq8, eq9
    
#      elif self.problem == 1:
    
#         xN2, xH2 = variables
        
#         self.xN2 = xN2
#         self.xH2 = xH2
#         self.fN2 = self.xN2 * self.P
#         self.fH2 = self.xH2 * self.P
        
#         self.fNH3 = np.sqrt ( self.K7 * self.fN2 * self.fH2**3 )
#         self.xNH3 = self.fNH3 / self.P
        
#         hn = ( 3 * self.xNH3 + 2 * self.xH2 ) / ( 2 * self.xN2 + self.xNH3 )
        
#         eq1 = 1 - self.xH2 - self.xNH3 - self.xN2
#         eq2 = self.HN - hn

#         print (f"eq1  {eq1  :.2g}" )
#         print (f"eq2  {eq2  :.2g}\n" )
        
#         return eq1, eq2
