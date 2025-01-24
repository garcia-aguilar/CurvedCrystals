
# -*- coding: utf-8 -*-
"""
Created on Feb 12th 2019

@author: garcia-aguilar

Shape object for energy calculations. 
Attributes are lenght-independent terms in the free energy. 

Assign dimensionless quantities and stretching coefficients. 
Otherwise, initialize to zero. 

"""

#import pandas as pd
import numpy as np

lattice_a = 5e-5
avgScarLength = 0.69

#________________________
#########################  
  # Geometry object
#########################    
#
class Shape:
    def __init__(self,gname,ico,list_stretchingCoefs,vol=4.*np.pi/3):
        self.type = gname
        self.geo = ico
      
      # Stretching energy coefficients
        self.c2 =  list_stretchingCoefs[0]
        self.c1 = list_stretchingCoefs[1]
        self.c0 = list_stretchingCoefs[2]
        
      # Energy dimensionless terms 
        self.bminDim = -self.c1/(2*self.c2)
        self.dimEsM = self.c0-self.c1*self.c1/(4*self.c2)
        self.areaRatioToSphere = 0.   # x(alpha) in model
        self.dimA = 0.
        self.hBend = 0.
        self.dimZ = 0.
        self.dimIntH = 0.
        self.V = vol
        
      # Effective charge density
        self.qVertex = np.pi/3.-self.bminDim
        self.qFace = 3*self.bminDim/5
        
     # Other geometry quantities
        self.asph = 0.
        
#        print 'Shape object created with name {0}'.format(self.type)
#        print '\t c2 = {0},   c1 = {1},   c0 = {2}'.format(self.c2,self.c1,self.c0)
#        print 'SET dimesionless: area (dimA), bending (dimH), buoyantTerm (dimBuo), asphericity (asph)\n\n'

    ###########################
        # Assign terms
    #########################                   
    def set_dimA(self,area):
        self.dimA = area
        self.areaRatioToSphere = self.dimA/(4*np.pi)
        self.OmegaDisl = np.log(4*np.pi*self.areaRatioToSphere/(lattice_a*lattice_a))
        self.OmegaScar = np.log(4*np.pi*self.areaRatioToSphere/(lattice_a*avgScarLength))
    def set_dimIntH(self,intH):
        self.dimIntH = intH
    def set_hBend(self,bending):
        self.hBend = bending
    def set_dimBuo(self,Z):
        self.dimZ = Z
    def set_asph(self,asph):
        self.asph = asph
         
    # Assign all above terms
    def assign_dimTerms(self,asphericity,dimArea,dimBending,dimBuoyancy,dimIntHCurvature):
        self.set_dimA(dimArea)
        self.set_hBend(dimBending)
        self.set_dimBuo(dimBuoyancy)
        self.set_asph(asphericity)
        self.set_dimIntH(dimIntHCurvature)
        
    def correct_dimTerms(self,initVol,V0=4.*np.pi/3):
        corrFactor = pow(V0/initVol,1./3)
        self.initVol = initVol
        self.dimEsM = self.dimEsM*pow(corrFactor,2.)
        self.dimA = self.dimA*pow(corrFactor,2.)        
        self.dimIntH = self.dimIntH*pow(corrFactor,-1.)
        self.dimZ = self.dimZ*pow(corrFactor,4.)
