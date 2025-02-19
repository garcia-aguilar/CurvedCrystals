
# -*- coding: utf-8 -*-
"""
Created on May  30th 2017

@author: garcia-aguilar

Geometry object, define attributes, get data from .c output and calculate energies

"""

import pandas as pd
import numpy as np

lattice_a = 5e-5

#________________________
#########################  
  # Geometry object
#########################    
#
class Geo:
    def __init__(self,gname):
        self.type = gname
    
    def defectNames(self):  
        numDefects = len(self.Rsd)                   #always read defects before ge
        nameList = ['d'+str(d+1) for d in range(numDefects)] 
        return nameList
    
    ###########################
        # Get Data
    #########################                     
    def read_GeData(self,fname):  
        #numDefect = len(self.Rsd)               #always read defects before ge
        defNames = self.defectNames()
        names_ge = ['i','x','y','z','area','H2','K','S','W','hx','hy','hz','nx','ny','nz']+defNames
        cols = range(len(names_ge))
        ge = pd.read_csv(fname, sep = '\t', names=names_ge, usecols=cols)    #  problems after update
        #up#ge = pd.read_csv(fname, sep ='\t')
        #up#ge = ge[ge.columns[0:27]]
        ge.columns = names_ge
        self.X = ge.x
        self.Y = ge.y
        self.Z = ge.z
        self.S = ge.S       
        self.W = ge.W
        self.H2 = ge.H2
        self.K = ge.K
        self.area = ge.area
        self.hx = ge.hx
        self.hy = ge.hy
        self.hz = ge.hz
        self.nx = ge.nx
        self.ny = ge.ny
        self.nz = ge.nz
        self.ls = [getattr(ge,d) for d in defNames]
        #self.ls = [ge.d1,ge.d2,ge.d3,ge.d4,ge.d5,ge.d6,ge.d7,ge.d8,ge.d9,ge.d10,ge.d11,ge.d12]        
        self.rho = np.array(self.S+self.K-self.W)
        self.rho_defs = np.array(self.S+self.K)
        self.radii = np.sqrt(self.X*self.X+self.Y*self.Y+self.Z*self.Z)
        
      # Calculate shifted buoyancy at vertex
        Zshift = self.Z.max()
        if self.Z.max() < 0:
            Zshift = -Zshift              
        Zfx = (self.Z-Zshift)*(self.Z-Zshift)*self.nz   
        self.Bu = -Zfx
        
      # Height as a "measure" of volume
        self.height = self.Z.max()-self.Z.min()
      
      # Width in X direction
        self.width = self.X.max()-self.X.min()
      
      # Sum area
        self.totArea = self.area.sum()
        #print 'total area for ',self.type,' is ',self.totArea
      # Sum volume
        v_i = (self.X*self.nx+self.Y*self.ny+self.Z*self.nz)/3.
        vol = v_i*self.area        
        self.totVolume = vol.sum()
        
      # Initial spherical droplet radius
        self.R = np.power(self.totVolume*3./(4*np.pi),1/3.)
        
        #CHECK#print 'total volume for ',self.type,' is ',self.totVolume
        def get_asphericity(self):
            mean_radius = self.radii.mean()
            std_Dev_radii = self.radii.std(ddof=0)
            asp = np.power(std_Dev_radii/mean_radius,2)
            return asp
        self.asphericity = get_asphericity(self)
        
      # Local curvature integrated in the area
        hi = np.sqrt(ge.H2)*ge.area
        self.localH = hi.sum()
        
    def read_GeData_Disclination(self,fname):  
        names_ge = ['i','x','y','z','area','H2','K','S','hx','hy','hz','nx','ny','nz',
                    'd1','d2','d3','d4','d5','d6','d7','d8','d9','d10','d11','d12']
        cols = range(26)
        ge = pd.read_csv(fname, sep = '\t', names=names_ge, usecols=cols)    #  problems after update
        #up#ge = pd.read_csv(fname, sep ='\t')
        #up#ge = ge[ge.columns[0:27]]
        ge.columns = names_ge
        self.X = ge.x
        self.Y = ge.y
        self.Z = ge.z
        self.S = ge.S 
        self.H2 = ge.H2
        self.K = ge.K
        self.area = ge.area
        self.hx = ge.hx
        self.hy = ge.hy
        self.hz = ge.hz
        self.nx = ge.nx
        self.ny = ge.ny
        self.nz = ge.nz
        self.ls = [ge.d1,ge.d2,ge.d3,ge.d4,ge.d5,ge.d6,ge.d7,ge.d8,ge.d9,ge.d10,ge.d11,ge.d12]        
        self.rho = np.array(self.S+self.K)
        self.radii = np.sqrt(self.X*self.X+self.Y*self.Y+self.Z*self.Z)
        
      # Calculate shifted buoyancy at vertex
        Zshift = self.Z.max()
        if self.Z.max() < 0:
            Zshift = -Zshift              
        Zfx = (self.Z-Zshift)*(self.Z-Zshift)*self.nz   
        self.Bu = -Zfx
        
      # Height as a "measure" of volume
        self.height = self.Z.max()-self.Z.min()
      
      # Sum area
        self.totArea = self.area.sum()
        #print 'total area for ',self.type,' is ',self.totArea
        # Sum volume
        v_i = (self.X*self.nx+self.Y*self.ny+self.Z*self.nz)/3.
        vol = v_i*self.area        
        self.totVolume = vol.sum()
        
      # Initial spherical droplet radius
        self.R = np.power(self.totVolume*3./(4*np.pi),1/3.)
        
        #print 'total volume for ',self.type,' is ',self.totVolume
        def get_asphericity(self):
            mean_radius = self.radii.mean()
            std_Dev_radii = self.radii.std(ddof=0)
            asp = np.power(std_Dev_radii/mean_radius,2)
            return asp
        self.asphericity = get_asphericity(self)
        
    def read_GeData_old(self,fname):
        names_ge = ['i','x','y','z','area','H2','K','S','hx','hy','hz','nx','ny','nz']
        cols = range(14)
        ge = pd.read_csv(fname, sep = '\t', names=names_ge, usecols=cols)    #  
        self.X = ge.x
        self.Y = ge.y
        self.Z = ge.z
        self.S = ge.S       
        self.H2 = ge.H2
        self.K = ge.K
        self.area = ge.area
#        self.hx = ge.hx
#        self.hy = ge.hy
#        self.hz = ge.hz
#        self.nx = ge.nx
#        self.ny = ge.ny
#        self.nz = ge.nz
#        self.rho = np.array(self.S+self.K)
#        self.radii = np.sqrt(self.X*self.X+self.Y*self.Y+self.Z*self.Z)
#        
#        
#        # Calculate shifted buoyancy at vertex
#        Zshift = self.Z.max()
#        if self.Z.max() < 0:
#            Zshift = -Zshift              
#        Zfx = (self.Z-Zshift)*(self.Z-Zshift)*self.nz   
#        self.Bu = -Zfx
#        
#        # Height as a "measure" of volume
#        self.height = self.Z.max()-self.Z.min()
#      
#        # Sum area
#        self.totArea = self.area.sum()
#        #print 'total area for ',self.type,' is ',self.totArea
#        # Sum volume
#        v_i = (self.X*self.nx+self.Y*self.ny+self.Z*self.nz)/3.
#        vol = v_i*self.area        
#        self.totVolume = vol.sum()
#        #print 'total volume for ',self.type,' is ',self.totVolume
#        
#        def get_asphericity(self):
#            mean_radius = self.radii.mean()
#            std_Dev_radii = self.radii.std(ddof=0)
#            asp = np.power(std_Dev_radii/mean_radius,2)
#            return asp
#        self.asphericity = get_asphericity(self)
        
        
    def read_field(self,fname):
        names_last = ['x','y','z','phi','err2']
        last = pd.read_csv(fname, sep = '\t', names=names_last)    #    
        self.phi = last.phi            
        self.phi2 = last.phi*last.phi
        self.err2 = last.err2
        
    def read_triangles_old(self,fname):
        triangles = []
        Dfile = open(fname,"r")
        lines = Dfile.readlines()
        for l in lines:
            triangles.append(np.fromstring(l, sep = '\t'))
        self.triangles = triangles
        
    def read_triangles(self,fname):
        triangles = []
        Dfile = open(fname,"r")     # read tri_*.mshd file
        lines = Dfile.readlines()
        num_vertices = int(lines[1])
        start_at = 7 + num_vertices    #file heading are 5 lines+"N_triangles"+N_t
        num_triangles = int(lines[start_at-1])      
        for l in lines[start_at:]:
            triangles.append(np.fromstring(l, sep = '\t'))
        self.triangles = triangles

    def read_defects(self,DefFile):
        self.Xsd = []
        self.Ysd = []
        self.Zsd = []
        self.Qsd = []
        self.QOff = []
        self.Rsd = []
        self.NearVs = []            
        
        Dfile = open(DefFile,"r")
        lines = Dfile.readlines()
        num_defect = int(lines[3])
        start_at = 6
        for i in range(num_defect):
            props = np.fromstring(lines[start_at+i], sep = '\t')
            self.Xsd.append(props[0])
            self.Ysd.append(props[1])
            self.Zsd.append(props[2])
            self.Qsd.append(props[3])
            self.QOff.append(props[4])
            self.Rsd.append(props[5]) 
            self.NearVs.append(int(props[6]))
    ######## ONLY CST rcore:
        self.core = np.array(self.Rsd).mean()
            
    def read_defects_prePolar(self,DefFile):
        self.Xsd = []
        self.Ysd = []
        self.Zsd = []
        self.Qsd = []
        self.Rsd = []
        self.NearVs = []            
        
        Dfile = open(DefFile,"r")
        lines = Dfile.readlines()
        num_defect = int(lines[3])
        start_at = 6
        for i in range(num_defect):
            props = np.fromstring(lines[start_at+i], sep = '\t')
            self.Xsd.append(props[0])
            self.Ysd.append(props[1])
            self.Zsd.append(props[2])
            self.Qsd.append(props[3])
            self.Rsd.append(props[4]) 
            self.NearVs.append(int(props[5]))
    ######## ONLY CST rcore:
        self.core = np.array(self.Rsd).mean()
        
    ###########################
        # Calculate VOI
    #########################   
    
    def int_bending(self):
        mult = self.H2*self.area*4.
        integral = mult.sum()
        return integral
        
    def int_stretch(self):
        mult = self.phi*self.phi*self.area
        integral = mult.sum()
        return integral       
       
    def int_buoyancy(self):        
        mult = self.Bu*self.area
        integral = mult.sum()
        return integral
    
    def entropy_term(self):
        bDisloc = getattr(self,"b",0.)
        entropy = bDisloc*np.log(self.totArea/lattice_a**2)
        return entropy
        
    def int_buoyancy_noShift(self):
        Zfx = self.Z*self.Z*self.nz   
        self.Bu = -Zfx
        mult = self.Bu*self.area
        integral = mult.sum()
        return integral
    
    def sum_phi(self):
        prod = self.phi*self.area
        return prod.sum()
        
        
    ###########################
        # Assign Energies
    #########################   
    def dim_stretch(self):
        s = self.int_stretch()/np.power(self.totVolume,2./3.)
        return s
    def dim_buoyancy(self):
        b = self.int_buoyancy()/np.power(self.totVolume,4./3.)
        return b
    def dim_area(self):
        return self.totArea/np.power(self.totVolume,2./3.)   
        
    def assign_energies(self):
        self.Es = self.int_stretch()
        self.Eb = self.int_bending()
        self.Ebuo =self.int_buoyancy()
    
    def assign_entropy(self):
        self.entropy = self.entropy_term()

        
    ###########################
        # Energy scaling
    #########################     
    # Elastic energy
    def elastic_energy(self,fvk):
        energy = self.int_bending()+fvk*self.dim_stretch()
        return energy
    # Calculate total energies
    def energy_in_size(self,Fp,gamp,SqBp, cubicRootVolume):
        """ total energy (2E/B), as a function of constant def parameters and scaling volume"""          
        #sl = np.power(volume,2/3.)     
        sl = cubicRootVolume*cubicRootVolume
        E = (Fp*sl*self.dim_stretch() +
             self.int_bending() +
             gamp*sl*self.dim_area() +
             SqBp*sl*sl*self.dim_buoyancy())
        return E 
    def energy_in_size2(self,Fp,Eop,SqBp, cubicRootVolume):
        """ total energy (2E/B), as a function of constant def parameters and scaling volume"""          
        #sl = np.power(volume,2/3.) 
        sl = cubicRootVolume*cubicRootVolume
        E = (Fp*sl*self.dim_stretch() +
             self.int_bending() +
             SqBp*Eop*sl*self.dim_area() +
             SqBp*sl*sl*self.dim_buoyancy())
        return E        
    def energy_in_param(self, FvK, gamma, SqBeta):
        """ total energy (2E/B), as a function of def scaling parameters"""        
        E = (FvK*self.dim_stretch() +
             self.int_bending() +
             gamma*self.dim_area() +
             SqBeta*self.dim_buoyancy())
        return E   
    def energy_in_param2(self, FvK, EoI, SqBeta):
        """ total energy (2E/B), as a function of def scaling parameters"""        
        E = (FvK*self.dim_stretch() +
             self.int_bending() +
             SqBeta*EoI*self.dim_area() +
             SqBeta*self.dim_buoyancy())
        return E      
    def energy_in_netParam(self,Ymod,Bmod,gamma,weight, volume):
        """ total energy (2E/B), as a function of constant charact parameters and scaling volume""" 
        Ymod = 1500.
        Bmod = 1.
        gamma = -100.
        weight = 1000.        
        sl = np.power(volume,2/3.)
        E = (Ymod*sl*self.dim_stretch() +
             Bmod*self.int_bending() +
             gamma*sl*self.dim_area() +
             weight*sl*sl*self.dim_buoyancy())
        return E        
    # Constant volume scaling
    def scaleTo(self,refGeo):
        volRatio = refGeo.totVolume/self.totVolume
        return volRatio
        
#________________________
#########################  
  # Functions / Geo's
#########################    
#    

def critical_FvK(geo1,geo2,fArray):
    volRatio12 = geo2.scaleTo(geo1)        
    Es1 = geo1.elastic_energy(fArray)
    Es2 = geo2.elastic_energy(fArray*np.power(volRatio12,2./3))
    FvK_critical_ix = np.argwhere(np.diff(np.sign(Es1-Es2))!=0)+1
    FvK_critical = fArray[FvK_critical_ix][0][0]
    if not FvK_critical:
        FvK_critical = 0    
        print "\n\t*** No energy intersection found in the given range"
    return FvK_critical
    
def an_crit_FvK(geo1,geo2):
    dBen = geo2.int_bending()-geo1.int_bending()
    dStr = geo1.dim_stretch()-geo2.dim_stretch()
    return dBen/dStr
    
def get_LMN(fvk):
    fact = np.power(4*np.pi/3.,2/3.)
    return fvk/fact

def crit_dimGamma(geo1,geo2,fvk,SqBeta):
    dS = geo1.dim_stretch()-geo2.dim_stretch()
    dH = geo1.int_bending()-geo2.int_bending()
    dZ = geo1.dim_buoyancy()-geo2.dim_buoyancy()
    dArea = geo2.dim_area()-geo1.dim_area()
    return (fvk*dS+dH+SqBeta*dZ)/dArea
def crit_gammaP(geo1,geo2,fvkp,SqBetap,radius):
    dS = geo1.dim_stretch()-geo2.dim_stretch()
    dH = geo1.int_bending()-geo2.int_bending()
    dZ = geo1.dim_buoyancy()-geo2.dim_buoyancy()
    dArea = geo2.dim_area()-geo1.dim_area()
    return (fvkp*dS+dH/(radius*radius)+SqBetap*dZ*radius*radius)/dArea
def crit_beta(geo1,geo2,fvk,gamma):
    dS = geo1.dim_stretch()-geo2.dim_stretch()
    dH = geo1.int_bending()-geo2.int_bending()    
    dArea = geo1.dim_area()-geo2.dim_area()
    dZ = geo2.dim_buoyancy()-geo1.dim_buoyancy()
    return (fvk*dS+dH+gamma*dArea)/dZ     
def crit_Eo_cstF(geo1,geo2):
    dZ = geo1.dim_buoyancy()-geo2.dim_buoyancy()
    dArea = geo2.dim_area()-geo1.dim_area()
    return dZ/dArea
def crit_gamma_cstF(geo1,geo2,SqBeta):
    dZ = geo1.dim_buoyancy()-geo2.dim_buoyancy()
    dArea = geo2.dim_area()-geo1.dim_area()
    return SqBeta*dZ/dArea
