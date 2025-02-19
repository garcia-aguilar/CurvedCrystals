#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 10:40:28 2018

@author: garcia-aguilar

*** Add screening dislocation cloud as a change in the effective disclination charge of disclination of vertices ***
Add opposite sign effective dislocations cores at the center of the faces of the icosahedron
Total charge consistent with topological constrains


** Look at stretching energy behaviour of the different Bs **


Generated mesh filesin: 
../Geo/

ATM, r0 read as input from inputParameters.txt, but still written on that file to the assumed r0 value.

Fields on Geos on plotFields.py


Fix r0 = 0.04
include geo class file 


"""


import pandas as pd
import os
import re
import numpy as np
#from mayavi import mlab
import imp
Import geo_obj as geo

#geo = imp.load_source('geo_obj', '/data1/Garcia/OilDroplets/Code/PyScripts/geo_obj.py')
#PlotMayavi = imp.load_source('Geo_plotting_with_Mayavi','/data1/Garcia/OilDroplets/Code/PyScripts/Plot/Geo_plotting_with_Mayavi.py')

##########################################################################
####   Specify  Data       ###############################################
##########################################################################

mainDir = '../'
geoDir = mainDir+'Geo/'
resDir = mainDir+'Output_examples/'
outDir = mainDir+'Results_scripts/Output/'  

energiesInFile = outDir+'energyDF_ALL_withB.dat'

fixed_r0 = 0.04
dist_method = 3


##########################################################################
####   Read Data       ###############################################
##########################################################################

## get Geos energy information
geosDF = pd.read_csv(energiesInFile,sep='\t') 
geosDF = geosDF.drop('geo',1)        # delete column with "geo object"


##########################################################################
####    Stretching energy in Param     ###############################################
########################################################################## 

##########################################################################
####   FUNCTIONS Minima in parameters    ################################

def iMin_for_iR(dataFrame, ico,R0):             # b in EsMin
    subDF = dataFrame.groupby(['R0','ico']).get_group((R0,ico)).copy()
    indx_min = subDF.Es.argmin()
    return indx_min

def iMin_EFree_for_iR(dataFrame, ico,R0):             # b in EFreeMin = Etot-TS
    subDF = dataFrame.groupby(['R0','ico']).get_group((R0,ico)).copy()
    indx_min = subDF.EFree.argmin()
    return indx_min

def find_EsMin_idxs_for_i(dataFrame,ico):
    list_ixs = []
    sortedDF = dataFrame.sort_values(by='R0').copy()
    listR0s = np.sort(sortedDF.groupby('R0').groups.keys())
    for R0 in listR0s:
        idx = iMin_for_iR(sortedDF,ico,R0)        
        list_ixs.append(idx)
    return list_ixs

def find_EFreeMin_idxs_for_i(dataFrame,ico):
    list_ixs = []
    sortedDF = dataFrame.sort_values(by='R0').copy()
    listR0s = np.sort(sortedDF.groupby('R0').groups.keys())
    for R0 in listR0s:
        idx = iMin_EFree_for_iR(sortedDF,ico,R0)        
        list_ixs.append(idx)
    return list_ixs

##########################################################################
####   Find average fitting dim coefficients from a DF    ##############
    
def get_avgFitCoeff_perIco(dataFrame,degree):       
    '''return [ico,c2,c1,c0]'''
    
    coefficients = {}    
    listIcos = np.sort(dataFrame.groupby('ico').groups.keys())    
    for ico in listIcos:
        subDF = dataFrame.groupby(['ico']).get_group(ico).copy()       
        listR0s = np.sort(subDF.groupby('R0').groups.keys())
        listC1s = []
        listC2s = []
        listC0s = []
        for R0 in listR0s:
            CDF = subDF.groupby(['R0']).get_group(R0).copy()    
            y = CDF.Es
            x = CDF.b            
            fit = np.polyfit(x, y, degree)  
            c1 = -fit[degree-1]/R0
            c2 = fit[degree-2]
            c0 = fit[degree]/R0/R0
            listC1s.append(c1)
            listC2s.append(c2)
            listC0s.append(c0)
        
        c1 = np.average(np.array(listC1s))
        c2 = np.average(np.array(listC2s))
        c0 = np.average(np.array(listC0s))
        coefficients.update({ico:[c2,c1,c0]})
    
    return coefficients
        
avgCoefficients_all = get_avgFitCoeff_perIco(geosDF,2)  

##########################################################################
####   Filter DataSet    ################################
# Not consider entries with the exact same Stretching energy (likely repeats of b=0 params or Check/Main)
#currDF = geosDF.drop_duplicates(subset=['Es'])

# Stick to sphere of radius R0 to study different configurations
#currDF = currDF.groupby(['R0','ico']).get_group((1.0,12)).copy()
    
#def exclude_some()

##########################################################################
####   PLOTS    #################
print '\n***------***\n'
import matplotlib.pyplot as plt

baseFont = 14

##########################################################################
####   Plot Es(B) for a specific size R, ico as key     #####################
    
def plot_Es_vsB_inIco(subDF,R0):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())

    plt.clf()
    plt.ylabel('Stretching E', fontsize=baseFont)
    plt.xlabel('B', fontsize=baseFont)
    plt.title('plot all results for R0={0} with ico as key'.format(R0),fontsize=baseFont+2)
    
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))
        CDF = subDF.groupby(['ico','R0']).get_group((currIco,R0)).copy()
        plt.plot(CDF.b,CDF.Es,'o',color =colorR,label = 'i='+str(currIco))
      
        # Plot b = 0 line
        SerieB = CDF.groupby('b').get_group(0.0).copy().iloc[0]
        EsB0 = SerieB.Es
        
        plt.axhline(y=EsB0,linestyle='--',color=colorR)    
        
    plt.legend(fontsize=baseFont)
    

#plt.figure()

df1= geosDF.groupby('ico').get_group(1).copy()
df5= geosDF.groupby('ico').get_group(5).copy()
df12= geosDF.groupby('ico').get_group(12).copy()
dficos = pd.concat([df1,df5])
currDF = pd.concat([df12,df5])
#plot_Es_vsB_inIco(currDF,0.75)    
#plt.figure()

##########################################################################
####    Es(B) polynomail fit R=1, ico as key     #####################

## np.polyfit, return p[] --- f(x) = p[0] x**d + p[1] x**d-1 +....+ p[d] 

def fit_Es_vsB(subDF, R0, ico, degree):
    #listIcos = np.sort(subDF.groupby('ico').groups.keys())
    #listIcos = [12]
    
    plt.clf()
    #plt.ylabel('Stretching E', fontsize=baseFont)
    #plt.xlabel('B', fontsize=baseFont)
    #plt.title('plot results for R0={0} and ico {1} with fit to degree {2}'.format(R0,ico,degree),fontsize=baseFont+2)
    plt.title('R={0}'.format(R0))
    #for currIco in listIcos:   
    colorR = str(0.8*((13-ico)/12.))     
    CDF = subDF.groupby(['ico','R0']).get_group((ico,R0)).copy()    
    y = CDF.Es
    x = CDF.b
    fit = np.polyfit(x, y, degree)
    polynomial = np.poly1d(fit)
    
    plt.plot(x, y, 'o',color ='b',label = 'data')
    xmid = (x.max()-x.min())/2
    #plt.plot(xmid,polynomial(xmid),'*',color='r')
    
    xcont = np.linspace(x.min(),x.max(),num=50)
    plt.plot(xcont, polynomial(xcont),'--',color = '0.5')
    #plt.plot(xcont, polynomial(xcont)-0.01*xcont,'--',color = '0.8')
    
    # Plot b = 0 line
    SerieB = CDF.groupby('b').get_group(0.0).copy().iloc[0]
    EsB0 = SerieB.Es
    
    plt.axhline(y=EsB0,linestyle='--',color='0.5')   
    
  # Coefficients:
    print 'Coefficients for i{0}_R{1} to degree {2} are '.format(ico,R0,degree)
    for n in reversed(range(degree+1)):
        print '\t {0} ->  {1}'.format(n,fit[degree-n])    
    return fit
    
#somefit = fit_Es_vsB(geosDF, 0.5, 12, 2)
#plt.text(0.05,0.0018,"${0:.3}\ b^2  {1:.3}\ b + E(b=0)$".format(somefit[0],somefit[1]),fontsize=baseFont-2)
 
    

def plot_fittedEs_vsB_inIco(subDF,R0,degree):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())

    plt.clf()
    plt.ylabel('Stretching E', fontsize=baseFont)
    plt.xlabel('B', fontsize=baseFont)
    plt.title('plot all results for R0={0} with ico as key'.format(R0),fontsize=baseFont+2)
    
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))
        #colorR = '0.1'
        CDF = subDF.groupby(['ico','R0']).get_group((currIco,R0)).copy()
        y = CDF.Es
        x = CDF.b
        plt.plot(x,y,'.',color =colorR,label = 'i='+str(currIco))
        
        # Fitting
        fit = np.polyfit(x, y, degree)
        polynomial = np.poly1d(fit)
        xcont = np.linspace(x.min(),x.max(),num=50)
        plt.plot(xcont, polynomial(xcont),'--',color = colorR)
        # chk
        xmid = (x.max()-x.min())/2
        plt.plot(xmid,polynomial(xmid),'*',color='r')
      
        # Plot b = 0 line
        SerieB = CDF.groupby('b').get_group(0.0).copy().iloc[0]
        EsB0 = SerieB.Es        
        plt.axhline(y=EsB0,linestyle='dotted',color=colorR)    
        
    plt.legend(fontsize=baseFont)
    
plotDF = geosDF.groupby('ico').get_group(12)
#plot_fittedEs_vsB_inIco(plotDF,0.5,2)
    


def plot_C1_fit_vsR(subDF,degree):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel('C1', fontsize=baseFont)
    plt.xlabel('Size R0', fontsize=baseFont)
    plt.title('Coefficient of linear term vs size per ico with fit to degree {0}'.format(degree),fontsize=baseFont+2)
    
    print 'Coefficient for b to degree {0} are'.format(degree)    
    for ico in listIcos:
        colorR = str(0.8*((13-ico)/12.))
        listR0s = np.sort(subDF.groupby('R0').groups.keys())
        listCs = []
        for R0 in listR0s:
            CDF = subDF.groupby(['ico','R0']).get_group((ico,R0)).copy()    
            y = CDF.Es
            x = CDF.b
            fit = np.polyfit(x, y, degree)                 
            listCs.append(fit[degree-1])
            print '\t for i{0}_R{1} -> {2}'.format(ico,R0,fit[degree-1])
            
        plt.plot(listR0s, listCs, 'o--', color=colorR, label='i='+str(ico))
        
    plt.legend(fontsize=baseFont)
            
#plot_C1_fit_vsR(geosDF,2)           

def plot_bmin_fit_vsR(subDF,degree):
#def plot_Esmin_fit_vsR(subDF,degree):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel('bmin', fontsize=baseFont)
    #plt.ylabel('Es Min', fontsize=baseFont)
    plt.xlabel('Size R', fontsize=baseFont)
    plt.title('Smin per ico with and fit to degree {0} '.format(degree),fontsize=baseFont+2)
    
    print '\nCalculated bmin using coefficients to degree {0} - along with set with lowest energy'.format(degree)    
    for ico in listIcos:
        colorR = str(0.8*((13-ico)/12.))
        listR0s = np.sort(subDF.groupby('R0').groups.keys())
        listBmin = []
        listSmin = []
        listC1s = []
        listC2s = []
        listC0s = []
        for R0 in listR0s:
            CDF = subDF.groupby(['ico','R0']).get_group((ico,R0)).copy()    
          # From Fit
            y = CDF.Es
            x = CDF.b            
            fit = np.polyfit(x, y, degree)  
            c1 = -fit[degree-1]/R0
            c2 = fit[degree-2]
            c0 = fit[degree]/R0/R0
            listC1s.append(c1)
            listC2s.append(c2)
            listC0s.append(c0)
            
          # From Data
            Smin = y.min()
            listSmin.append(Smin)
            indx_min = y.argmin()
            bmin = CDF.b.loc[indx_min]
            listBmin.append(bmin)
            
            bmin_seff = 5*np.pi*R0/24.
            
            if ico==12:
                print '\t for i{0}_R{2} -> {1}  (b={3})  [beff={4}]'.format(ico,c1*R0/(2.*c2),R0,bmin,bmin_seff)
            else:
                print '\t for i{0}_R{2} -> {1}  (b={3}) '.format(ico,c1*R0/(2.*c2),R0,bmin,bmin_seff)
    
      # From Fit
        c1 = np.average(np.array(listC1s))
        c2 = np.average(np.array(listC2s))
        c0 = np.average(np.array(listC0s))
        xcont = np.linspace(0,listR0s[-1],num=50)  
        fittedValues = c1*xcont/(2.*c2)         #bmin
        #fittedValues = (c0 - c1*c1/4/c2)*xcont*xcont
        plt.plot(xcont, fittedValues,'--',color = colorR, label='i='+str(ico))              
      # From Data
        calcValues = listBmin   #bmin
        #calcValues = listSmin
        plt.plot(listR0s, calcValues, 'o', color=colorR, label='i='+str(ico))  
        
        print '\nCoefficients to degree {0}'.format(degree)  
        print 'for i{0} ->\n\tC1 = {1}\n\tC2= {2}\n\tC0= {3}'.format(ico,listC1s,listC2s, listC0s)
        print '\tAVGS for i{3}_{4} (2,1,0) = [{0},\t{1},\t{2}]'.format(c2,c1,c0,ico,fixed_r0)
        print "R^2-factor for ico{1} \t {0:.3E}\n\n".format(c0-(c1*c1/4/c2),ico)
    #plt.legend(fontsize=baseFont)
    #plt.text(0.5,1.,"$c_3\ b^3+c_2\ b^2 - c_1\ R\ b + E(b=0)$\n"+r"$bmin=\frac{c_1 R}{2c_2}$",fontsize=baseFont)
    #plt.text(0.25,0.01,"$c_2\ b^2 - c_1\ R\ b + E(b=0)$\n"+r"$Emin=[c_0-\frac{c_1^2 }{4c_2}]R^2$",fontsize=baseFont)


#plt.figure()  
sub = geosDF.groupby('ico').get_group(1)          
#plot_bmin_fit_vsR(sub,2)  


##########################################################################
####    Bmin for Estre+Entropy compare to fit expression, ico as key     ###########



##########################################################################
####    Is the polynomial the "same" for R !=1 ?  YES :) understood from analytics   #####################
    


##########################################################################
####   Plot Bmin (R) ico as key     #####################

def plot_Bmin_vsR_inIco(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel('B at minEs', fontsize=baseFont)
    plt.xlabel('R0', fontsize=baseFont)
    plt.title('B value at min stretching with ico as key',fontsize=baseFont+2)
    
    for currIco in listIcos:
        #colorR = str(0.8*((13-currIco)/12.))        
        Bmins = []
        CDF = subDF.groupby(['ico']).get_group(currIco).copy()
        R0s = np.sort(CDF.groupby('R0').groups.keys())
        for R0 in R0s:
            ixMin = iMin_for_iR(CDF,currIco,R0)
            bmin = CDF.b.loc[ixMin]
            Bmins.append(bmin)
            
        #plt.plot(R0s,Bmins,'o',color =colorR,label = 'i='+str(currIco))
        plt.plot(R0s,Bmins,'o', label = 'i='+str(currIco))
      
    plt.legend(fontsize=baseFont)
    
#plt.figure()
#plot_Bmin_vsR_inIco(currDF)


##########################################################################
##########################################################################
####    Effective density charge       #######################################
########################################################################## 

def plot_Seff_inMin_vsR_inIco(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel('Effective defect charge density at minEs', fontsize=baseFont)
    plt.xlabel('R', fontsize=baseFont)
    plt.title('Seff at vertices (o) and faces (x) at min stretching',fontsize=baseFont+2)
    
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))        
        Seff_vertices = []
        Seff_faces = []
        CDF = subDF.groupby(['ico']).get_group(currIco).copy()
        R0s = np.sort(CDF.groupby('R0').groups.keys())
        for R0 in R0s:
            ixMin = iMin_for_iR(CDF,currIco,R0)
            Svert = CDF.sVert.loc[ixMin]
            Sface = CDF.sFace.loc[ixMin]
            Seff_vertices.append(Svert)
            Seff_faces.append(Sface)
            
        #plt.plot(R0s,Bmins,'o',color =colorR,label = 'i='+str(currIco))
        plt.plot(R0s,Seff_vertices,'o', color=colorR, label = 'i='+str(currIco))
        plt.plot(R0s,Seff_faces, 'x', color=colorR) 
    
    plt.axhline(y=np.pi/8., linestyle='--',color='0.5')
    plt.axhline(y=np.pi/3., linestyle='dotted',color='0.7')
    plt.axhline(y=0, linestyle='dotted',color='0.7')
      
    plt.legend(fontsize=baseFont)    
    
#plt.figure()
#plot_Seff_inMin_vsR_inIco(geosDF)


##########################################################################
##########################################################################
####    dimensionless stretching energy       ###############################################
########################################################################## 

##########################################################################
####   Plot dimEs(R) for bMin, ico as key     #####################
    
def check_dimEs_forB0_vsR_inIco(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel(r'Es/V$^{2/3}$ (only Disclination)', fontsize=baseFont)
    plt.xlabel('R0', fontsize=baseFont)
    plt.title('Dimensionless stretching at B=0 for each geo',fontsize=baseFont+2)
    
    for currIco in listIcos:
        #colorR = str(0.8*((13-currIco)/12.))        

        CDF = subDF.groupby(['ico','b']).get_group((currIco,0.0)).copy()
        CDF = CDF.sort_values(by=['R0'])
        
        plt.plot(CDF.R0, CDF.dimEs, 'o', label='i='+str(currIco))
               
#        R0s = np.sort(CDF.groupby('R0').groups.keys())
#        
#        
#        for R0 in R0s:
#            
#            ixMin = iMin_for_iR(CDF,currIco,R0)
#            dimEs = CDF.dimEs.loc[ixMin]
#            #bmin = CDF.b.loc[ixMin]
#            dimEss.append(dimEs)
#            
#        #plt.plot(R0s,Bmins,'o',color =colorR,label = 'i='+str(currIco))
#        plt.plot(R0s,dimEss,'o', label = 'i='+str(currIco))
#      
    SerieB = subDF.groupby(['b','ico','R0']).get_group((0.0,12,1.0)).copy()
    dimEsB0 = SerieB.dimEs*1.0   
    dimEsB0 = 0.017805          #**************######   
    plt.axhline(y=dimEsB0,xmin=0.,xmax=2.25,linestyle='--',color='0.5')   

    plt.legend(fontsize=baseFont)

#plt.figure()    
#check_dimEs_forB0_vsR_inIco(currDF)    
##########################################################################
####   Plot dimEs(R) for bMin, ico as key     #####################
    
def plot_dimEs_vsR_inIco(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel(r'Es/V$^2/3$ (@Bmin)', fontsize=baseFont)
    plt.xlabel('R0', fontsize=baseFont)
    plt.title('Dimensionless stretching at min B for each geo',fontsize=baseFont+2)
    
    for currIco in listIcos:
        #colorR = str(0.8*((13-currIco)/12.))        
        dimEss = []
        CDF = subDF.groupby(['ico']).get_group(currIco).copy()
        R0s = np.sort(CDF.groupby('R0').groups.keys())
        for R0 in R0s:
            ixMin = iMin_for_iR(CDF,currIco,R0)
            dimEs = CDF.dimEs.loc[ixMin]
            #bmin = CDF.b.loc[ixMin]
            dimEss.append(dimEs)
            
        #plt.plot(R0s,Bmins,'o',color =colorR,label = 'i='+str(currIco))
        plt.plot(R0s,dimEss,'o', label = 'i='+str(currIco))
      
    SerieB = subDF.groupby(['b','ico']).get_group((0.0,12)).copy().iloc[0]
    EsB0 = SerieB.dimEs    
    plt.axhline(y=EsB0,linestyle='--',color='0.5')   

    plt.legend(fontsize=baseFont)

#plt.figure()
#plot_dimEs_vsR_inIco(currDF)



##########################################################################
####    Compare with analytical sphere     #####################
##########################################################################

#---------------------------------------------------------
# Calculate coefficients for an analytical sphere
c2Analytical = 0.103143
c1Analytical = 0.136819
c0Analytical = 0.0527043

#c2Analytical = 0.0899127
#c1Analytical = 0.121058
#c0Analytical = 0.0407481
#---------------------------------------------------------
colorList = ['r','purple','purple','purple','purple','purple','purple','g','purple','purple','purple','blue']

def plot_fittedEs_vsB_withSphere(subDF,R0,degree):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())

    plt.clf()
    plt.ylabel('Stretching E', fontsize=baseFont)
    plt.xlabel('B', fontsize=baseFont)
    plt.title('plot energy for fitted coeffs for R0={0} vs analSphere'.format(R0),fontsize=baseFont+2)
    
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))
        CDF = subDF.groupby(['ico','R0']).get_group((currIco,R0)).copy()
        y = CDF.Es
        x = CDF.b
        plt.plot(x,y,'o',color =colorR,label = 'i='+str(currIco))
        
        # Fitting
        fit = np.polyfit(x, y, degree)
        polynomial = np.poly1d(fit)
        xcont = np.linspace(x.min(),x.max(),num=50)
        plt.plot(xcont, polynomial(xcont),'--',color = colorR)
        # chk
        xmid = (x.max()-x.min())/2
        plt.plot(xmid,polynomial(xmid),'*',color='r')
        
        # Avg Fitting
        if also_avg:
            c2Avg = avgCoefficients_all[12][0] 
            c1Avg = avgCoefficients_all[12][1]
            c0Avg = avgCoefficients_all[12][2]
            Es_avg = c2Avg*xcont*xcont-c1Avg*R0*xcont+c0Avg*R0*R0
            plt.plot(xcont,Es_avg,linestyle='dotted',color = 'g',label='avg')
        
        # Analytical
        Es_analy = c2Analytical*xcont*xcont-c1Analytical*R0*xcont+c0Analytical*R0*R0
        plt.plot(xcont,Es_analy,linestyle='solid',color = '0.4',label='analytical')
        # translated to found =0 result
        Es_shifted = c2Analytical*xcont*xcont-c1Analytical*R0*xcont+fit[degree]
        plt.plot(xcont,Es_shifted,linestyle='solid',color = '0.8',label='shifted')
      
        # Plot b = 0 line
        SerieB = CDF.groupby('b').get_group(0.0).copy().iloc[0]
        EsB0 = SerieB.Es        
        plt.axhline(y=EsB0,linestyle='dotted',color=colorR)    
        
    plt.legend(fontsize=baseFont)

#also_avg = 1    
#plot_fittedEs_vsB_withSphere(df12,1.0,2)

def plot_fittedBmin_vsR_withSphere(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())

    plt.clf()
    plt.ylabel('Bmin', fontsize=baseFont)
    plt.xlabel('R0', fontsize=baseFont)    
    plt.title('Bmin vs analyticalSphere',fontsize=baseFont+2)
    
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))   
        colorR =colorList[12-currIco]
             
        CDF = subDF.groupby(['ico']).get_group(currIco).copy()
        
        # Numerical 
        Bmins = []
        R0s = np.sort(CDF.groupby('R0').groups.keys())        
        for R0 in R0s:
            ixMin = iMin_for_iR(CDF,currIco,R0)
            bmin = CDF.b.loc[ixMin]
            Bmins.append(bmin)
            
        #plt.plot(R0s,Bmins,'.',color =colorR,label = 'i='+str(currIco))
        #plt.plot(R0s,Bmins,'o', label = 'i='+str(currIco))
                
        Rcont = np.linspace(R0s.min()*0.2,R0s.max(),num=50)        

        # Averages
        c1Avg = avgCoefficients_all[currIco][1]
        c2Avg = avgCoefficients_all[currIco][0]
        B_avg = c1Avg*Rcont/(2*c2Avg)
        #plt.plot(Rcont, B_avg, linestyle='dotted',color='0.4',label='avg') \
        plt.plot(Rcont, B_avg, linestyle='-',color=colorR,label='avg')
       
        # Analytical
        B_analyt = c1Analytical*Rcont/(2*c2Analytical)
        #plt.plot(Rcont, B_analyt, linestyle='solid',color='0.8',label='analytical')
    #plt.legend(fontsize=baseFont)
 
#plot_fittedBmin_vsR_withSphere(geosDF)

def plot_fittedEsmin_vsR_withSphere(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())

    plt.clf()
    plt.ylabel('Es_min', fontsize=baseFont)
    plt.xlabel('R0', fontsize=baseFont)    
    plt.title('Minimum stretching energy in size vs analyticalSphere',fontsize=baseFont+2)
            
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))  
        #colorR =colorList[12-currIco]
              
        CDF = subDF.groupby(['ico']).get_group(currIco).copy()
        
        # Numerical 
        Esmins = []
        R0s = np.sort(CDF.groupby('R0').groups.keys())        
        for R0 in R0s:
            ixMin = iMin_for_iR(CDF,currIco,R0)
            Esmin = CDF.Es.loc[ixMin]
            Esmins.append(Esmin)
            
        plt.plot(R0s,Esmins,'.',color =colorR,label = 'i='+str(currIco))
        #plt.plot(R0s,Bmins,'o', label = 'i='+str(currIco))
                
        Rcont = np.linspace(R0s.min()*0.02,R0s.max(),num=50)        

        # Averages
        c0Avg = avgCoefficients_all[currIco][2]
        c1Avg = avgCoefficients_all[currIco][1]
        c2Avg = avgCoefficients_all[currIco][0]
        Es_avg = (c0Avg-c1Avg*c1Avg/(4*c2Avg))*Rcont*Rcont
        plt.plot(Rcont, Es_avg, linestyle='dotted',color='0.4',label='avg')     
        #plt.plot(Rcont, Es_avg, linestyle='-',color=colorR,label='avg')
       
        # Analytical
        Es_analyt = (c0Analytical-c1Analytical*c1Analytical/(4*c2Analytical))*Rcont*Rcont
        #plt.plot(Rcont, Es_analyt, linestyle='solid',color='0.8',label='analytical')
    plt.legend(fontsize=baseFont)
 
#plot_fittedEsmin_vsR_withSphere(geosDF)
    
 ##########################################################################
##########################################################################
####    Effective density charge       #######################################
########################################################################## 

def plot_Seff_inMin_vsR_inIco(subDF):
    listIcos = np.sort(subDF.groupby('ico').groups.keys())
    
    plt.clf()
    plt.ylabel('Effective defect charge density at minEs', fontsize=baseFont)
    plt.xlabel('R', fontsize=baseFont)
    plt.title('Seff at vertices (o) and faces (x) at min stretching',fontsize=baseFont+2)
    
    for currIco in listIcos:
        colorR = str(0.8*((13-currIco)/12.))        
        Seff_vertices = []
        Seff_faces = []
        CDF = subDF.groupby(['ico']).get_group(currIco).copy()
        R0s = np.sort(CDF.groupby('R0').groups.keys())
        for R0 in R0s:
            ixMin = iMin_for_iR(CDF,currIco,R0)
            Svert = CDF.sVert.loc[ixMin]
            Sface = CDF.sFace.loc[ixMin]
            Seff_vertices.append(Svert)
            Seff_faces.append(Sface)
            
        #plt.plot(R0s,Bmins,'o',color =colorR,label = 'i='+str(currIco))
        plt.plot(R0s,Seff_vertices,'o', color=colorR, label = 'i='+str(currIco))
        plt.plot(R0s,Seff_faces, 'x', color=colorR) 
    
    plt.axhline(y=np.pi/8., linestyle='--',color='0.5')
    plt.axhline(y=np.pi/3., linestyle='dotted',color='0.7')
    plt.axhline(y=0, linestyle='dotted',color='0.7')
    
    # Analytical
    qVAnalytical = np.pi/3-c1Analytical/(2*c2Analytical)
    qFAnalytical = 3*c1Analytical/(10*c2Analytical)
    plt.axhline(y=qVAnalytical, linestyle='solid',color='0.3',label='analyt vertices')
    plt.axhline(y=qFAnalytical, linestyle='solid',color='0.3',label='analyt faces')
    
      
    plt.legend(fontsize=baseFont)    
    
#plt.figure()
#plot_Seff_inMin_vsR_inIco(geosDF)
    

##########################################################################
##########################################################################
##########################################################################
####    Analytical expressions for obtained coefficients      #####################
##########################################################################

#Correct to having negative c_1 and keeping analytical expressions

## Use Shape Class
shape = imp.load_source('shape_class', '/data1/Garcia/OilDroplets/Code/PyScripts/shape_class.py')    
    
    
# Use Dictionary   
def get_stretchingInfo(all_Coefficients):
    stretchDict = {}
    for ico in all_Coefficients:     
        all_Coefficients[ico][1]=-all_Coefficients[ico][1]
        c0 = all_Coefficients[ico][2]
        c1 = all_Coefficients[ico][1]
        c2 = all_Coefficients[ico][0]    
        bmin = -c1/(2*c2)
        dimE = c0-c1*c1/(4*c2)
        qEffVertex = np.pi/3.-bmin
        qEffFace = 3*bmin/5        
        stretchDict.update({ico:[[c2,c1,c0],bmin,dimE,[qEffVertex,qEffFace]]})
    return stretchDict

###############################################################
####    Data     #####################
############################################################
    
numCs = get_avgFitCoeff_perIco(geosDF,2)   #{ico:[c2,-c1,c0]}
anCs_s = [c2Analytical,-c1Analytical,c0Analytical]   

shapeDict = get_stretchingInfo(numCs)


sphere = shape.Shape('sphere',12,shapeDict[12][0])
aSphere = shape.Shape('AnSphere',12,anCs_s)
iSharp = shape.Shape('iSharp',1,shapeDict[1][0])
iRound = shape.Shape('iRound',5,shapeDict[5][0])

###############################################################
####    Effective charge density    #####################
############################################################
    
def plot_qEff_vsIco(listShapes):
    xData = np.arange(len(listShapes))+1
    
  # Numerical
    yDataF = []
    yDataV = []
    for shape in listShapes:
        yDataF.append(shape.qFace)
        yDataV.append(shape.qVertex)
        
  # Analytical
    qF_a = aSphere.qFace
    qV_a = aSphere.qVertex
        
    plt.clf()
    plt.axhline(qF_a,0,4,linestyle='dashed',color='0.7')
    plt.axhline(qV_a,0,4,linestyle='dotted',color='0.7')
    plt.plot(xData,yDataV,'o',label='Sources')
    plt.plot(xData,yDataF,'o',label='Sinks')

    plt.legend()

#plot_qEff_vsIco([sphere,iRound,iSharp])

#... y otros... check Plots/Manuscript/ManuscriptPlots.py

##########################################################################
##########################################################################
 ##########################################################################
##########################################################################
####    Using results from *dat     #######################################
########################################################################## 

def read_coefficients_dat(r0,ico,coefDatFile):
    """ Read averaged coefficients from coefDatFile
       for some core size r0 and for some shape ico.
    Return a list [c2, c1,c0] """
    avgCoefLineFormat = '>>>  AVGS for i'+str(ico)+'_r'+str(r0)+'.*= \[(?P<coeffs>[0-9 , .]*)\]'
    #avgCoefLineFormat = '>>>  AVGS for i'+str(ico)+'_r'+str(r0)+'.* = \[(?P<coeffs>[0-9 .]+)\]'
    
    #\[(?P<c2>[0-9]+\.[0-9]+), (?P<c1>[0-9]+\.[0-9]+), (?P<c0>[0-9]+\.[0-9]+)\]'
    structr = re.compile(avgCoefLineFormat)  
    for line in open(coefDatFile,'r'):
        #print line
        match = structr.match(line)
        if match:
            #print "line matched! *****"
            #print line
            arrayCoeffs = np.fromstring(match.groups()[0], sep = ',')
            arrayCoeffs[1] = -arrayCoeffs[1]
            #arrayCoeffs = []
            #for val in listDists:
                #print val
            return arrayCoeffs
        else:
            pass
    print "No match on file\n\n",coefDatFile
    return np.array([])

cs= read_coefficients_dat(0.04,12,mainDir+'coefficients.dat')       

