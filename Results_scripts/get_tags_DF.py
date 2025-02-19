#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:17:28 2019

@author: garcia-aguilar

*** Add screening dislocation cloud as a change in the effective disclination charge of disclination of vertices ***
Add opposite sign effective dislocations cores at the center of the faces of the icosahedron
Total charge consistent with topological constrains

** Read all results in resDir specified, get data, calculate analytical stretch and assemble into energy.dat file **

Generated mesh files
/data1/Garcia/OilDroplets/Calc/Defects/Dislocations/OnFaces/Geo/

ATM, r0 read as input from inputParameters.txt, but still written on that file to the assumed r0 value.

Fix r0 = 0.04

include geo class file 

"""


import pandas as pd
import os
import re
import numpy as np
#from mayavi import mlab
import imp
import geo_obj as geo

#geo = imp.load_source('geo_obj', '/data1/Garcia/OilDroplets/Code/PyScripts/geo_obj.py')

##########################################################################
####    Read Data       ###############################################
##########################################################################


mainDir = '../' 
geoDir = mainDir+'../Geo/'
resDir = mainDir+'Output_examples/'   #directory with results from numerical integration
outDir = mainDir+'Results_scripts/Output/'  

energiesOutFile = 'energyDF_r0.04'
tagsOutFile = 'current_tags_r0.04'

fixed_r0 = 0.04

###################
fileName = 'last'

def get_tags_in_dir(directory,fileName,tagFormat):
    tag_structure = re.compile(fileName+tagFormat)          
    tags = []
    paramList = []   
    for someFile in os.listdir(directory):
        match = tag_structure.match(someFile) 
        if match:
            tag = match.group()[len(fileName):-len('.dat')]
            #print match.group()
            matchTuple = match.groups()
            tags.append(tag)            
            paramList.append([int(matchTuple[0]),float(matchTuple[1]),float(matchTuple[2]),tag])        #matches last_i*_R*_b*.dat            
    print "Number of result files for  currModel =", len(paramList)
    return paramList

##########################################################################
####    Get this one       ###############################################
Tag_Pairs = '_i(?P<geo>[0-9]+)_R(?P<R0>[0-9]*\.[0-9]*)_b(?P<b>[0-9]*\.[0-9.]*).*\.dat'            #matches last_i*_R*_b*.dat
Params_Pairs = get_tags_in_dir(resDir, 'last', Tag_Pairs)


## Data variables
allGeosDF = pd.DataFrame()
allGeos = []
allTags = []

##########################################################################
####    Read Data / Calculate Geo Values       ##########################

def getData(ParamList, rDir):           # Returns geos DF, list and tags
    global allGeosDF
    geosDF = pd.DataFrame()
    readGeosList = [] 
    readTagsList = []
    for resIndx,params in enumerate(ParamList):
        r0 = fixed_r0
        m = params[0]
        tag = params[3]         
        c = 0.0
        R0 = params[1]
        b = params[2]   
        
        mshTag = 'ico{0}'.format(m)
        resTag = tag       
        
        dataDict = dict()
        g = geo.Geo(resTag)        
        g.read_defects(rDir+'defs'+resTag+'.dat')        
        g.read_GeData(rDir+'ge'+resTag+'.dat')
        g.read_field(rDir+'last'+resTag+'.dat')
        g.read_triangles(geoDir+'tri_'+mshTag+'_3.mshd')
        g.core = r0
        g.b = b
        g.c = c
        g.R0 = R0       
        
        g.assign_energies()
        g.assign_entropy()        

            
        # Add DataFrame info
        dataDict['geo'] = g
        dataDict['tag'] = resTag
        dataDict['ico'] = m
        dataDict['R0'] = R0
        dataDict['r0'] = r0
        dataDict['b'] = b
        dataDict['c'] = c
        dataDict['asph'] = g.asphericity
        dataDict['Es'] = g.Es
        dataDict['dimEs'] = g.dim_stretch()
        dataDict['Eb'] = g.Eb
        dataDict['A'] = g.totArea
        dataDict['dimA'] = g.dim_area()
        dataDict['Ebuo'] = g.Ebuo
        dataDict['dimEbuo'] = g.dim_buoyancy()
        dataDict['entropy'] = g.entropy
        dataDict['V'] = g.totVolume   
        dataDict['sVert'] = g.Qsd[0]
        dataDict['sFace'] = g.Qsd[-1]      
        dataDict['intLocalH'] = g.localH
        frame = pd.DataFrame.from_dict(data = [dataDict], orient='columns')
        geosDF = pd.concat([geosDF,frame], axis = 0 )
        
        readGeosList.append(g)
        readTagsList.append(resTag)
        
        allGeosDF = pd.concat([allGeosDF,frame], axis = 0) 
        allGeos.append(g)   
        allTags.append(resTag) 
        
    print "Data read"
    return geosDF, readGeosList, readTagsList

##############################################################################
##############################################################################
    
someDF, geosSome, tagsSome = getData(Params_Pairs,resDir)

allGeosDF.index = range(len(allGeosDF))    
mIndxList = allGeosDF.groupby('ico').groups.keys()


##########################################################################
####    Check Data       ###############################################
##########################################################################
# CHECK
print "A total of",len(allGeos)," results read"

# Split into lists per ico#
burgerListsPerGeos = []
listsPerGeos = []
asphericities = []
for icoNum in mIndxList:        
    sorted_DF = allGeosDF.groupby('ico').get_group(icoNum).sort_values(by='b')
    sorted_DF = sorted_DF.reset_index()    
    bSerie = sorted_DF.b
    burgerListsPerGeos.append(bSerie)
    listsPerGeos.append(sorted_DF.geo)
    asphericities.append(sorted_DF.iloc[0].asph)

print "Asphericities = ",asphericities    
###END  TODO !!!!!!!!!!!!!! checking #

##########################################################################
####    Save Data       ###############################################
##########################################################################

## Save Geos energy information
def save_each_Size(listRs):
    for R0 in listRs:
        partialDF = allGeosDF.groupby('R0').get_group(R0).copy()
        partialDF.to_csv(outDir+energiesOutFile+'_R{0}.dat'.format(R0), sep='\t')
    

##allDF = pd.concat([pairsDF,checkDF],axis=0)
allGeosDF.to_csv(outDir+energiesOutFile+'.dat',sep='\t')
 
## Save list of read tags
with open(outDir+tagsOutFile+'.txt', 'w') as tagFile:
    for tag in allTags:
        tagFile.write(tag+"\n")
        
#save_each_Size([0.25,0.5,0.75,1.0,2.0])