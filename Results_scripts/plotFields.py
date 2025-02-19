#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 10:40:28 2018

@author: garcia-aguilar

*** Add screening dislocation cloud as a change in the effective disclination charge of disclination of vertices ***
Add opposite sign effective dislocations cores at the center of the faces of the icosahedron
Total charge consistent with topological constrains

** Plot surface fields, by reading specified tags **


Generated mesh files: 
../Geo/


Tri_files in:
../Geo/

get r0 from Tag or the default r0=0.04

include geo class file 
include file with mesh plotting functions
"""


import pandas as pd
import os
import re
import numpy as np
import imp
import geo_obj as geo

##########################################################################
####   Specify  Data       ###############################################
##########################################################################

mainDir = '../' 
geoDir = mainDir+'../Geo/'
resDir = mainDir+'Output_examples/'   #directory with results from numerical integration
outDir = mainDir+'Results_scripts/Output/'  

energiesInFile = outDir+'energyDF_All_withB.dat'
tagsInFile = outDir+'current_tags_withB.txt'    # this does include the other r0. #Output for get_tags_DF.py
energiesOtherInFile = 'energyDF_r0_All.dat'

default_r0 = 0.04
dist_method = 3

Tag_Res = '_i(?P<geo>[0-9]+)_R(?P<R0>[0-9]*\.[0-9]*)_b(?P<b>[0-9]*\.[0-9.]*)'            #matches last_i*_R*_b*
Tag_r0 = '_r(?P<r>[0-9]*\.[0-9]*)_i(?P<geo>[0-9]+)_R(?P<R0>[0-9]*\.[0-9]*)_b(?P<b>[0-9]*\.[0-9.]*)'            #matches last_r*_i*_R*_b*

tag_structure = re.compile(Tag_Res)
tagr0_structure = re.compile(Tag_r0)

##########################################################################
####   Which tags?       ###############################################

# Ex params for give_param_byHand
R0s= [0.25,0.5,1.0]
b_numbers = [0.,0.02,0.5]


def read_Tags_fromFile(tagFile):
    with  open(tagFile, 'r') as openFile:
        tags = [line.strip() for line in openFile]
    return tags

def give_param_byHand():
    tags = []
    for b in b_numbers:
        for R0 in R0s:
            tag = '_i12_R{0}_b{1}'.format(R0,b)
            tags.append(tag)
    return tags

whichTags = []
tagsFile = tagsInFile
#whichTags = give_param_byHand()
fileTags = read_Tags_fromFile(tagsFile)

##########################################################################
####   Read Data       ###############################################
##########################################################################

## get Geos energy information
geosDF = pd.read_csv(energiesInFile,sep='\t') 
geosDF = geosDF.drop('geo',1)        # delete column with "geo object"

# SPHERE
sp01DF = geosDF.groupby('ico').get_group(12)
#spOthersDF = pd.read_csv(mainDir+energiesOtherInFile,sep='\t') 
#spOthersDF = spOthersDF.drop('geo',1)        # delete column with "geo object"
#spDF = pd.concat([sp01DF,spOthersDF], axis=0)


def read_geo_fromTag(tag, mainResultsDir,geosList):
    
    disclinationOnly = False 
    effSModel = False

    if tag_structure.match(tag):
        match = tag_structure.match(tag)
        print tag,'matches tag_structurePair'
        matchTuple = match.groups()
        r0 = default_r0
        m = int(matchTuple[0])        
        R0 = float(matchTuple[1])
        b = float(matchTuple[2])
        c = 0.0
        effSModel = True
    elif tagr0_structure.match(tag):
        match = tagr0_structure.match(tag)
        print tag,'matches tag_structure other r0'
        matchTuple = match.groups()
        r0 = float(matchTuple[0])
        m = int(matchTuple[1])        
        R0 = float(matchTuple[2])
        b = float(matchTuple[3])
        c = 0.0
        effSModel = True
    else:
        print "no matching tag for  "+tag
        return
    
    mshTag = 'ico{0}_{1}'.format(m,dist_method)
    resTag = tag        
    resultsDir = mainResultsDir
    
    g = geo.Geo(resTag)        
    if not os.path.isfile(resultsDir+'defs'+resTag+'.dat'):            
        print 'results for tag '+tag+" not found in "+resultsDir
        return
    
    if effSModel:
        g.read_defects(resultsDir+'defs'+resTag+'.dat') 
    else:
        g.read_defects_prePolar(resultsDir+'defs'+resTag+'.dat')
        
    if disclinationOnly:
        g.read_GeData_Disclination(resultsDir+'ge'+resTag+'.dat')
    else:
        g.read_GeData(resultsDir+'ge'+resTag+'.dat')
    g.read_field(resultsDir+'last'+resTag+'.dat')
    g.assign_energies()            
    g.read_triangles(geoDir+'tri_'+mshTag+'.mshd')
    g.R0 = R0
    g.core = r0    
    g.b = b
    g.c = c    
    print 'geo read'
    geosList.append(g)
    geosDic[g]=tag
    return g

## Data variables
geosDic = {}


##########################################################################
####   READ       ###############################################

def find_tags_withPartOfTag(tagFormat, tagList):
    newList = []
    tagPiece = re.compile(tagFormat)
    for tag in tagList:
        if tagPiece.search(tag):
            newList.append(tag)
    return newList

allGeos = []
whichTags = fileTags[0:3]
    
nameString = '.*i12.*R1.0.*b0.[.123]$'
nameString = '.*i12_.*R1.0.*b0.0$'
#whichTags = find_tags_withPartOfTag(nameString,fileTags)

#for tag in whichTags:
#    read_geo_fromTag(tag,resDir,allGeos)
    


##########################################################################
####   Test       ###############################################
def find_geo(someTag,geosList):    
    print '\n'
    for g, tag in geosDic.items():
        if tag == someTag:
            print 'Es = ',g.Es
            print 'Type = ',g.type
            return g
        else:
            print someTag,'not found in dictionary'
            print '\t Reading from\n\t\t',resDir
            read_geo_fromTag(someTag,resDir)
            g = geosList[-1]
            print 'Es = ',g.Es
            print 'Type = ',g.type
            return g


##########################################################################
####    Plot Geos  - Test     ###############################################
########################################################################## 

##MinB ATM (from run based on analyt results )
#            #       _i12_R0.25_b0.164
#            #       _i12_R0.5_b0.327
#            #      _i12_R0.75_b0.491
#            #      _i12_R1.0_b0.655
#            #     _i12_R2.0_b1.31 
minEs_12_tags = ['_i12_R0.25_b0.164','_i12_R0.5_b0.327','_i12_R0.75_b0.491','_i12_R1.0_b0.655','_i12_R2.0_b1.31']
minEs_R1_r0_tags = ['_r0.01_i12_R1.0_b0.655','_i12_R1.0_b0.655','_r0.08_i12_R1.0_b0.655']
#readTags = minEs_12_tags+['_i12_R1.0_b0.0','_i12_R1.0_b1.047']
#
theseTags = ['_i5_R1.0_b0.0','_i5_R1.0_b0.176','_i5_R1.0_b0.4']
#theseTags = ['_i12_R0.5_b0.0','_i12_R0.5_b0.327','_i12_R0.5_b0.7']

zeroTags = ['_i1_R1.0_b0.0','_i5_R1.0_b0.0','_i12_R1.0_b0.0']
minTagsR1 = ['_i1_R1.0_b0.02','_i5_R1.0_b0.176','_i12_R1.0_b0.655']
minTags = ['_i1_R2.0_b0.045','_i5_R2.0_b0.35','_i12_R2.0_b1.31','_i1_R1.0_b0.02','_i5_R1.0_b0.176','_i12_R1.0_b0.655','_i1_R0.5_b0.013','_i5_R0.5_b0.088','_i12_R0.5_b0.327']
lotsTag = ['_i1_R1.0_b0.05','_i5_R1.0_b0.4','_i12_R1.0_b1.047']
sphereNdTags = ['_i12_R1.0_b0.0','_i12_R1.0_b0.655','_i12_R1.0_b1.35']
ico1NdTags = ['_i1_R1.0_b0.0','_i1_R1.0_b0.02','_i1_R1.0_b0.05']
ico5NdTags = ['_i5_R1.0_b0.0','_i5_R1.0_b0.176','_i5_R1.0_b0.4']
theseTags = zeroTags+minTagsR1+lotsTag

#theseTags = ['_i5_R1.0_b0.0','_i5_R1.0_b0.2','_i5_R1.0_b0.5']
#theseTags = ['_r0.01_i12_R1.0_b0.0','_i12_R1.0_b0.0','_r0.08_i12_R1.0_b0.0']
#theseTags = minEs_R1_r0_tags

plotGeos = []
for tag in theseTags:
    read_geo_fromTag(tag,resDir,plotGeos)

def plot(geos):    
    from mayavi import mlab
    PlotMayavi = imp.load_source('Geo_plotting_with_Mayavi','/data1/Garcia/OilDroplets/Code/PyScripts/Plot/Geo_plotting_with_Mayavi.py')
    
    mlab.clf()
#    #PlotMayavi.plot_geos_dist_to_d(singleGeo,0, 'type','wireframe')
#    #PlotMayavi.plot_geos_Density(singleGeo,'type','surface')
#    ##PlotMayavi.plot_geos_Density_sameScale(allGeos[13:18],'asphericity','surface')
#    #PlotMayavi.plot_geos_sigma(singleGeo,'core','surface')
     #PlotMayavi.plot_polarPair_Conf(geos,'type','surface',12)        #last argument 1->numDefects
#    #PlotMayavi.show_geos_sigma_sameScale(compareGeos,'type','surface')
#    #PlotMayavi.show_geos_sigma_sameScale(geos_vsS,'type','wireframe')
    #PlotMayavi.plot_geos_sigma(geos,'type','surface')
#    #PlotMayavi.plot_geos_sigma_withSeff(singleGeo,'type','surface')
#    #PlotMayavi.plot_geos_sigma_withDefectN(geos,'type','surface')
#    ##PlotMayavi.plot_geos_dist_to_d_sameScale(allGeos[5:12],2,'type','surface')
#    ##PlotMayavi.plot_geos_dist_to_d(allGeos[0:5],2,'type','surface')
#    #
#    ##PlotMayavi.plot_geos_W(allGeos[0:5],'type','surface')
    #PlotMayavi.plot_geos_sigma(geos,'type','surface')
   # PlotMayavi.plot_geos_sigma_sameScale(geos,'type','surface')
#    #PlotMayavi.plot_geos_Vorticity(singleGeo,'type','surface')
    #PlotMayavi.plot_geos_sigma_sameScale(geos,'type','surface')
    PlotMayavi.plot_geos_S(geos,'type','wireframe')
    #PlotMayavi.plot_geos_S_withSeff(geos,'type','surface')
    #PlotMayavi.plot_geos_S_sameScale_withSeff(geos,'type','surface')
    #PlotMayavi.plot_geos_S_sameScale(geos,'type','surface')
    #
    mlab.show()


#geos_vsS = [allGeos[0],allGeos[2],allGeos[1],allGeos[3]]
#singleGeo = [allGeos[1]]
#compareGeos = [find_geo('_r0.04_i12_b100.0_c0.04'),find_geo('_r0.04_i12_b0.0_c1.0')]

#compareGeos = [plotGeos[5],plotGeos[3],plotGeos[6]]
compareGeos = plotGeos

plot(compareGeos[:3])


##########################################################################
####   B Study Plots       ###############################################
########################################################################## 

