#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Modified on Mon Aug 27 18:55:05 2018

@author: garcia-aguilar

write "def.tmp" file for c++ input, with previous info:
    - crystal coordination number "p"
    - Number of defects 
    - Defect info: (x,y,z) ; q (int) ; core_r
Include also the sign of the extra charge density induced by dislocation clouds around the disclinations.
    - qOff = +/- 1   [according to the disclination polar pair configuration]
    - qOff = 0       [to recover previous model with onlyDisclination]

NOW also include extra 20 defects in the centers of the faces of the icosahedron   

** Icosahedral geometry with one face against z = ztop = Rin

** Icosahedral geometry with rounded edge of eps = factor*edge on the x-coordinate

** Spherical geometry with golden ratio permutation of inscribed icosahedron
 

OilDroplets/Code/interface_sizeDisl_effS_File/Run/writeDefects_effS.py

"""

import numpy as np
#import config as cf


## In/Out directories
ouDir = "/data1/Garcia/OilDroplets/Code/interface_sizeDisl_effSFaces_File/Run/tst/Geos/"
#ouDir = "/home/garcia-aguilar/"
#ouDir = cf.inGeoDir

## Based on the raw mesh sphere R, write geo parameters
R0 = 1.0

## Crystal coordination
c = 6

## Disclination core radius  - Keep to fixed value according to BNT2000
core_ratio = 0.01          #[in R0 units]   # or list if for diff defects **Leave inside fxs
core_radius = core_ratio*R0

## Now qOff of the same sign for all vertices
 # List allows different weights on each vertex, and different weights for faces (eventually edges?)
qOff_vertices = [0., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.]

faceFraction = 12/20.       #Disclination/exta dislocation cores
qOff_faces = np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1])*faceFraction



triLists_onFace = [[3,2,1],[1,2,6],[1,6,7],[1,7,4],[1,4,3],
                   [3,5,2],[2,5,8],[2,8,6],[6,8,12],[6,12,7],[7,12,10],[7,10,4],[4,10,9],[4,9,3],[3,9,5],
                   [11,10,12],[11,12,8],[11,8,5],[11,5,9],[11,9,10]]

def write_line(fi,r, qDiscl, qOff, core_radius):
        fi.write(str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t\t"+str(qDiscl)+"\t\t"+str(qOff)+"\t\t"+str(core_radius)+"\n")
        
def find_faceCenter(r1,r2,r3):      #r1 np position array
    rMidEdge = (r2+r3)/2.
    rCenter = (r1+2*rMidEdge)/3.
    return rCenter

#Icosahedral symm  trans
def write_12Ico_centered(fname,R0, rcore,configDisclination, rcoreFace,configFaces):    
  #Defects
    print 'r_core',rcore    
    num_defects = 12+20         # 12 vertices and 20 faces    
    charge = 1.                  # or list if different charges    
    chargeFace = 0.               # no disclination there
      
  #Geometry
    factor_ico = np.power((16/5.)*(np.pi/(3+np.sqrt(5))),1/3.)      # cst Volume
    aa=factor_ico*R0
    print 'edge ',aa    

    ycenter = aa*np.cos(np.pi/6)/3.     # z axis along icosahedral center
    dihA = np.arccos(-np.sqrt(5)/3.)
    Rin = np.sqrt(3)*(3.0+np.sqrt(5))*aa/12
    print 'halfheight - Rin',Rin
    a3 = np.pi/6.
    a6 = np.pi/3.
   
  #File
    f = open(fname,"w")
    f.write("$ Crystal coordination $\n")
    f.write(str(c)+"\n")
    f.write("$ Number of defects $\n")
    f.write(str(num_defects)+"\n")
    f.write("$ Defect info $\n")
    f.write("x\t\ty\t\tz\t\tq\t\tqOff\t\tcore_r\n")
    
  #small top
    z1 = Rin                                        # geo centered in (0,0,0)    
    ytop = aa*np.cos(a3)
    r1 = np.array([aa/2,0.-ycenter,z1])
    write_line(f, r1, charge, configDisclination[1], rcore)      #V1
    r2 = np.array([-aa/2,0.-ycenter,z1])
    write_line(f, r2, charge, configDisclination[2], rcore)     #V2
    r3 = np.array([0,ytop-ycenter,z1])
    write_line(f, r3, charge, configDisclination[3], rcore)       #V3
    
  #big top       
    L = 2*aa*np.sin(dihA/2)*np.cos(a3)   
    z2 = -aa*np.sin(a6)*np.cos(dihA-np.pi/2.)+Rin   # geo centered in (0,0,0)   
    y_a = aa*np.cos(a3)*np.sin(dihA-np.pi/2.)
    y_b = np.cos(a3)*(L-aa*np.sin(dihA-np.pi/2.))    
    r4 = np.array([L/2,y_b-ycenter,z2])
    write_line(f, r4, charge, configDisclination[4], rcore)      #V4
    r5 = np.array([-L/2,y_b-ycenter,z2])
    write_line(f, r5,charge, configDisclination[5], rcore)      #V5
    r6 = np.array([0,-y_a-ycenter,z2])
    write_line(f, r6, charge, configDisclination[6], rcore)       #V6
    
  #big bottom    
    z3 = -z2                                        # geo centered in (0,0,0)    
    y_c = L/np.cos(a3)-aa*np.cos(a3)*np.sin(dihA-np.pi/2)   
    y_d = aa*np.cos(a3)*np.sin(dihA-np.pi/2.)+L*(np.sin(a3)-1)/(2.*np.cos(a3))
    r7 = np.array([L/2,-y_d-ycenter,z3])
    write_line(f, r7, charge, configDisclination[7], rcore)     #V7\
    r8 = np.array([-L/2,-y_d-ycenter,z3])
    write_line(f, r8, charge, configDisclination[8], rcore)     #V8
    r9 = np.array([0,y_c-ycenter,z3])
    write_line(f, r9, charge, configDisclination[9], rcore)        #V9
    
  #small bottom   
    z4 = -Rin                                       # geo centered in (0,0,0)   
    ybot = aa*(np.cos(a3)-np.tan(a3))  
    ybotT = aa*np.tan(a3)  
    r10 = np.array([aa/2,ybotT-ycenter,z4])
    write_line(f, r10, charge, configDisclination[10], rcore) #V10
    r11 = np.array([-aa/2,ybotT-ycenter,z4])
    write_line(f, r11, charge, configDisclination[11], rcore)#V11
    r12 = np.array([0,-ybot-ycenter,z4])
    write_line(f, r12, charge, configDisclination[12], rcore)    #V12    
    
    rs = [[0],r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12]
   
  #Faces
    for t, triangle in enumerate(triLists_onFace):
        rCenter = find_faceCenter(rs[triangle[0]],rs[triangle[1]],rs[triangle[2]])        
        write_line(f, rCenter, chargeFace, configFaces[t+1], rcoreFace)  
    
    f.close()

#Icosahedral symm  trans
def write_12Ico_roundingEps(fname,epsFactor,shiftInt,R0, rcore,configDisclination, rcoreFace,configFaces): 
  #Defects
    print 'r_core',rcore    
    num_defects = 12+20         # 12 vertices and 20 faces    
    charge = 1.                  # or list if different charges    
    chargeFace = 0.               # no disclination there  

      
  #Geometry
    factor_ico = np.power((16/5.)*(np.pi/(3+np.sqrt(5))),1/3.)
    aa = factor_ico*R0
    eps = epsFactor*aa
    print 'edge ',aa   
    
    #edFac = 1./3*np.sqrt((3+np.sqrt(5))/2)/np.sin(np.pi/5)/np.sin(2*np.pi/5)
    #aa = aa-eps*edFac/2
    
  #Rounding offsets 
    ycenter = aa*np.cos(np.pi/6)/3.         # z axis along icosahedral center
    #ySft = -eps/(np.sqrt(3)*4.)*shiftInt    # PREVIOUS  (?check)
    ySft_top = eps*2/(aa*np.sqrt(3))*shiftInt
    ySft_bot = eps/(np.sqrt(3))*shiftInt
    xSft = eps*shiftInt
    
   # print '\n\nEdge',aa,' yS_T',ySft_top,' yS_B',ySft_bot,' xS', xSft,'\n\n'    

    dihA = np.arccos(-np.sqrt(5)/3.)
    Rin = np.sqrt(3)*(3.0+np.sqrt(5))*aa/12
    print 'halfheight - Rin', Rin
    a3 = np.pi/6.
    a6 = np.pi/3.
  
  #File
    f = open(fname,"w")
    f.write("$ Crystal coordination $\n")
    f.write(str(c)+"\n")
    f.write("$ Number of defects $\n")
    f.write(str(num_defects)+"\n")
    f.write("$ Defect info $\n")
    f.write("x\t\ty\t\tz\t\tq\t\tcore_r\n")
    
  #small top
    #z1 = 0.0                           # in z<0
    z1 = Rin                            # geo centered in (0,0,0)    
    ytop = aa*np.cos(a3) - ySft_top
    r1 = np.array([aa/2-xSft,0.-ycenter+ySft_bot,z1])
    write_line(f, r1, charge, configDisclination[1], rcore)      #V1
    r2 = np.array([-aa/2+xSft,0.-ycenter+ySft_bot,z1])
    write_line(f, r2, charge, configDisclination[2], rcore)     #V2
    r3 = np.array([0,ytop-ycenter,z1])
    write_line(f,r3, charge, configDisclination[3], rcore)                     #V3
    
  #big top      
    L = 2*aa*np.sin(dihA/2)*np.cos(a3)  
    #z2 = -aa*np.sin(a6)*np.cos(dihA-np.pi/2.)         # in z<0    
    z2 = -aa*np.sin(a6)*np.cos(dihA-np.pi/2.)+Rin      # geo centered in (0,0,0)#   
    y_a = aa*np.cos(a3)*np.sin(dihA-np.pi/2.) + ySft_top   
    y_b = np.cos(a3)*(L-aa*np.sin(dihA-np.pi/2.)) -ySft_bot    #***CHK
    r4 = np.array([L/2-xSft,y_b-ycenter,z2])
    write_line(f,r4, charge, configDisclination[4], rcore)             #V4
    r5 = np.array([-L/2+xSft,y_b-ycenter,z2])
    write_line(f,r5, charge, configDisclination[5], rcore)            #V5
    r6 = np.array([0,-y_a-ycenter,z2])
    write_line(f, r6, charge, configDisclination[6], rcore)                   #V6
    
  #big bottom
    #z3 = -2*Rin-z2                     # in z<0
    z3 = -z2                            # geo centered in (0,0,0)    
    y_c = L/np.cos(a3)-aa*np.cos(a3)*np.sin(dihA-np.pi/2) - ySft_top   
    y_d = aa*np.cos(a3)*np.sin(dihA-np.pi/2.)+L*(np.sin(a3)-1)/(2.*np.cos(a3)) + ySft_bot   ##****CHK
    r7 = np.array([L/2-xSft,-y_d-ycenter,z3])
    write_line(f, r7, charge, configDisclination[7], rcore)            #V7
    r8 = np.array([-L/2+xSft,-y_d-ycenter,z3])
    write_line(f, r8, charge, configDisclination[8], rcore)           #V8
    r9 = np.array([0,y_c-ycenter,z3])
    write_line(f, r9, charge, configDisclination[9], rcore)                    #V9
    
  #small bottom
    #z4 = -2*Rin                        # in z<0
    z4 = -Rin                           # geo centered in (0,0,0)  
    ybot = aa*(np.cos(a3)-np.tan(a3)) - ySft_top
    ybotT = aa*np.tan(a3)  - ySft_bot
    r10 = np.array([aa/2-xSft,ybotT-ycenter,z4])
    write_line(f, r10, charge, configDisclination[10], rcore)         #V10
    r11 = np.array([-aa/2+xSft,ybotT-ycenter,z4])
    write_line(f, r11, charge, configDisclination[11], rcore)        #V11
    r12 = np.array([0,-ybot-ycenter,z4])
    write_line(f, r12, charge, configDisclination[12], rcore)                 #V12        
    
    rs = [[0],r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12]
   
  #Faces
    for t, triangle in enumerate(triLists_onFace):
        rCenter = find_faceCenter(rs[triangle[0]],rs[triangle[1]],rs[triangle[2]])        
        write_line(f, rCenter, chargeFace, configFaces[t+1], rcoreFace)  
        
    f.close()      
      
   


############################################################################
####### Set defects at ~ rounding circle 
#############################

##*## Ico old Rounding
#*# set:
ouFile = "tst_defRound_faces.tmp"
inCalc = ouDir
#*# write:
##write_12Ico_roundingEps(inCalc+ouFile, 0.05, 1., R0, 0.04, qOff_vertices, 0.04, qOff_faces)
    

##*## Ico Centered (evolver)
#*# set:
ouFile = "def_faces.tmp"
inCalc = ouDir
#*# write:
write_12Ico_centered(inCalc+ouFile, R0, 0.04, qOff_vertices, 0.04, qOff_faces)

