#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Modified on Mon Jan 27 17:49:05 2020

@author: garcia-aguilar

-------- for HEXAGONAL shape ----------------

write "def.tmp" file for c++ input, with previous info:
    - crystal coordination number "p"
    - Number of defects 
    - Defect info: (x,y,z) ; q (int) ; core_r
Include also the sign of the extra charge density induced by dislocation clouds around the disclinations.
    - qOff_sink =  -1   [to screen positive disclination charge]
    - qOff = 1*numSources/numSinks       

HERE: include only 2 sinks at the hexagonal face centers   

** Hexagonal geometry with hexagonal face against +dz*a/2n
 

OilDroplets/Code/PyScripts/writeDefects_hex.py

"""

import numpy as np

## In/Out info
ouDir = "/data1/Garcia/OilDroplets/Geometries/Evolver/Hexagonal/"

## Based on the raw mesh sphere R, write geo parameters
R0 = 1.0

## Evolver geometry 
h_hex = 1/10.  #in side a-units
factor_hex = 8.*np.pi/(9*np.sqrt(3)*h_hex)
side = np.power(factor_hex,1/3.)*R0
coordx = np.sqrt(3)*side/2
coordy = side/2
height = h_hex*side

## Crystal coordination
c = 6

## Disclination core radius  - Keep to fixed value 
core_ratio = 0.04          #[in R0 units]   # or list if for diff defects **Leave inside fxs
rcore = core_ratio*R0

##
 # List allows different weights on each vertex, and different weights for dislocation sinks
qOff_vertices = [0., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.]

faceFraction = 12/20.       #Disclination/exta dislocation cores
qOff_faces = np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1])*faceFraction



def write_line(fi,r, qDiscl, qOff, core_radius):
	fi.write(str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t\t"+str(qDiscl)+"\t\t"+str(qOff)+"\t\t"+str(core_radius)+"\n")
 
def add_hexagon(f, r, q, qOff, rcore):
    vertex = np.array([1,1,1])
    for j in [-1,1]:
        for k in [-1,1]:
		    vertex = np.array([j, k, 1.])
		    write_line(f, r*vertex, q, qOff, rcore)           #[+-a Sin60, +-a Cos60, z]   
        vertex = np.array([0.,j*side,r[2]])
        write_line(f, vertex, q, qOff, rcore)				  #[0, +-a, z]
  
def write_Hex(fname, numDiscl=6):

  ## Defect configuration   --- assume symmetry of disclinations and eventual dislocation sinks	
    charge = 12./numDiscl 
    num_defects = numDiscl+2		# assume 2 sinks                
  
  # Open file and add heading
    f = open(fname,"w")
    f.write("$ Crystal coordination $\n")
    f.write(str(c)+"\n")
    f.write("$ Number of defects $\n")
    f.write(str(num_defects)+"\n")
    f.write("$ Defect info $\n")
    f.write("x\t\ty\t\tz\t\tq\t\tqOff\t\tcore_ratio\n")

  # Disclinations
    qOff_source = -1.          #Assume source symmetry in all disclination cores    
    coordz = 0.
    coord = [coordx,coordy,coordz]       
    if numDiscl != 6:         # If 12 disclinations, instead of 6, add the additional ones
		coord[2] = height/2.
		add_hexagon(f, coord*np.array([1.,1.,-1.]), charge, qOff_source, rcore)
    add_hexagon(f, coord, charge, qOff_source, rcore)
        
  # Dislocations sinks (add two at each hexagonal face center)
    qOff_sink = 1.*numDiscl/2.     		#Disclination/exta dislocation sinks -- assume sink symmetry
    write_line(f, [0.,0., height/2], 0., qOff_sink, rcore)
    write_line(f, [0.,0.,-height/2], 0., qOff_sink, rcore)
 
    f.close()


############################################################################
####### Set defects 
#############################

##*## Hexagonal prism or 2Pyramid with q=2 disclinations
#*# set:
ouFile = "defs_hex_faces.tmp"
inCalc = ouDir
#*# write:
write_Hex(inCalc+ouFile)


##*## Hexagonal prism or 2Pyramid with q=1 disclinations
#*# set:
ouFile = "def_hex12_faces.tmp"
inCalc = ouDir
#*# write:
write_Hex(inCalc+ouFile, 12)
    


