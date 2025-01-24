# CurvedCrystals

## General outline

Numerical integration on a frozen geometry, given by a triangulated surface. 

The crystal geometry to be used is stored in .msh or .mshd file (they have very different formats but both can be used in the numerical integrator).
Files .msh contain the information of vertices positions and the triplets forming the faces of the discretization, in a format similar to the output from Evolver (**see for example ico5.msh**). Files .mshd contain the mesh together with the distances between each vertex and each defect pre-defined on the mesh (**see for example tri_ico5.mshd**). 

The calculation of this distance can be done within the numerical integration too, when only a .msh file is provided. Do take not on which calculation method to define in the configuration file because more accurate methods are also considerably slower, so it becomes impractical. For crystal geometries where the Euclidean distance is a good approximation to the geodesic distance for the vertices closest to the defect, slower methods can be used, and it is not necessary to genereate a .mshd file. 

### In general, the steps are:
1. Create initial crystal geometry
2. Create defect file
3. Run numerical integration of stretching field
4. Calculate energy values, or read the geometry+fields

## Generate initial geometry 

## Calculation of geodesic distances : .msh -> .mshd file
(Not necessary but makes calculations faster and more accurate)
#### Input files

#### Output files

## Generate defect file

## Run numerical integration

#### Input files

#### Output files

## Get results:
Either DataFrame with energy calculations OR read out fields on python on a Shape object

#### Input / output files:



#### Scripts:


