# CurvedCrystals

## General outline

Numerical integration on a frozen geometry, given by a triangulated surface. 

The crystal geometry to be used is stored in .msh or .mshd file (they have different formats but both can be used in the numerical integrator).
Files .msh contain the information of vertices positions and the triplets forming the faces of the discretization, in a format similar to the output from Evolver (**see for example ico5.msh**). Files .mshd contain the mesh together with the distances between each vertex and each defect pre-defined on the mesh (**see for example tri_ico5.mshd**). 

The calculation of this distance can be done within the numerical integration too, when only a .msh file is provided. Do take notd on which calculation method to define in the configuration file because more accurate methods are also considerably slower, so it becomes impractical. For crystal geometries where the Euclidean distance is a good approximation to the geodesic distance for the vertices closest to the defect, slower methods can be used, and it is not necessary to genereate a .mshd file. 

### In general, the steps are:
1. Create initial crystal geometry
2. Create defect file
3. Run numerical integration of stretching field
4. Calculate energy values, or read the geometry+fields

## Generate initial geometry 
The example geometry (ico5.msh) was created using Surface Evolver. The starting configuration was a perfect icosahedron, which was refined to the desired mesh resolution. The geometry was left to relax under volume constrain, in order to round the sharp edges. Other regularization steps were done in order to avoid obstuse angles in the triangulation. 

## Calculation of geodesic distances : .msh -> .mshd file
(Not necessary but makes more accurate calculations run faster)

The code to calculate distances between vertices on a discretized mesh and predefined defects at positions (x,y,z) is implemented in C++, and found in Code/Calculate_distances/
Several algorithms are posible: 
1. arc length on a circle (for the spherical geometries)
2. Geodesic distante using Dijkstra algorithm
3. Exact geodesic distance 
5. Euclidean distance

See https://code.google.com/archive/p/geodesic/ for more information on geodesic distances 2 and 3. All these are implemented in the numerical integrator and can be used here. However, the most accurate 3. Exact algorithm is computationally intense. Since these distances need to be calculated only once for a particular geometry (including defect geometry), it is best to perform this calculation apart to generate the **.mshd** file, and use that one as an input to the numerical integrator. 

#### Input files
- Input meshed geometry (.msh file)
- Defect positions (def.tmp file). This will also be an input file for the numerical integration

#### Output files
- New geometry file with the same name but different suffix (.mshd), which includes the calculated distances
- File with defect information (see for example "defs_ico5_3.dat"), including a column with the vertex index closest to the input defect position, and therefore used for the numerical geodesic distance. This file is only relevant for checking the calculation method or the input distances.
  
## Generate defect file

## Run numerical integration

#### Input files

#### Output files

## Get results:
Either DataFrame with energy calculations OR read out fields on python on a Shape object

#### Input / output files:



#### Scripts:


