# CurvedCrystals

## General outline

Numerical integration of the Poisson equation for the surface stress on a curved crystal with topological dislocations and disclinations (see Eq.(7) in [García-Aguilar,Fonda, Giomi. PRE 2020](https://doi.org/10.1103/PhysRevE.101.063005)). In the continuum model for curved crystals developed, these defects are expressed effectively as point "charges" on the surface. For the numerical integration, the geometry of the crystal is fixed and given by a triangulated mesh. The defects are considered independently of the mesh, as a collection of points in space, wich an associated core size (note this is for the numerical approximation of the effective charge density) and charge. 

The numerical integrator to find the stress field on the frozen surface, given the pre-defined defects in written in C++. It requires a few input files and results in various output files; these are described below. Python scripts were used to run the integrator for various parameter, and later to read and process the output files. 

Example files are provided for the sharp icosahedron in the PRE, tagged with "ico5". 

### In general, the steps are:
1. Create initial frozen geometry, as a discretized triangulated mesh (.msh).
2. Create defect file.
3. (OPTIONAL) Calculate the distances from each vertex in the mesh to the pre-determined defect position. Save a new geometry file that includes those distances (.mshd)
4. Run numerical integration of stretching field
5. Calculate energy values, or read the geometry+fields
   
The frozen geometry to be used is stored in .msh or .mshd file (they have different formats but both can be used in the numerical integrator).
Files .msh contain the information of vertices positions and the triplets forming the faces of the discretization, in a format similar to the output from Evolver (**see for example ico5.msh**). Files .mshd contain the mesh together with the distances between each vertex and each defect pre-defined on the mesh (**see for example tri_ico5_3.mshd**). 

The calculation of this distance can be done within the numerical integration too, when only a .msh file is provided. Do take note on which calculation method to use (set in the configuration file) because more accurate methods are also considerably slower, so it becomes impractical. For crystal geometries where the Euclidean distance is a good approximation to the geodesic distance for the vertices closest to the defect, slower methods can be used, and it is not necessary to genereate a .mshd file. 

## Generate initial geometry 
The example geometry (ico5.msh) was created using Surface Evolver. The starting configuration was a perfect icosahedron, which was refined to the desired mesh resolution. The geometry was left to relax under volume constrain, in order to round the sharp edges. Other regularization steps were done in order to avoid obstuse angles in the triangulation. 

## Generate defect file
The defects on the surface are pre-defined through an input file (see **defs_faces.tmp**). The format is simple to read from opening the file. However, to link to the Evolver geometries, the python scripts **writeDefects_\*.py** can be used, to approximate the defect positions close to the input frozen geometries. This is relevant because when calculating the distance between each vertex in the mesh to the defect, the center of the defect is approximated to the nearest mesh vertex. The new positions are stored in one of the output files for the numerical integration and/or the calculation of the geodesic distances, (see for example **defs_ico5_R1.0_b0.06.dat**), which can be used for further calculations on the surface.

For the numerical integration, the effect of each "point" defect is approximated by a Gaussian function centered in the defect core, with standard deviation equal to the core_size given in the run configuration file.  

## (OPTIONAL) Calculation of geodesic distances : .msh -> .mshd file
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
- File with defect information (see for example **"defs_ico5_3.dat"**), including a column with the vertex index closest to the input defect position, and therefore used for the numerical geodesic distance.
  
## Run numerical integration

#### Input files

#### Output files

## Get results:
Either DataFrame with energy calculations OR read out fields on python on a Shape object

#### Input / output files:



#### Scripts:


