""" Python configuration file for parameters """
""" Specific for run file with three nested loops for three varying parameters. At the moment implemented for changing 
droplet size, dislocation flux (dislStrenght or  burgerNumber is various parts of the code), and dislocation extend (which 
was not used in the PRE study)."""



# Mesh Geometry
inGeoDir = '../Geo/'
ico = 5
mshName = 'ico'+str(ico)+'.msh'
mshName = 'tri_ico'+str(ico)+'_3.mshd'     # "3" here means the distance calculation was done previously using Exact algorithm
distance_flag = 4 		#(1. Circle, 2. Dijkstra, 3. Exact, 4. FromFile, 5. Euclidean)
mshFile = inGeoDir + mshName
R0_mesh = 1.0		#prel value


# Droplet Geometry
range_R_init = [1.0]		#loop over droplet sizes *** keep > core_sizeRatio
latticeConstant = 5E-5			# [10mu]from work in oil/water interfaces ***look


# Defects
facesNumber = 20
checkConf = ''			# '' for main conf
defName = 'def_faces.tmp'
defFile = inGeoDir + defName
core_sizeRatio = 0.04           		#single value  core_radius=core_sizeRatio*R_init =0.04
dislStrength = 1.0   
range_dislStrength = [0.08]       # in [1/lattic_a] units in c++
dislExtend = 0.0	
range_dislExtend = [0.0]


# Integration parameters
#num_it = 500000
num_it = 500        # TST
#run_time = 100.		// change to DT***
time_step = 5E-3
integration_flag = 2  #(1. Euler, 2. RK-2steps, 3. RK-4steps)


# Whole Sim
N_runs = 1  # init, while no param study
N_runs = len(range_R_init)*len(range_dislExtend)*len(range_dislStrength)#* anyOtherParamList


# Init
init_stress = 0.     
range_stress = 0.2
seed = 17
convergence_error_stop = 1E-23




# Output Results
preTag = '_i'+str(ico)+''			# init value
outDir = '../Results/'   ## 


# Program to use
progDir = '../interface_sizeDisl_effSFaces_File/bin/Release/'
progVersion = 'interface_effSFaces'
prog = progDir+progVersion
