""" Python run file-> change parameters at config.py 
  
 Parameter scan on initial spherical droplet radii Dislocation density strength (dislStrength or burgerNumber)  and dislocation extend (or at most only three param) """


import config as cf
import os
import subprocess
import sys
import shutil



# Simulation parameters
NIt = cf.num_it
#Rtime = cf.run_time
tStep = cf.time_step
int_flag = cf.integration_flag
Nruns = cf.N_runs


# Init
init_stress = cf.init_stress
range_stress = cf.range_stress
seed = cf.seed
conv_stop = cf.convergence_error_stop


# Geometry
mshFile = cf.mshFile
mshName =  cf.mshName
mshRadius = cf.R0_mesh
sphereRadii = cf.range_R_init
lattice_a = cf.latticeConstant


distanceFlag = cf.distance_flag


# Defects
defectFile = cf.defFile
defectName = cf.defName
core_ratio = cf.core_sizeRatio	#Keep core size/R constant for all defects for now  
burgerNumbers = cf.range_dislStrength		# in [1/lattic_a] units
c_extends = cf.range_dislExtend
NFaces = cf.facesNumber




# Check files exist or print error
def check_if_file_exists(someFile):
	if os.path.isfile(someFile):		
		pass
	else:
		print "*** -> NOT found : "+someFile
		print "  -- Check file is in directory or check config.py "
check_if_file_exists(defectFile)
check_if_file_exists(mshFile)




#writeDefects.write()


# Input mesh / defs
inDir = cf.inGeoDir
### TODO: cp temporarily mshFile and defectFile in ouDir


# Output Results
pretag = cf.preTag     # initial val
out = cf.outDir 


# Program to use
Prog = cf.prog


# Start
print "Using :",Prog
print "\n+++ START runs for: ",cf.mshFile
print "    Running ",Nruns," calculations"


## first loop
runCounter = 0
for radius in sphereRadii:
	print "\n ***** Start droplet SIZE "
	print "*********\n\t R =", radius
	### second loop
	for burgerNumber in burgerNumbers:
		print "\n********* Start Dislocation B-number"
		print "*************\n\t\t b =", burgerNumber
		#### third loop
		for c_extend in c_extends:	
			tag = pretag+'_R'+str(radius)+'_b'+str(burgerNumber)
        #tag = pretag+'_R'+str(radius)+'_b'+str(burgerNumber)+'_c'+str(c_extend)
			
			print "\n\t\t\t Tag = ", tag


			# Write def.tmp file where to read defects approx positions
			#written = writeDefects.write()
			#if written is False:
			#	exit()


			# Write input file for ./interface
			def write_input():
				fname = "inputParameters.txt"
				f = open(fname,"w")
				f.write("MeshFile         -          {0}\n".format(mshFile))
				f.write("MeshName         -          {0}\n".format(mshName))	
				f.write("DefectFile       -          {0}\n".format(defectFile))
				f.write("DistanceMethod   -          {0}\n".format(distanceFlag))
				f.write("(1. Great-circle, 2. Dijkstra, 3. Exact, 4. ReadFromFile, 5. Euclidean)\n")
				f.write("DropletRadius	.	MeshRadius	-	{0}	.	{1}\n".format(radius,mshRadius))
				f.write("LatticeConstant	-		{0}\n".format(lattice_a))
				f.write("CoreRatio       -          {0}\n".format(core_ratio))
				f.write("DislocBurger   .  DislocExtend  -   {0}  .  {1}\n".format(burgerNumber,c_extend))
				f.write("Init_stress  .  Range   -   {0}  .  {1}\n".format(init_stress,range_stress))
				f.write("Seed	      -          {0}\n".format(seed))
				#f.write("N_steps      .  T_run   -   {0}  .  {1}\n".format(NIt,Rtime))	// change to DT**
				f.write("N_steps      .  time_step   -   {0}  .  {1:.2E}\n".format(NIt,tStep))
				f.write("Integration      -          {0}\n".format(int_flag))
				f.write("(1. Euler, 2. RK-2steps, 3. RK-4steps)\n")
				f.write("ConvergenceStop  _	        {:.3E}\n".format(conv_stop))
				f.close()  


			write_input()


			print "\n-> ____________________________________________\n"
			def execute_CppProgram(program):
				with open(out+'screen_cat_out'+tag+'.dat', "w") as oFile:
					oFile.write("Using :  "+program)
			# run the simulation
				CProcess = subprocess.Popen(program, bufsize=1, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
				stdoutProc = []
				textOut = ''
				while True:
					line = CProcess.stdout.readline()
					stdoutProc.append(line)
					textOut = textOut+line
					print line,
					if line == '' and CProcess.poll() != None:
						break
						
				(COut, CErr) = CProcess.communicate()
				
				if CProcess.returncode != 0:
					print "C error :\n\t"+ str(CErr)		
				
				with open(out+'screen_cat_out'+tag+'.dat', "w") as oFile:
					oFile.write(textOut)




			execute_CppProgram(Prog)






			print "\t\t\t C++ done"
			print "____________________________________________ <-\n\n"




			## Loop over file_name_list and move those that exist. 	******* TODO
			def move_C_outFiles():	
				out_file_names = ['ge','hi','last','defs','ou']


				for OutName in out_file_names:
					currentFile = OutName+'.dat'	
					if os.path.isfile(currentFile):			# os.path.exists for dir
						targetFile = out+OutName+tag+'.dat'
						if os.path.isfile(targetFile):
							print '\tThere is already a '+OutName+tag+'.dat  file'	
							targetFile = out+OutName+'-old'+tag+'.dat'
							print '  \t--> renaming this to '+OutName+'-old'+tag+'.dat'				
						os.rename(currentFile, targetFile)
					else:
						print "\t*** -> NO "+OutName+'.dat  file'
						
				 # and s-t*dat files:
				currentFiles = os.listdir('.')
				for anyFile in currentFiles:	
					if (anyFile.startswith("s-") and anyFile.endswith(".dat")):   #at least one s-t* file
						targetDir = out+'TimeData'+tag
						if os.path.exists(targetDir):
							print '\tThere is already a TimeData'+tag, 'directory'			
							print '  \t--> renaming this to TimeData-old'+tag
							os.rename(targetDir, out+'TimeData-old'+tag)		
							os.makedirs(targetDir)			
						else:						
							os.makedirs(targetDir)
							
						for someFile in currentFiles:
							if (someFile.startswith("s-") and someFile.endswith(".dat")):
								os.rename(someFile,targetDir+'/'+someFile)		
						break


			move_C_outFiles()


			print "\n********* END param set =..  "+ tag
			runCounter = runCounter+1
			print " \t- - - - Still "+str(Nruns-runCounter)+" runs  of "+str(Nruns)+"- - - -"
		  #End third loop


		print "\n********* END Burger = "+ str(burgerNumber)	
		## End second loop
		
	print "\n******* END Radius = "+ str(radius)	
	## End first loop


print "\n+++ END runs for: ",cf.mshName,"\n+++++++++++++++++++\n"                


