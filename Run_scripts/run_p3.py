""" Python run file-> change parameters at config.py 
  
 Parameter scan on initial spherical droplet radii Dislocation density strength (dislStrength or burgerNumber)  and dislocation extend (or at most only three param) 

# Python 3 version of the original Python 2 code
"""
import config as cf
import os
import subprocess
import sys
import shutil

# Simulation parameters
NIt = cf.num_it
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
mshName = cf.mshName
mshRadius = cf.R0_mesh
sphereRadii = cf.range_R_init
lattice_a = cf.latticeConstant
distanceFlag = cf.distance_flag

# Defects
defectFile = cf.defFile
defectName = cf.defName
core_ratio = cf.core_sizeRatio
burgerNumbers = cf.range_dislStrength
c_extends = cf.range_dislExtend
NFaces = cf.facesNumber

# Check files exist or print error
def check_if_file_exists(someFile):
    if not os.path.isfile(someFile):
        print("*** -> NOT found : " + someFile)
        print("  -- Check file is in directory or check config.py ")

check_if_file_exists(defectFile)
check_if_file_exists(mshFile)

# Input mesh / defs
inDir = cf.inGeoDir

# Output Results
pretag = cf.preTag
out = cf.outDir

# Program to use
Prog = cf.prog

# Start
print("Using:", Prog)
print("\n+++ START runs for:", cf.mshFile)
print("    Running", Nruns, "calculations")

# First loop
runCounter = 0
for radius in sphereRadii:
    print("\n ***** Start droplet SIZE ")
    print("*********\n\t R =", radius)
    for burgerNumber in burgerNumbers:
        print("\n********* Start Dislocation B-number")
        print("*************\n\t\t b =", burgerNumber)
        for c_extend in c_extends:
            tag = f"{pretag}_R{radius}_b{burgerNumber}"
            print("\n\t\t\t Tag =", tag)

            # Write input file for ./interface
            def write_input():
                fname = "inputParameters.txt"
                with open(fname, "w") as f:
                    f.write(f"MeshFile         -          {mshFile}\n")
                    f.write(f"MeshName         -          {mshName}\n")
                    f.write(f"DefectFile       -          {defectFile}\n")
                    f.write(f"DistanceMethod   -          {distanceFlag}\n")
                    f.write("(1. Great-circle, 2. Dijkstra, 3. Exact, 4. ReadFromFile, 5. Euclidean)\n")
                    f.write(f"DropletRadius\t.\tMeshRadius\t-\t{radius}\t.\t{mshRadius}\n")
                    f.write(f"LatticeConstant\t-\t\t{lattice_a}\n")
                    f.write(f"CoreRatio       -          {core_ratio}\n")
                    f.write(f"DislocBurger   .  DislocExtend  -   {burgerNumber}  .  {c_extend}\n")
                    f.write(f"Init_stress  .  Range   -   {init_stress}  .  {range_stress}\n")
                    f.write(f"Seed\t      -          {seed}\n")
                    f.write(f"N_steps      .  time_step   -   {NIt}  .  {tStep:.2E}\n")
                    f.write(f"Integration      -          {int_flag}\n")
                    f.write("(1. Euler, 2. RK-2steps, 3. RK-4steps)\n")
                    f.write(f"ConvergenceStop  _\t        {conv_stop:.3E}\n")

            write_input()

            print("\n-> ____________________________________________\n")

            def execute_CppProgram(program):
                with open(out + "screen_cat_out" + tag + ".dat", "w") as oFile:
                    oFile.write("Using :  " + program)
                CProcess = subprocess.Popen(program, bufsize=1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                stdoutProc = []
                textOut = ""
                for line in CProcess.stdout:
                    stdoutProc.append(line)
                    textOut += line
                    print(line, end="")

                (COut, CErr) = CProcess.communicate()

                if CProcess.returncode != 0:
                    print("C error :\n\t" + str(CErr))

                with open(out + "screen_cat_out" + tag + ".dat", "w") as oFile:
                    oFile.write(textOut)

            execute_CppProgram(Prog)

            print("\t\t\t C++ done")
            print("____________________________________________ <-\n\n")

            def move_C_outFiles():
                out_file_names = ['ge', 'hi', 'last', 'defs', 'ou']
                for OutName in out_file_names:
                    currentFile = OutName + '.dat'
                    if os.path.isfile(currentFile):
                        targetFile = out + OutName + tag + '.dat'
                        if os.path.isfile(targetFile):
                            print(f'\tThere is already a {OutName}{tag}.dat file')
                            targetFile = out + OutName + '-old' + tag + '.dat'
                            print(f'  \t--> renaming this to {OutName}-old{tag}.dat')
                        os.rename(currentFile, targetFile)
                    else:
                        print(f"\t*** -> NO {OutName}.dat file")

                currentFiles = os.listdir('.')
                for anyFile in currentFiles:
                    if anyFile.startswith("s-") and anyFile.endswith(".dat"):
                        targetDir = out + 'TimeData' + tag
                        if os.path.exists(targetDir):
                            print(f'\tThere is already a TimeData{tag} directory')
                            print(f'  \t--> renaming this to TimeData-old{tag}')
                            os.rename(targetDir, out + 'TimeData-old' + tag)
                        os.makedirs(targetDir, exist_ok=True)

                        for someFile in currentFiles:
                            if someFile.startswith("s-") and someFile.endswith(".dat"):
                                os.rename(someFile, targetDir + '/' + someFile)
                        break

            move_C_outFiles()

            print("\n********* END param set =..  " + tag)
            runCounter += 1
            print(f" \t- - - - Still {Nruns - runCounter} runs  of {Nruns}- - - -")

        print("\n********* END Burger =", burgerNumber)
    print("\n******* END Radius =", radius)

print("\n+++ END runs for:", cf.mshName, "\n+++++++++++++++++++")
