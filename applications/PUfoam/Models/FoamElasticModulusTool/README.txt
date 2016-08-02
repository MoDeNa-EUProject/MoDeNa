Foam Elastic Modulus:
Having microstructure containing bubbles with specific size distribution, elastic response of polyurethane foams can be affected by bubble size distribution. Bubble size distribution could be described by two values: mean value, MU and standard deviation, SIGMA (Assuming bubble size follows the normal distribution). As a result, ''Foam Elastic Modulus'' has been coded  in a way that it receives mean value and standard deviation of bubble size distribution in foam microstructure and returns the elastic modulus of polyurethane foam. 
The folder ''Foam Elastic Modulus'' contains ''FoamElasticModulus.py''. By default, the directory which has the ''Foam Elastic Modulus.py'', is used as a working directory. Following steps made the ''FoamElasticModulus.py'' ready to be executed:

1- The folder includes the surce code of ''SpherePackFB'' that is called by ''Foam Elastic Modulus.py'' during execution. Consecuently, in advance to use the ''FoamElasticModulus.py'', ''SpherePackFB'' must be compiled. To compile ''SpherePackFB'', it is recommended to use the Dev-C++ compiler. “SpherPackFB.dev” could be imported by Dev-C++ and compiled. By compiling the “SpherPackFB.dev”, new file of “SpherePackFB.exe” is constructed. The ''FoamElasticModulus.py'' calls the “SpherePackFB.exe” using “wine”. “wine” is a application that makes it possible to execute the executable files in Windows (e.g. SpherePackFB.exe) in Linux as well. The “wine” could be installed via UBUNTU software center. Enough care must be taken that to have the “SpherePackFB.exe” in the same directory as ''Foam Elastic Modulus.py''. 

2- ''FoamElasticModulus.py''  calls the ABAQUS in order to run the input file generated. To make the ABAQUS able to use the ''FoamElasticModulus.py'' directory as its own directory (normally ABAQUS uses “commands” folder in home directory as a working directory), it is recommended to use folloeing commands in order to make the ABAQUS executable in any directory on your PC: 
sudo ln -s /usr/local/abaqus/Commands/abq<yourversion> /usr/bin/abq<yourversion>

3- It  is also necessary to introduce the version of ABAQUS. To do that, line 16 of the ''FoamElasticModulus.py'' must be modified based on user's  ABAQUS version. 

4- Dependencies: The following softwares (or modules) have to be installed before running the ''FoamElasticModulus.py'':

I. Python 2.7.10.
II. Python modules:
. os — Miscellaneous operating system interfaces
. os.path — Common pathname manipulations
. itertools — Functions creating iterators for efficient looping
. NumPy  fundamental package for scientific computing with Python.
. math — Mathematical functions
. SciPy Stack
. time — Time access and conversions
. commands — Utilities for running commands
. glob — Unix style pathname pattern expansion
. Random (Generate pseudo-random numbers)
III. IDLE (Integrated DeveLopment Environment or Integrated Development and Learning Environment)  for Python (using python 2.7.10). Since  incompatibility between some python IDE's (e.g. Pycharm) and ABAQUS has been seen, it is highly recommended to use IDLE to run ''Foam Elastic Modulus''.
IV. NEPER (e.g. NEPER 2.0.3). NEPER could be downloaded through: ''http://neper.sourceforge.net/ ''
V. Gmsh 2.10.1  (http://gmsh.info/)
VI. ABAQUS FEA software package (e.g. ABAQUS 6.13).  
