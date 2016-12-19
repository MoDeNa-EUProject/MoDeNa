@ingroup mod_elasticmodulus

# Foam Elastic Modulus

Having microstructure containing bubbles with specific size distribution, elastic response of polyurethane foams can be affected by bubble size distribution, density of foam and strut_content. Bubble size distribution could be described by two values: mean value, MU and standard deviation, SIGMA (Assuming bubble size follows the normal distribution).
In addition, rho and strut_content parameters stand for density and strut content, respectively.
As a result, ''Foam Elastic Modulus'' has been coded  in a way that it receives mean value and standard deviation of bubble size distribution as well as rho and strut_content in foam microstructure and returns the elastic modulus of polyurethane foam. 
The folder ''FoamElasticModulus'' contains ''AbaqusSimulation.py''.
By default, the directory which has the ''AbaqusSimulation.py'', is used as a working directory. Following steps made the ''FoamElasticModulus.py'' ready to be executed:

1. The folder includes the surce code of ''spherepack'' that is called by ''AbaqusSimulation.py'' during execution. Consecuently, in advance to use the ''AbaqusSimulation.py'', ''spherepack'' tool must be compiled. To compile ''spherepack'', you need to run the cmake in Terminal as:

```bash
cmake .
```

And then:

```bash
make
```

By compiling the “spherpack”, new file of “spherepack.exe” is constructed.
Enough care must be taken that to have the “spherepack.exe” in the same directory as ''AbaqusSimulation.py''. 


2. ''AbaqusSimulation.py .py''  calls the ABAQUS in order to run the input file generated. To make the ABAQUS able to use the ''AbaqusSimulation.py'' directory as its own directory (normally ABAQUS uses “commands” folder in home directory as a working directory), it is recommended to use folloeing commands in order to make the ABAQUS executable in any directory on your PC: 

```bash
sudo ln -s /usr/local/abaqus/Commands/abq<yourversion> /usr/bin/abq<yourversion>
```

3. It  is also necessary to introduce the version of ABAQUS. To do that, lines 42 of the ''AbaqusSimulation.py'' must be modified based on user's  ABAQUS version.  Also, line 43 is asking for number of cores available on computing machine (By default is 20).




