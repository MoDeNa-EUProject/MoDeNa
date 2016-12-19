@ingroup mod_surfacetension

Surface Tension
===============

## Scope of this module

The surface tension module returns the surface tension of a fluid system in equilibrium with coexisting 
liquid and vapor phases. The module consists of two main components: 

* The detailed model
* The surrogate model

A PC-SAFT based density functional theory is used as the detailed model, the details of which are outlined in \cite gross_2009 .
However, here we use a modification which is computationally less demanding and can be applied to systems of more than two components.
The detailed model takes temperature as input.
Pressure is set to one bar in all calculations.
The output of the model is the surface tension of the system.

Since surface tension shows an almost perfectly linear temperature dependence, the surrogate model is
a simple linear equation in temperature: $\\gamma(T) = A + BT$ where A and B are adjustable parameters.

A nonlinear least-squares algorithm is used to determine the optimal values of these parameters in order to correlate the results of the detailed model as closely as possible. 

[publication]: http://scitation.aip.org/content/aip/journal/jcp/131/20/10.1063/1.3263124

## Prerequisites {#surfacetensionPrerequisite}

In order to compile and run the module a fortran compiler, preferably gfortran, as well as Portable, Extensible Toolkit for Scientific Computation ([PETSc]) 3.4.5 need to be installed.
[PETSc] 3.4.5 requires specific versions of MPI, BLAS, Lapack and Scalapack.
In order to ensure compatibility, PETSc should be configured to automatically download and install the correct versions, see section Installing PETSc.

As there is no backward compatibility of different versions of PETSc, older as well as newer version of PETSc will most likely not work.

[PETSc]: (https://www.mcs.anl.gov/petsc/)


## Installing PETSc {#surfacetensionPETSc}

PETSc should be configured with the following options:

```bash
./configure --with-cc=gcc --with-fc=gfortran --download-f-blas-lapack --download-mpich --download-scalapak=yes
```

Furthermore, the variables PETSC_DIR and PETSC_ARCH need to be set.

## Compiling and running the detailed model code {#surfacetensionCompileAndRun}

### Compile {#surfacetensionCompile}

```
cd src
make DFT
```

### Run {#surfacetensionRun}

```
cd src
make run_mf
```

## Input / Output {#surfacetensionIO}


| Input       | Symbol | Unit |
| ----------- | ------ | ---- |
| Temperature | T      | K    |


| Output          | Symbol | Unit |
| --------------- | ------ | ---- |
| Surface Tension | ST     | mN/m |


## Workflow {#surfacetensionWorkflow}

On initialisation of the module, the surrogate model parameters are adjusted to the initial data points which have to be set in the module prior to program execution.
Subsequently, as long as the input temperatures do not fall outside this initial temperature intervall, only the surrogate model is evaluated at a model call.
Once the module is called with a temperature outside this intervall, new data points are generated according to the *out of bounds strategy* which also needs to be specified in the module and the parameters of the surrogate model are readjusted. 




