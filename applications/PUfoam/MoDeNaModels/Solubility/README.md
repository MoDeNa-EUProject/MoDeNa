 
# Solubility

## Scope of this module
The solubility module calculates an effective Henry coefficient of a gaseous component in a mixture.
The module consists of two main components: 

* The detailed model
* The surrogate model

The [PCSAFT](http://pubs.acs.org/doi/abs/10.1021/ie0003887) equation of state is used as the detailed model.
The detailed model takes temperature and the liquid composition as input values and subsequently performs a
bubble point calculation in order to obtain the composition of the coexisting vapor phase and returns an effective 
Henry coefficient as result. This coefficient is defined as $H_i = \frac{y_i p}{x_i} $ where $x_i$ and $y_i$ denote
the molar fractions of component i in the liquid and vapor phase, respectively, and p denotes the equilibrium pressure.
Which components are present in the system is defined using index sets.

The surrogate model is a simple exponential function with three adjustable parameters, A, B and C:
$H_i(T) = A \cdot \text{exp} (B(1/T - 1/C)) $. 

A nonlinear least-squares algorithm is used to determine the optimal values of these parameters in order to correlate the 
results of the detailed model as closely as possible. 




## Prerequisites
In order to compile and run the module a fortran compiler, preferably gfortran, needs to be installed. A makefile to compile the
detailed model code is provided. 

## Compilation and Execution of the detailed model

* Compilation: make 
* Execution: ./pcsaft

## Input / Output

Module inputs:

* Temperature (Kelvin)
* Composition of liquid phase (Molar fractions)

Module outputs:

* Effective Henry coefficient (bar)


## Workflow

On initialisation of the module, the surrogate model parameters are adjusted to the initial data points which have to be set in the module prior to 
program execution. Subsequently, as long as the input temperatures do not fall outside this initial temperature intervall, only the surrogate model
is evaluated at a model call. Once the module is called with a temperature outside this intervall, new data points are generated according to
the *out of bounds strategy* which also needs to be specified in the module and the parameters of the surrogate model are readjusted. 




