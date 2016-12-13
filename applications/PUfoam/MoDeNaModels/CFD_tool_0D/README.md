# README
@ingroup app_foaming
@brief This is a 0D example of the CFD code with the calls to the database using the MoDeNa interface.

## Description
This is a zero dimensional prototype of the application of MoDeNa interface library to connect different 
modeling scales. The first prototype shows the successful connections of nano, meso and macro sclaes for the simulation of foaming process.
The macro-scale code builds the RHS of 20 ODEs within the `QmomKinetics( const state_type &y , state_type &dydt , double t )` function including:
- dydt[0] : XW, conversion of the blowing reaction
- dydt[1] : XOH, conversion of the gelling reaction
- dydt[2] : T, temperature of the foam, K
- dydt[3] : L_l, weight fraction of the blowing agent in the liquid
- dydt[4] : L_g, weight fraction of the blowing agent in the gas
- dydt[5] : CO2_l, weight fraction of CO2 in the liquid
- dydt[6] : CO2_g, weight fraction of CO2 in the gas
- dydt[7] : m0, moment of order zero of the BSD
- dydt[8] : m1, moment of order one of the BSD
- dydt[9] : m2, moment of order two of the BSD
- dydt[10]: m3, moment of order three of the BSD
- dydt[11]: EG_NCO, Concentration of NCO end groups
- dydt[12]: EG_OH, Concentration of OH end groups
- dydt[13]: H2O, Water concentration
- dydt[14]: CO2, CO2 Concentration
- dydt[15]: PENTANE, Cylcopentane concentration
- dydt[16]: POLYMER, Dummy concentration of polymer
- dydt[17]: POLYMERBLOW, Second dummy concentration of polymer
- dydt[18]: UREA, Concentration of urea end groups
- dydt[19]: R_1_temp, Temperature

Throughout the computation of RHSs the inputs of the surrogate models are set and the outputs are retirieved using modena library.
It should be remembred that the surrogate models should be initialized on the local machine before running this code. 
The `write_kinetics(y, t)` function appends the results into the different text files that can be later plotted.
This function also computes the density of the foam based on the moments of BSD and converts the moments based 
on the unit volume of the foam.
The `main(int argc, char **argv)` function starts with the initialization of the solution. Then, a 'stepper' has been 
used to solve the ODEs. Further details on the integration method can be found [here.](http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/steppers.html)


## How to run independently?

1. Get Boost C++ Library
	* Skip this step if you have boost at /usr/local/ (I have used boost_1_57_0)
	* If not, get boost:
	\verbatim sudo apt-get install libboost-dev
	\endverbatim
2. Compile the project:
	* \verbatim cmake . \endverbatim
3. Make the executatble from the main code
	* \verbatim make \endverbatim
4. Set the environmental variable for the shared libraries by:
	* \verbatim export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./eigen \endverbatim
5. Run the executable by:
	* \verbatim ./QmomKinetics \endverbatim
6. Plot the results by software your choice.

### Note:
	In modifications of the solver are required,
	first run the cleanupScript.sh (./cleanupScript) to remove the old results
	then compile the modified code (make), 
	run the executable (./QmomKinetics) and plot the results (gnuplot gnuplot_script.gnu).

-------------

@authors Mohsen Karimi, Daniele Marchisio, Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
