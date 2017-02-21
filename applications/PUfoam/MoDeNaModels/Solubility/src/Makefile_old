 FC=     gfortran 
#  FC=     ifort

# ---- compiling options for ifort: debuging ----
# FFLAGS= -r8 -O3 -132 -warn all -check bounds -traceback 
# FFLAGS= -r8 -O3 -132 -warn all -check bounds 
# FFLAGS= -r8 -O3 -132 -check bounds 
# ---- compiling options for ifort ----
# FFLAGS= -r8 -O3
# FFLAGS= -r8 -O3 -132 -opt-report 

# ---- compiling options for gfortran ----
# FFLAGS= -O3 -fdefault-real-8 -Wall -Wextra -std=f2003 
# the -std=f90 is removed to allow the definition of allocatable variables in type definitions
FFLAGS= -O3 -fdefault-real-8 -Wall -Wextra
# ---- use for debugging only ----
# FFLAGS= -O3 -fdefault-real-8 -Wall -Wextra -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow # -fbacktrace

LFLAGS= $(FFLAGS)

EXE=    pcsaft 

SRC= mod_basic_variables.f90  hbrd.f90  levenberg_marquardt.f90  2D_optimizer.f90  cg_descent.f90  utilities.f90  eos_enthalpy.f90  eos_const.f90  minimize.f90  starting_value.f90  DFT-utilities.f90  eos_polar.f90  DFT-nMF2.f90  DFT-nMF-units.f90  DFT-nMF-mixtures.f90  DFT-MF.f90  eos_f_contributions.f90 eos_numerical.f90  eos_critical_renorm.f90  eos_critical_renorm_mix.f90  eos_polar_second_deriv.f90  eos.f90  input_output.f90  para_input.f90  pure_par_fit.f90  kij-fitting.f90 BubblePoint.f90  main_prog.f90  phase_equilib.f90  gc_method.f90  force_field_fit.f90  rdf_hs.f90  pt1.f90  visco-fit.f90  crit_point_mixtures.f90  matrix_inversion2.f90 

#
## No need to edit below this line
#
#
#.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core *.mod *genmod.f90


