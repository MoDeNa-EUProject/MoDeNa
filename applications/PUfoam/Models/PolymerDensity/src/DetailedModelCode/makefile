 

#Source files
SOURCE = modules.f90 \
	module_solve_nonlinear.f90 \
	getting_started_subroutines.f90 \
	Numeric_subroutines.f90 \
	VLE_subroutines.f90 \
	VLE_main.f90\
	main.f90 \




#Object files
OBJECT = $(SOURCE:%.f90=%.o)



#define target for non-PETSc files
%.o: %.f90
	gfortran -c -fdefault-real-8 $< -o $@

EOS: $(OBJECT)
	gfortran -o PCSAFT_Density $(OBJECT) 

clean:
	rm *.o *.mod



