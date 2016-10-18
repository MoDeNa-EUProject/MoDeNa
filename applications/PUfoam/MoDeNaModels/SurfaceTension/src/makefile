 
#include /usr/ITT/mhofer/Documents/Diss/NumericalMethods/Libraries/Petsc/petsc-3.4.4/conf/variables
#include /usr/ITT/mhofer/Documents/Diss/NumericalMethods/Libraries/Petsc/petsc-3.4.4/conf/rules

#include ~/NumLib/PETSc/petsc-3.4.4/conf/variables
#include ~/NumLib/PETSc/petsc-3.4.4/conf/rules


include /${PETSC_DIR}/conf/variables
include /${PETSC_DIR}/conf/rules



SOURCE = Modules.F90 \
	mod_PETSc.F90 \
	mod_DFT_FMT.F90 \
	mod_DFT_FMT_d.F90 \
	mod_DFT_CHAIN.F90 \
	mod_DFT_CHAIN_d.F90 \
	mod_DFT_DISP_WDA.F90 \
	mod_DFT_DISP_WDA_d.F90 \
	module_solve_nonlinear.F90 \
	getting_started_subroutines.F90 \
	Helfer_Routinen.F90 \
	Numeric_subroutines.F90 \
	Spline_Integration_d.F90 \
	VLE_main.F90 \
	VLE_subroutines.F90 \
	crit_point_mixtures.F90 \
	Function.F90 \
	AD_Routines.F90 \
	InitialGuess.F90 \
	SolverSetup.F90 \
	Main.F90 \

#Object files
OBJECT = $(SOURCE:%.F90=%.o)

#define target for non-PETSc files
%.o: %.F90
	${PETSC_FCOMPILE} -fdefault-real-8 -c $< -o $@

DFT: $(OBJECT)
	-${FLINKER} -o PCSAFT_SurfaceTension $(OBJECT) ${PETSC_SNES_LIB} 

#--------------------------------------------------------------------------

#Anzahl Prozessoren
NP = 1

#Initial profile: 0: normal, 1-3: add perturbation to regular initial profile 
PERT = 0

#DFT Settings
NGRID     = 800
CUTOFF    = 9.0
BOXSIZE   = 300.0 

#Solver Settings
ITS_SNES = 15
ITS_KSP  = 10
E_REL    = 1e-08

#Its for Anderson Mixing
ITS_SNES_ANDERSON = 100
#Its for Picard Iterations
ITS_SNES_PICARD = 100

#Toleranzen
ATOL_SNES = 1e-08
RTOL_SNES = 1e-08
STOL_SNES = 1e-08

ATOL_KSP = 1e-08
RTOL_KSP = 1e-04

#Dämpfungsfaktoren
DAMP = 0.3
DAMP_LBFGS = 1.0
DAMP_ANDERSON = 0.01
DAMP_PICARD = 0.01

#------------------------------------------------------------------------------
#1) Inexact-Newton type solvers (iterative solver for linear system)
#------------------------------------------------------------------------------

#matrix-free, numerical approximation of directional derivatives (choose between trust region and line search)
run_mf:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-erel $(E_REL) \
	-init_pert $(PERT) \
	-snes_type newtonls  -snes_converged_reason  \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-ksp_gmres_cgs_refinement_type refine_always \
	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
	-jac 0 \
	-pc_type none \



#matrix-free, AD-calculation directional derivatives (choose between trust region and line search)
run_ad:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type newtonls  -snes_converged_reason  \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-ksp_gmres_cgs_refinement_type refine_always \
	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
	-jac 1 \
#	-pc_type none

#finite-difference approximation of Jacobi matrix (-> PC can be used!) (choose between trust region and line search)
run_fd:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension -snes_fd -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type newtonls  -snes_converged_reason  \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-ksp_gmres_cgs_refinement_type refine_always \
	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
	-jac 2 \
	-dm_mat_type aij \
#	-pc_type ilu

#Build complete Jacobi matrix with AD (-> PC can be used!) (choose between trust region and line search)
run_anad:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type newtonls  -snes_converged_reason  \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-ksp_gmres_cgs_refinement_type refine_always \
	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
	-jac 3 \



#------------------------------------------------------------------------------
#2) Newton type solvers (direct solver for linear system)
#------------------------------------------------------------------------------

#finite-difference approximation of Jacobi matrix (-> PC can be used!) (choose between trust region and line search)
run_fdmumps:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type newtonls  -snes_converged_reason  \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-jac 2 \
	-dm_mat_type aij \
	-pc_type lu  -pc_factor_mat_solver_package mumps \

#Build complete Jacobi matrix with AD (-> PC can be used!) (choose between trust region and line search)
run_anadmumps:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type newtonls  -snes_converged_reason  \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-jac 3 \
	-pc_type lu  -pc_factor_mat_solver_package mumps \


#------------------------------------------------------------------------------
#3) Quasi-Newton type solvers (secant updates for approximation to inverse Jacobian)
#------------------------------------------------------------------------------

#limitd memory quasi newton with BFGS updates -> choose restart parameter
run_LBFGS:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type qn -snes_qn_type lbfgs \
	-snes_converged_reason \
	-snes_linesearch_type cp -snes_linesearch_damping $(DAMP_LBFGS) \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-jac -1 \

#limitd memory quasi newton with BFGS updates -> choose restart parameter
#calculate initial jacobian with AD, then use BFGS updates in later iterations

#anpassen damping, jac 0/1

run_LBFGSinitJac:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type qn -snes_qn_type lbfgs -snes_qn_scale_type jacobian \
	-snes_converged_reason \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP_LBFGS) \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-jac 0 \


#use a direct solver for linear system

#limitd memory quasi newton with BFGS updates -> choose restart parameter
#calculate initial jacobian with AD, then use BFGS updates in later iterations
run_LBFGSinitJac2:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type qn -snes_qn_type lbfgs -snes_qn_scale_type jacobian \
	-snes_converged_reason \
	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP_LBFGS) \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
	-jac 3 \



#------------------------------------------------------------------------------
#4) Simple Fixpoint Iterations
#------------------------------------------------------------------------------


#Anderson mixing
run_Anderson:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type anderson \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES_ANDERSON) \
	-snes_converged_reason \
	-snes_linesearch_type l2 -snes_linesearch_damping $(DAMP_ANDERSON) -snes_linesearch_monitor \
	-jac -1 \


#Picard Iteration
run_Picard:
	-@${MPIEXEC} -n $(NP) ./PCSAFT_SurfaceTension   -snes_monitor_short -ksp_monitor_short \
	-nx $(NGRID) \
	-rc $(CUTOFF) \
	-box $(BOXSIZE) \
	-init_pert $(PERT) \
	-snes_type nrichardson \
	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES_PICARD) \
	-snes_converged_reason \
	-snes_linesearch_type l2 -snes_linesearch_damping $(DAMP_PICARD) -snes_linesearch_monitor \
	-jac -1 \




# # #--------------------------------------------------------------------------
# # 
# # #Anzahl Prozessoren
# # NP = 1
# # 
# # 
# # #DFT Settings
# # NGRID     = 800
# # CUTOFF    = 9.0
# # BOXSIZE   = 10.0 
# # 
# # #Solver Settings
# # ITS_SNES = 40
# # ITS_KSP  = 20
# # E_REL    = 1e-11
# # 
# # #Its for Anderson Mixing
# # ITS_SNES_ANDERSON = 250
# # #Its for Picard Iterations
# # ITS_SNES_PICARD = 250
# # 
# # #Toleranzen
# # ATOL_SNES = 1e-08
# # RTOL_SNES = 1e-08
# # STOL_SNES = 1e-08
# # 
# # ATOL_KSP = 1e-08
# # RTOL_KSP = 1e-04
# # 
# # #Dämpfungsfaktoren
# # DAMP = 0.5
# # DAMP_LBFGS = 0.5
# # DAMP_ANDERSON = 0.5
# # DAMP_PICARD = 0.01
# # 
# # #------------------------------------------------------------------------------
# # #1) Inexact-Newton type solvers (iterative solver for linear system)
# # #------------------------------------------------------------------------------
# # 
# # #matrix-free, numerical approximation of directional derivatives (choose between trust region and line search)
# # runNLsolver_mf:
# # 	-@${MPIEXEC} -n $(NP) ./NLsolver -snes_monitor_short -ksp_monitor_short \
# # 	-nx $(NGRID) \
# # 	-rc $(CUTOFF) \
# # 	-box $(BOXSIZE) \
# # 	-erel $(E_REL) \
# # 	-snes_type newtonls  -snes_converged_reason  \
# # 	-snes_atol $(ATOL_SNES) -snes_rtol $(RTOL_SNES) -snes_stol $(STOL_SNES) -snes_max_it $(ITS_SNES) \
# # 	-ksp_atol  $(ATOL_KSP) -ksp_rtol $(RTOL_KSP) -ksp_max_it $(ITS_KSP) -ksp_gmres_restart 50 \
# # 	-snes_linesearch_type bt -snes_linesearch_damping $(DAMP) -snes_linesearch_monitor \
# # 	-snes_max_fail 1 -snes_max_linear_solve_fail 100 \
# # 	-ksp_gmres_cgs_refinement_type refine_always \
# # 	-snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
# # 	-jac 0 \
# # 	-pc_type none \
