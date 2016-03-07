

!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!> \file Main.F90
!! \brief Main program  
!!
!! THIS CODE WAS WRITTEN AT 
!! UNIVERSITY OF STUTTGART,
!! INSTITUTE OF TECHNICAL THERMODYNAMICS AND THERMAL PROCESS ENGINEERING
!! BY
!! JOACHIM GROSS
!! JONAS MAIRHOFER
!!
!!
!! This program calculates surface tensions using the a Density Functional 
!! Theory based on the PC-SAFTequation of state.
!! The contributions to the Helmholtz energy functional are calculated
!! as folows:
!! Hard Sphere Contribution: White-Bear Version of the Fundamental Measure Theory
!! Chain Formation: TPT1
!! Dispersion is treated in a  weighted density approximation
!! Associative and polar contributions are treated in a local density approximation
!!
!!
!! The input parameters are read from the file "in.txt" which has to
!! be in the same directory as the executable.
!!
!! The input file must have the following format:
!! Line1:       Value of temperature in Kelvin
!! Line2:       Number of components present in the system (ncomp)
!! Line3        Name of component 1
!! ...
!! Line3+ncomp Name of component ncomp
!! Line3+ncomp+1 Molar (overall) concentration of component 1
!! ...
!! Line3+2ncomp  Molar (overall) concentration of component ncomp
!!
!! For a binary system, these molar fractions are only treated as an initial guess and may be set to e.g. 0.
!!
!!
!! So far, pressure is set to 1bar in all calculaions
!!
!!
!!If you would like to use this code in your work, please cite the 
!!following publications:
!!
!!Gross, Joachim, and Gabriele Sadowski. "Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules." Industrial & engineering chemistry research 40.4 (2001): 1244-1260.
!!Gross, Joachim, and Gabriele Sadowski. "Application of the perturbed-chain SAFT equation of state to associating systems." Industrial & engineering chemistry research 41.22 (2002): 5510-5515.
!!Gross, Joachim, and Jadran Vrabec. "An equation-of-state contribution for polar components: Dipolar molecules." AIChE journal 52.3 (2006): 1194-1204.
!!Gross, Joachim. "A density functional theory for vapor-liquid interfaces using the PCP-SAFT equation of state." The Journal of chemical physics 131.20 (2009): 204705.
!!Klink, Christoph, and Joachim Gross. "A density functional theory for vapor-liquid interfaces of mixtures using the perturbed-chain polar statistical associating fluid theory equation of state." Industrial & Engineering Chemistry Research 53.14 (2014): 6169-6178.
!!
!!
!! In order to run this code, PETSc 3.4.4 has to be installed
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 


Program DFT

Use PARAMETERS, Only: np,nc,KBOL
Use BASIC_VARIABLES, Only: ncomp,t,p,ensemble_flag,num,parame,eos,pol,xif,compna
Use EOS_VARIABLES, Only: dhs,mseg,eta,dd_term,dq_term,qq_term
Use EOS_CONSTANTS, Only: ap,bp,dnm
Use mod_DFT, Only: box,dzp,fa,zp,fa_disp,ab_disp,pbulk
Use VLE_VAR, ONLY: tc,pc,rhob,density
USE DFT_FCN_MODULE, ONLY: chemPot_total

!PETSc modules
Use PetscManagement
Use f90moduleinterfaces
Use Global_x      !(snes,ngrid,ngp,x,r,timer)


Implicit None

#include <finclude/petscsys.h>

!> ---------------------------------------------------------------------- 
!!Variables 
!! ----------------------------------------------------------------------

PetscErrorCode     :: ierr
type (userctx)     :: user
PetscBool          :: flg
INTEGER            :: i
character(80)      :: filename=''
REAL               :: zges(np)
REAL               :: dhs_save(nc)
REAL               :: cif(nc)
REAL               :: psi_factor

!external subroutines associated with Solver
external FormInitialGuess
external FormFunction


!> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  Initialize program
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,user%rank,ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD,user%num_procs,ierr)

!> ---------------------------------------------------------------------- 
!!Read information from inputfile "in.txt"
!! ----------------------------------------------------------------------
 
      filename='./in.txt'
      CALL file_open(filename,77)       ! open input file
      READ (77,*) t                     ! read temperature
      READ (77,*) ncomp                 ! read number of components in the system
      Do i = 1,ncomp                    ! read component names
        READ (77,*) compna(i)
      End Do
      Do i = 1,ncomp                    ! read component overall molar concentrations
       READ (77,*) xif(i)
     End Do
 
!       !calculate molar fractions from molar concentrations
!       xif(1:ncomp) =cif(1:ncomp) / sum(cif(1:ncomp))
      
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GENERAL SIMULATION SET UP
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  num = 0                     ! (num=0: using analytical derivatives of EOS)
                              ! (num=1: using numerical derivatives of EOS)
                              ! (num=2: White's renormalization group theory)
  IF ( num /= 0 ) CALL SET_DEFAULT_EOS_NUMERICAL

  eos = 1
  pol = 1
  
  p = 1.000e05
  
  CALL para_input             ! retriev pure comp. parameters

  ensemble_flag = 'tp'        ! this specifies, whether the eos-subroutines 
                              ! are executed for constant 'tp' or for constant 'tv'

  
  
!> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  Start phase equilibrium calculation and determine critical point
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  CALL EOS_CONST (ap, bp, dnm)

      dd_term = 'GV'              ! dipole-dipole term of Gross & Vrabec (2006)
      qq_term = 'JG'              ! quadrupole-quadrupole term of Gross (2005)
      dq_term = 'VG'              ! dipole-quadrupole term of Vrabec & Gross (2008)

  CALL VLE_MIX(rhob,density,chemPot_total,user) !user is passed so user%rank is known and only node 0 prints out the reslts of the VLE calculation

  !determine critical point
    num = 0
    dhs_save = dhs               !needed because subroutine CRIT_POINT_MIX changes the value
    CALL CRIT_POINT_MIX(tc,user) !of the global variable dhs! In cases where Tc doesnt converge 
    dhs = dhs_save               !dhs can be NAN after CRIT_POINT_MIX!!
    !chance to overwite NAN results 
!     IF(user%rank == 0) THEN
!         WRITE (*,*) 'Give final value of Tc:'
!         READ (*,*) tc
!     END IF
!     CALL MPI_Bcast(tc,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr)



   
!>-------------------------------------------
!!DFT Set Up 
!!-------------------------------------------

   num = 1
   CALL SET_DEFAULT_EOS_NUMERICAL     

   !set default values (overwritten from makefile)
   box   = 45.0                !box length [A]
   ngrid = 614                 !grid points

   !overwrite box size and ngrid if options are set in makefile
   call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-box',box,flg,ierr)
   call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nx',ngrid,flg,ierr)


   dzp   = box / REAL(ngrid)   !grid spacing [A]
   Allocate(fa(ncomp))
   fa(1:ncomp) = NINT( parame(1:ncomp,2) / dzp  ) !grid points per sigma [-]
   
   !FOR DISP
   ALLOCATE(fa_disp(ncomp),ab_disp(ncomp))
   psi_factor = 1.5              
   fa_disp(1:ncomp) = NINT( psi_factor * REAL(fa(1:ncomp)) ) 
   Do i=1,ncomp
     if( MOD(fa_disp(i),2) /= 0 ) fa_disp(i) = fa_disp(i) + 1            
   End Do
 
   
   ab_disp(1:ncomp) = psi_factor * dhs(1:ncomp)      


   ngp = 2 * maxval(fa_disp(1:ncomp)) + 5   !number of ghost points (+5 kann eig weg)

   Allocate(zp(-ngp:ngrid+ngp))      

   Do i=-ngp,ngrid+ngp 
      zp(i) = REAL(i) * dzp          !z-distance from origin of grid points [A]
   End Do

   pbulk = ( p * 1.E-30 ) / ( KBOL*t* rhob(1,0) ) * rhob(1,0) !p/kT


   !update z3t, the T-dependent quantity that relates eta and rho, as eta = z3t*rho
   CALL PERTURBATION_PARAMETER
  
!>-------------------------------------------
!!Evaluate Initial Guess,Set up solver,solve system 
!!-------------------------------------------
     call SolverSetup



End Program DFT 
