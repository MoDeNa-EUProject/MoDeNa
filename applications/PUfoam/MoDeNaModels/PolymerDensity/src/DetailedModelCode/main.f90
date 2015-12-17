
!>
!!THIS CODE WAS WRITTEN AT 
!!UNIVERSITY OF STUTTGART,
!!INSTITUTE OF TECHNICAL THERMODYNAMICS AND THERMAL PROCESS ENGINEERING
!!BY
!!JOACHIM GROSS AND JONAS MAIRHOFER
!!


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This program calculates densities of a liquid phase which is in equilibrium with a gas phase
!! using the PC-SAFT equation of state.
!! The input parameters are read from the file "in.txt" which has to
!! be in the same directory as the executable.
!!
!! The input file must have the following format:
!! Line1:       Value of temperature in Kelvin
!! Line2:       Number of components present in the system (ncomp)
!! Line3        Name of component 1
!! ...
!! Line3+ncomp Name of component ncomp
!! Line3+ncomp+1 Molar (overall) fraction of component 1
!! ...
!! Line3+2ncomp  Molar (overall) fraction of component ncomp
!!
!! For a binary system, these molar fractions are only treated as an initial guess and may be set to e.g. 0.
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
!!
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 

PROGRAM PC_SAFT
!
! ----------------------------------------------------------------------
 USE BASIC_VARIABLES
 USE EOS_CONSTANTS
 USE DFT_MODULE
 USE EOS_VARIABLES!, ONLY: dd_term, qq_term, dq_term
 IMPLICIT NONE

!> ---------------------------------------------------------------------- 
!!Variables 
!! ----------------------------------------------------------------------

 REAL                                  :: tc,pc,chemPot_total(nc)
 REAL                                  :: rhob(2,0:nc),density(np)
 CHARACTER (LEN=50)                    :: filename
 INTEGER                               :: compID,i
 REAL                                  :: cif(nc)  

!> ---------------------------------------------------------------------- 
!!Read information from inputfile "in.txt"
!! ----------------------------------------------------------------------

  filename='./in.txt'
  CALL file_open(filename,77)
  READ (77,*) t
  READ (77,*) ncomp                 ! read number of components in the system
  Do i = 1,ncomp                    ! read component names
      READ (77,*) compna(i)
  End Do
  Do i = 1,ncomp                    ! read component overall molar concentrations
       READ (77,*) xif(i)
  End Do
 
 
 
  !calculate molar fractions from molar concentrations
  !xif(1:ncomp) = cif(1:ncomp) / sum(cif(1:ncomp))


 
!> ---------------------------------------------------------------------- 
!!General simulation set up
!! ----------------------------------------------------------------------

  num = 1                     ! (num=0: using analytical derivatives of EOS)
                              ! (num=1: using numerical derivatives of EOS)
                              ! (num=2: White's renormalization group theory)
  IF ( num /= 0 ) CALL SET_DEFAULT_EOS_NUMERICAL
 

  eos = 1
  pol = 1
 
  p = 1.000e05 
 
  CALL para_input            ! retriev pure comp. parameters
 
  ensemble_flag = 'tp'        ! this specifies, whether the eos-subroutines 
                              ! are executed for constant 'tp' or for constant 'tv'
 

!> ----------------------------------------------------------------------
!!Start phase equilibrium calculation
!! ----------------------------------------------------------------------

      CALL EOS_CONST (ap, bp, dnm)

      dd_term = 'GV'              ! dipole-dipole term of Gross & Vrabec (2006)
      qq_term = 'JG'              ! quadrupole-quadrupole term of Gross (2005)
      dq_term = 'VG'              ! dipole-quadrupole term of Vrabec & Gross (2008)
 


      CALL VLE_MIX(rhob,density,chemPot_total,compID)

         
END PROGRAM PC_SAFT 


