!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This module contains constants and upper boundaries for the number of
!! components and phases in the system
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module PARAMETERS

  implicit none
  save

  integer, parameter                            :: nc = 6
  integer, parameter                            :: np = 3
  integer, parameter                            :: nsite = 5

  real, parameter                               :: PI   = 3.141592653589793
  real, parameter                               :: RGAS = 8.31441
  real, parameter                               :: NAv  = 6.022045E23
  real, parameter                               :: KBOL = RGAS / NAv
  real, parameter                               :: TAU  = PI / 3.0 / SQRT(2.0)

End Module PARAMETERS


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This module contains parameters and variables that define the 
!! thermodynamic state of the system as well as simulation parameters
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module BASIC_VARIABLES

  use PARAMETERS, only: nc, np, nsite
  implicit none
  save

! ----------------------------------------------------------------------
! basic quantities defining the mixture
! ----------------------------------------------------------------------
  integer                                       :: ncomp
  integer                                       :: nphas

  real                                          :: t
  !real						:: tc
  real                                          :: p
  real, dimension(np)                           :: dense
  !real, dimension(np)				:: rhob

  real, dimension(np, nc)                       :: xi
  real, dimension(np, nc)                       :: lnx
  real, dimension(nc)                           :: xiF

  real, dimension(nc)                           :: mm
  real, dimension(np, nc, nsite)                :: mxx

  real                                          :: alpha


! ----------------------------------------------------------------------
! pure component parameters
! ----------------------------------------------------------------------
  real, dimension(nc, 25)                       :: parame = 0.0
  real, dimension(nc)                           :: chiR
  character*30, dimension(nc)                   :: compna
  real, dimension(nc, nc)                       :: kij, lij
  real, dimension(nc, nc)                       :: E_LC, S_LC
  real, dimension(nc)                           :: LLi, phi_criti, chap


! ----------------------------------------------------------------------
! thermodynamic quantities
! ----------------------------------------------------------------------
  real, dimension(np)                           :: densta
  real, dimension(0:nc*np+6)                    :: val_init, val_conv

  real                                          :: h_lv
  real, dimension(np)                           :: cpres, enthal, entrop, gibbs, f_res

  real, dimension(np)                           :: dp_dz, dp_dz2


! ----------------------------------------------------------------------
! choice of EOS-model and solution method
! ----------------------------------------------------------------------
  integer                                       :: eos, pol
  integer                                       :: num
  character (LEN=2)                             :: ensemble_flag
  character (LEN=10)                            :: RGT_variant


! ----------------------------------------------------------------------
! for input/output
! ----------------------------------------------------------------------
  integer                                       :: outp, bindiag
  real                                          :: u_in_T, u_out_T, u_in_P, u_out_P


! ----------------------------------------------------------------------
! quantities defining the numerical procedure
! ----------------------------------------------------------------------
  integer                                       :: n_unkw

  real                                          :: step_a, acc_a !, acc_i
  real, dimension(nc)                           :: scaling
  real, dimension(3500)                         :: plv_kon
  real, dimension(2, 3500)                      :: d_kond

  character*3, dimension(10)                    :: it, sum_rel
  character*3                                   :: running


End Module BASIC_VARIABLES




!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This module contains parameters and variables to identify thermodynamic
!! properties
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module EOS_VARIABLES

  use PARAMETERS, only: nc, nsite, PI, KBOL, TAU, NAv
  use BASIC_VARIABLES, only: ncomp, eos, t, p, parame, E_LC, S_LC, chiR, &
                             LLi, phi_criti, chap, kij, lij, ensemble_flag
  implicit none
  save

! ----------------------------------------------------------------------
! basic quantities defining the mixture (single phase)
! ----------------------------------------------------------------------
  real                                     :: x(nc)
  real                                     :: eta_start
  real                                     :: eta
  real                                     :: rho

! ----------------------------------------------------------------------
! thermodynamic quantities
! ----------------------------------------------------------------------
  real                                     :: fres
  real                                     :: lnphi(nc)
  real                                     :: pges
  real                                     :: pgesdz
  real                                     :: pgesd2
  real                                     :: pgesd3

  real                                     :: h_res
  real                                     :: cp_res
  real                                     :: s_res

! ----------------------------------------------------------------------
! quantities of fluid theory
! ----------------------------------------------------------------------
  real                                     :: gij(nc,nc)
  real                                     :: mx(nc,nsite)

  real                                     :: mseg(nc)
  real                                     :: dhs(nc)
  real                                     :: uij(nc,nc)
  real                                     :: sig_ij(nc,nc)
  real                                     :: vij(nc,nc)

  real                                     :: um
  real                                     :: order1
  real                                     :: order2
  real                                     :: apar(0:6)
  real                                     :: bpar(0:6)

  real                                     :: z0t
  real                                     :: z1t
  real                                     :: z2t
  real                                     :: z3t

  integer                                  :: nhb_typ(nc)
  real                                     :: ass_d(nc,nc,nsite,nsite)
  real                                     :: nhb_no(nc,nsite)
  real                                     :: dij_ab(nc,nc)

! ----------------------------------------------------------------------
! auxilliary quantities
! ----------------------------------------------------------------------
  real                                     :: tfr
  integer                                  :: phas

  character (LEN = 2)                      :: dd_term, qq_term, dq_term

  real                                     :: densav(3), denold(3)
  real                                     :: density_error(3)

  real                                     :: alpha_nematic
  real					   :: alpha_test(2)

End Module EOS_VARIABLES



!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This module contains parameters and variables of the PC-SAFT EoS
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module EOS_CONSTANTS

  use PARAMETERS, only: nc
  implicit none
  save

  real, dimension(0:6,3)                   :: ap, bp
  real, dimension(4,9)                     :: dnm

  real, dimension(28)                      :: c_dd, n_dd, m_dd, k_dd, o_dd
  real, dimension(nc,nc,0:8)               :: qqp2, qqp4, ddp2, ddp4, dqp2, dqp4
  real, dimension(nc,nc,nc,0:8)            :: qqp3, ddp3, dqp3

End Module EOS_CONSTANTS 

!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This module contains parameters and variables that control the 
!! numerical approximation of derivaties of thermodynamic quantities
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module EOS_NUMERICAL_DERIVATIVES

  use EOS_VARIABLES, only: dd_term, qq_term, dq_term

  implicit none
  save

  character (LEN=9)                        :: ideal_gas      ! (yes, no)
  character (LEN=9)                        :: hard_sphere    ! (CSBM, no)
  character (LEN=9)                        :: chain_term     ! (TPT1, no)
  character (LEN=9)                        :: disp_term      ! (PC-SAFT, CK, no)
  character (LEN=9)                        :: hb_term        ! (TPT1_Chap, no)
  character (LEN=9)                        :: LC_term        ! (MSaupe, no)
  character (LEN=9)                        :: branch_term    ! (TPT2, no)
  character (LEN=9)                        :: II_term        ! (......., no)
  character (LEN=9)                        :: ID_term        ! (......., no)

  character (LEN=9)                        :: subtract1      ! (1PT, 2PT, no)
  character (LEN=9)                        :: subtract2      ! (ITTpolar, no)

  character (LEN=9)                        :: save_eos_terms(10)

End Module EOS_NUMERICAL_DERIVATIVES


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! Module STARTING_VALUES
!! This module contains parameters and variables for a phase stability
!! analyis as part of a flash calculation.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

 Module STARTING_VALUES

  use PARAMETERS, only: nc
  implicit none
  save

  REAL, DIMENSION(nc)                   :: rhoif, rhoi1, rhoi2, grad_fd
  REAL                                  :: fdenf

 End Module STARTING_VALUES





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module DFT_MODULE
!
! This module contains parameters and variables for DFT calculations.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 Module DFT_MODULE

  use PARAMETERS, only: nc
  implicit none
  save

  INTEGER, PARAMETER                    :: NDFT = 4000
 !!integer		:: discret
 REAL 			:: box_l_no_unit
  INTEGER, PARAMETER                    :: r_grid = 200
  INTEGER                               :: kmax, den_step
  LOGICAL                               :: shift, WCA, MFT
  REAL                                  :: rc, rg, dzr, tau_cut,dzp
  REAL                                  :: d_hs, dhs_st, z3t_st
  REAL                                  :: z_ges
  REAL, DIMENSION(r_grid)               :: x1a
  REAL, DIMENSION(NDFT)                 :: x2a
  REAL, DIMENSION(r_grid,NDFT)          :: ya, y1a, y2a, y12a
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub
  REAL                                  :: fres_temp

  REAL, DIMENSION(r_grid)               :: x1a_11
  REAL, DIMENSION(NDFT)                 :: x2a_11
  REAL, DIMENSION(r_grid,NDFT)          :: ya_11, y1a_11, y2a_11, y12a_11
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub_11
  REAL, DIMENSION(r_grid)               :: x1a_12
  REAL, DIMENSION(NDFT)                 :: x2a_12
  REAL, DIMENSION(r_grid,NDFT)          :: ya_12, y1a_12, y2a_12, y12a_12
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub_12
  REAL, DIMENSION(r_grid)               :: x1a_22
  REAL, DIMENSION(NDFT)                 :: x2a_22
  REAL, DIMENSION(r_grid,NDFT)          :: ya_22, y1a_22, y2a_22, y12a_22
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub_22

 End Module DFT_MODULE


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module ........... 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
Module rdf_variables

  implicit none
  save

  real, dimension(0:20)                         :: fac(0:20)

End Module rdf_variables


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains the variables that are needed in the core DFT_FCN
! They are not passed directly to DFT_FCN because the nonlinear solver 
! needs a certain calling structure: fcn(x,fvec,n) 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
Module DFT_FCN_MODULE

use PARAMETERS, only: nc
use DFT_MODULE, only: ndft
  implicit none

INTEGER					:: irc(nc),irc_j,ih,fa(nc)
 REAL, DIMENSION(-NDFT:NDFT)            :: zp
 REAL, DIMENSION(-NDFT:NDFT)            :: f_tot
 REAL, DIMENSION(-NDFT:NDFT,2)          :: dF_drho_tot
 REAL					:: rhob_dft(2,0:nc)
 REAL                                   :: my0(nc)

End Module DFT_FCN_MODULE








!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module Module_Heidemann_Khalil
!
! This module ....
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 Module Module_Heidemann_Khalil

  implicit none
  save

  real                                  :: error_condition2

 End Module




