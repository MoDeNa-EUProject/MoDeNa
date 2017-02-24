!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains parameters and variables ....
! ................ 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module PARAMETERS

  implicit none
  save

  integer, parameter                            :: nc = 10
  integer, parameter                            :: np = 3
  integer, parameter                            :: nsite = 5

  real, parameter                               :: PI   = 3.141592653589793
  real, parameter                               :: RGAS = 8.31441
  real, parameter                               :: NAv  = 6.022045E23
  real, parameter                               :: KBOL = RGAS / NAv
  real, parameter                               :: TAU  = PI / 3.0 / SQRT(2.0)

End Module PARAMETERS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains parameters and variables ....
! ................ 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

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
  real                                          :: p
  real, dimension(np)                           :: dense

  real, dimension(np, nc)                       :: xi
  real, dimension(np, nc)                       :: lnx
  real, dimension(nc)                           :: xiF
  ! real, dimension(nc)                           :: rhoi_f
  real, dimension(nc)                           :: my_f

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

  real, dimension(np,nc)                        :: rhoi_cal
  real, dimension(np)                           :: p_cal, z_cal
  real, dimension(np,nc)                        :: my_cal
  real, dimension(np)                           :: gibbs, f_res
  real                                          :: h_lv
  real, dimension(np)                           :: enthal, entrop, cpres
  real, dimension(np)                           :: speed_of_sound

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
  logical                                       :: check_stability_of_phases = .false.

End Module BASIC_VARIABLES




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains parameters and variables ....
! ................ 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

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
  real                                     :: speed_sound

! ----------------------------------------------------------------------
! quantities of fluid theory
! ----------------------------------------------------------------------
  real                                     :: gij(nc,nc)
  real                                     :: mx(nc,nsite)

  real                                     :: mseg(nc)
  real                                     :: dhs(nc)
  real                                     :: uij(nc,nc)
  real                                     :: sig_ij(nc,nc)

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
  real                                     :: alpha_test(2)

End Module EOS_VARIABLES



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains parameters and variables ....
! ................ 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module EOS_CONSTANTS

  use PARAMETERS, only: nc
  implicit none
  save

  real, dimension(0:6,3)                   :: ap, bp

  real, dimension(28)                      :: c_dd, n_dd, m_dd, k_dd, o_dd
  real, dimension(nc,nc,0:4)               :: qqp2, qqp4, ddp2, ddp4, dqp2, dqp4
  real, dimension(nc,nc,nc,0:4)            :: qqp3, ddp3, dqp3

End Module EOS_CONSTANTS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains parameters and variables ....
! ................ 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

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

  character (LEN=9)                        :: save_eos_terms(12)

End Module EOS_NUMERICAL_DERIVATIVES


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module contains parameters and variables ....
! ................ 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module FITTING_RGT_PARAMETERS

  INTEGER                               :: critfit
  REAL                                  :: tcf, pcf, rcf
  REAL                                  :: LLfit, phifit, chapfit

End Module FITTING_RGT_PARAMETERS
