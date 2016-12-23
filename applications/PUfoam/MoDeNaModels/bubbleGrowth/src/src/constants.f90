!> @file      bubbleGrowth/src/src/constants.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_bubbleGrowth
!! @brief     Stores physical constants and floating point precision.
!! @details
!! Defines only unchangable parameters.
module constants
    use,intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: &
        pi=3.1415926535897932384626433832795028841971693993751058209749445923&
            &078164062862089986280348253421170679_dp,&  !<pi
        sigmaB=5.67037321e-8_dp,&       !<Stefan-Boltzmann constant
        kb=1.380648813e-23_dp,&         !<Boltzmann constant
        Rg=8.31446218_dp,&              !<gas constant
        NA=Rg/kb                        !<Avogadro constant
end module constants
