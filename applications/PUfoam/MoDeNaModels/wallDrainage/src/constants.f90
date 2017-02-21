!> @file      wallDrainage/src/constants.f90
!! @ingroup   src_mod_wallDrainage
!! @author    Pavel Ferkl
!! @brief     Physical constants.
!! @details
!! Stores physical constants and type definitions.
module constants
    use,intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi=3.141592653589793238462_dp !< pi
    real(dp), parameter :: sigmaB=5.67037321e-8_dp       !< Stefan-Boltzmann constant
    real(dp), parameter :: kb=1.380648813e-23_dp         !< Boltzmann constant
    real(dp), parameter :: Rg=8.31446218_dp              !< gas constant
    real(dp), parameter :: NA=Rg/kb                      !< Avogadro constant
end module constants
