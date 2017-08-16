!> @file      foamAging/src/src/constants.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_foamAging
!! @brief     Physical constants and other parameters.
!! @details
!! Names of gases and heat capacity need to be loaded first.
module constants
    use,intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: &
        pi=3.1415926535897932384626433832795028841971693993751058209749445923&
            &078164062862089986280348253421170679_dp,&  !<pi
        sigmaB=5.67037321e-8_dp,&       !<Stefan-Boltzmann constant
        kb=1.380648813e-23_dp,&         !<Boltzmann constant
        Rg=8.31446218_dp,&              !<gas constant
        NA=Rg/kb,&                      !<Avogadro constant
        eulergamma=0.5772156649015328606&
            &065120900824024310421_dp,& !<Euler-Mascheroni constant
        c0=299792458e0_dp,&             !<speed of light in vacuum
        hPc=6.6260695729e-34_dp,&       !<Planck constant
	    C2=0.014387752_dp             !<constant in Planck's law [m*K]
    complex(dp), parameter :: iu=(0.0e0_dp,1.0e0_dp)       !<imaginary constant
    real(dp), dimension(4) :: & !oxygen, nitrogen, carbon dioxide, cyclopentane
        Tc=(/154.58_dp,126.19_dp,304.17_dp,511.7_dp/),& !<critical temperature
        pc=(/5.043e6_dp,3.396e6_dp,7.386e6_dp,45.1e5_dp/),& !<critical pressure
        Mg=(/32e-3_dp,28e-3_dp,44e-3_dp,70e-3_dp/),& !<molar mass
        Tb=(/90.19_dp,77.36_dp,194.75_dp,322.4_dp/),& !<normal boiling point
        cpg=(/0,0,0,0/) !<heat capacity, needs to be filled in conductivity
end module constants
