!> @file
!! stores parameters and commonly used variables as globals
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
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
	    C2=0.014387752_dp,&             !<constant in Planck law [m*K]
        struttol=1e-3_dp                !<tolerance for strut content
                                        !<foam has no struts if fs<struttol
                                        !<foam has no walls if fs>1-struttol
    complex(dp), parameter :: iu=(0.0e0_dp,1.0e0_dp)       !<imaginary constant
    logical  :: &
        wdist,&               !<use wall thickness distribution
        testMode,&            !<true disables calculation of radiation
        numcond               !<calcualte effective conductivity numerically
    character(len=80) :: &
        structureName         !<name of the file with morphology
    integer  :: &
        mfi,&                 !<main file index
        nrays,&               !<number of testing rays
        nz,&                  !<spatial discretization
        nbox,&                !<number of gray boxes
        morph_input           !<morphology input
                              !<1=wall thickness, 2=strut content,
                              !<3=strut diameter, 4==strut content (alternative)
                              !<(3 is recommended others can
                              !<have multiple solutions)
    real(dp) :: lambda,&      !<wavelength
        temp1,temp2,&         !<temperatures at boundaries
        kappa2,&              !<absorption coefficient of solid
        cond1,cond2,&         !<conductivities
        n1,n2,&               !<real part of refractive indices
        k1,k2,&               !<imaginary part of refractive indices
        rhog,rhos,&           !<densities of gas and solid
        emi1,emi2,&           !<wall emittances
        tmean,&               !<mean temperature
        eqc,&                 !<equivalent conductivity
        eqc_ross,&            !<Rosseland equivalent conductivity
        effc,&                !<effective conductivity (only conduction)
        effc_num,&            !<effective conductivity from numerical simulation
        kgas,&                !<gas conductivity
        ksol,&                !<solid conductivity
        krad,&                !<radiative conductivity
        gcontr,&              !<contribution of gas
        scontr,&              !<contribution of solid
        rcontr,&              !<contribution of radiation
        planckextcoeff,&      !<Planck extinction coefficient
        rossextcoeff,&        !<Rosseland extinction coefficient
        albedo,&              !<scattering albedo
        effn,&                !<effective index of refraction (real part)
        por,&                 !<porosity
        rhof,&                !<foam density
        dcell,&               !<cell size
        dwall,&               !<wall thickness
        fs,&                  !<strut content
        dstrut,&              !<strut diameter
        dfoam,&               !<foam thickness
        wsdev,&               !<wall thickness standard deviation
        tm(10)                !<for time measurements
    real(dp), dimension(:), allocatable :: &
        lambdan,nwl,&       !<wavelength and real part of refractive index
        lambdak,kwl,&       !<wavelength and imaginary part of refractive index
        lambdagas,acgas,&   !<wavelength and absorption coefficient of gas
        lambdabox,&         !<boundaries of gray boxes
        trextcoeffbox,&     !<transport extinction coefficients of gray boxes
        albedobox,&         !<scattering albedos of gray boxes
        abscoeffbox,&       !<absorption coefficients of gray boxes
        scattcoeffbox,&     !<scatering coefficients of gray boxes
        fbepbox             !<fraction of blackbody radiation in the box
end module constants
