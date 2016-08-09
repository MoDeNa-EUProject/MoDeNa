!> @file
!! global variables
!! @author    Pavel Ferkl
!! @ingroup   wall_drain
module globals
    use constants
    implicit none
    character(len=1024) :: &
        viscosityModel,&    !viscosity model (constant, fromFile)
        growthRateModel !growth rate model (constantGrowth, fromFile)
    integer :: &
        maxts,& !maximum inner time steps
        its,&   !outer time steps (how many times is output written)
        meshpoints,&    !number of discretization points in space
        int_method  !type of integration method (see MF in ODEPACK)
    real(dp) :: &
        strutFilmParameter,&  !parameter determining strut-film boundary
            !film ends where h(r)>strutFilmParameter*h(r=0)
        timestep,&  !time step (how often are values written)
        hi,&    !film thickness at center
        rd,&    !computational domain size (radius)
        rc,&    !total window size (radius)
        rc0,&    !initial total window size (radius)
        dr,&    !mesh points spacing
        mu,&    !viscosity
        gam,&   !surface tension
        dstr,&  !where strut starts in initial domain
        ndp,&   !disjoining pressure constant
        mdp,&   !disjoining pressure constant
        cdp,&   !disjoining pressure constant
        hdp,&   !disjoining pressure constant
        bdp,&   !disjoining pressure constant, set to zero for no disjoining
            !pressure
        gr,&    !growth rate
        int_reltol,&    !relative tolerance of integrator
        int_abstol,&    !absolute tolerance of integrator
        initialTime,&   !starting time
        ae_tol  !tolerance of algebraic equation solver
end module globals
