!> @file      wallDrainage/src/globals.f90
!! @ingroup   src_mod_wallDrainage
!! @author    Pavel Ferkl
!! @brief     Global variables.
!! @details
!! Stores variables shared between several modules.
module globals
    use constants
    implicit none
    character(len=1024) :: &
        viscosityModel,&    !< viscosity model (constant, fromFile)
        growthRateModel !< growth rate model (constantGrowth, fromFile)
    integer :: &
        maxts,& !< maximum inner time steps
        its,&   !< outer time steps (how many times is output written)
        meshpoints,&    !< number of discretization points in space
        int_method  !< type of integration method (see MF in ODEPACK)
    real(dp) :: &
        strutFilmParameter,&  !< parameter determining strut-film boundary
        timestep,&  !< time step (how often are values written)
        hi,&    !< film thickness at center
        rd,&    !< computational domain size (radius)
        rc,&    !< total window size (radius)
        rc0,&    !< initial total window size (radius)
        dr,&    !< mesh points spacing
        mu,&    !< viscosity
        gam,&   !< surface tension
        dstr,&  !< where strut starts in initial domain
        ndp,&   !< disjoining pressure constant
        mdp,&   !< disjoining pressure constant
        cdp,&   !< disjoining pressure constant
        hdp,&   !< disjoining pressure constant
        bdp,&   !< disjoining pressure constant, set to zero for no DP
        gr,&    !< growth rate
        int_reltol,&    !< relative tolerance of integrator
        int_abstol,&    !< absolute tolerance of integrator
        initialTime,&   !< starting time
        ae_tol  !< tolerance of algebraic equation solver
end module globals
