!global variables
module globals
    use constants
    implicit none
    integer :: &
        maxts,& !maximum inner time steps
        its,&   !outer time steps (how many times is output written)
        meshpoints,&    !number of discretization points in space
        int_method  !type of integration method (see MF in ODEPACK)
    real(dp) :: &
        timestep,&  !time step (how often are values written)
        hi,&    !film thickness at center
        rs,&    !radius of initial strut curvature
        rd,&    !computational domain size (radius)
        rc,&    !total window size (radius)
        rc0,&    !initial total window size (radius)
        dr,&    !mesh points spacing
        s,& !film thickness derivative at outer domain boundary
        q,& !flux into domain from strut
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
