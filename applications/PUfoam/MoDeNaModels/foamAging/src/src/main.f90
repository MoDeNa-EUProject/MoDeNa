!> @file      foamAging/src/src/main.f90
!! @ingroup   src_mod_foamAging
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @brief     Main executable program.
!! @details
!! Simulation of gas diffusion through heterogeneous material. Originally
!! developed for permeation of H2 and N2 through membranes. Later rewritten for
!! foams and different type of gases.


!> Calls the appropriate subroutine, depending on what we want to simulate.
!!
!! Changes to this program are for testing purposes only.
program foam_diffusion
    use integration, only: integrate
	implicit none
    print*, 'Welcome to Foam aging'
    call integrate
    print*, 'Program Foam aging finished'
end program
