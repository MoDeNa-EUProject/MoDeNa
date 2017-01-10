!> @file      wallDrainage/src/main.f90
!! @ingroup   src_mod_wallDrainage
!! @author    Pavel Ferkl
!! @brief     Main executable program.
!! @details
!! Simulation of film drainage in foam between growing bubbles.
!! Predicts wall thickness profile and strut size and shape.
!! Drainage is caused by capillary forces.
!! Wall and strut formation is caused by drainage and bubble growth (wall
!! stretching).

!> Calls the appropriate subroutine, depending on what we want to simulate.
!!
!! Changes to this program are for testing purposes only.
!TODO add condition for film breakage
program drainage
    use integration, only: preprocess,integrate
    implicit none
    call preprocess
    call integrate
end program
