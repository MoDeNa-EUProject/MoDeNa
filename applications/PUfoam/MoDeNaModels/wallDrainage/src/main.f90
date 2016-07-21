!> @file
!! simulation of film drainage in foam between growing bubbles
!! predicts wall thickness profile and strut size and shape
!! drainage is caused by capillary forces
!! additional stretching is caused by bubble growth
!! @author    Pavel Ferkl
!! @ingroup   wall_drain
!! @page deps Dependencies
!! @section dep_wall_drain  Dependencies of Wall drainage model
!! - bspline-fortran
!! - fson
!TODO add condition for film breakage
program drainage
    use integration, only: preprocess,integrate
    implicit none
    call preprocess
    call integrate
end program
