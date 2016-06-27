!simulation of film drainage between two bubbles
!TODO: calculate minimum film thickness at each time step and its
!distance from center
!TODO: add condition for film breakage
!TODO: connect to bubble growth model
!TODO: realistic viscosity evolution
program drainage
    use integration, only: preprocess,integrate
    implicit none
    call preprocess
    call integrate
end program
