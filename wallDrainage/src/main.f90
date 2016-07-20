!simulation of film drainage in foam between growing bubbles
!TODO add condition for film breakage
!TODO connect to bubble growth model
!TODO realistic viscosity evolution
program drainage
    use integration, only: preprocess,integrate
    implicit none
    call preprocess
    call integrate
end program
