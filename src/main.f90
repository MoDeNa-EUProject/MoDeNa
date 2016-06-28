!simulation of film drainage between two bubbles
!TODO add condition for film breakage
!TODO connect to bubble growth model
!TODO realistic viscosity evolution
!TODO calculate strut content
program drainage
    use integration, only: preprocess,integrate
    implicit none
    call preprocess
    call integrate
end program
