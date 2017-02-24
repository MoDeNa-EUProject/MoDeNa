!> @file      bubbleGrowth/src/src/foaming_globals.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_bubbleGrowth
!! @brief     Stores global variables for Foaming application.
!! @details
!! Common variables with the application for the simulation of a box with
!! several growing bubbles.
module foaming_globals_m
  implicit none
  logical :: &
    firstrun,&       !< running for the first time
    printout = .true.  ! print intermediate results during integration in bblgr
  integer :: bub_inx             !< bubble index
  real(8) :: tstart,           & !< starting time of the direct simulations
             tend,             & !< end time of the direct simulations
             bub_pres,         & !< bubble pressure
             surface_tension     !< surface tension
  real(8), dimension(:,:), allocatable :: &
    bub_rad,& !< first column - time; rest radius of each bubble every timestep
    etat,& !< first column - time, second column viscosity for each timestep
    port,& !< first column - time, second column porosity for each timestep
    init_bub_rad !< first column - time, second column radius for each timestep
end module foaming_globals_m
