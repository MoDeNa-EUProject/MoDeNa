module foaming_globals_m

  ! use kind_defs_m
  ! use constants

  implicit none

  logical :: firstrun=.true.       ! running for the first time
  integer :: bub_inx              ! bubble index
  real(8) :: TSTART,           & ! starting time of the direct simulations
              TEND,             & ! end time of the direct simulations
              bub_pres,         & ! bubble pressure
              surface_tension     ! surface tension
  real(8), dimension(:,:), allocatable :: bub_rad,& ! first column represents
                                                     ! time and the n the radius
                                                     ! of each bubble every timestep
                                           etat,& ! first column represents
                                                  ! time, second column viscosity
                                                  ! for each timestep
                                          port,& ! first column represents
                                                 ! time, second column porosity
                                                 ! for each timestep
                                         init_bub_rad ! first column represents
                                                ! time, second column radius
                                                ! for each timestep


end module foaming_globals_m
