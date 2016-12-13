!> @file
!! Calculate dynamics diffusion in the foam structure.
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
!! @page deps Dependencies
!! @section dep_foam_aging  Dependencies of Foam aging model
!! - fson
!c  30.5.2012 - MV (michal.vonka@seznam.cz)
!c	2.7.2012 - MV, remaking to partial pressures of H2 and N2
!c   6.3.2015 - MV, application to PU foams solved by CO2 penetrating air
!   3.4.2015 PF (pavel.ferkl@vscht.cz), calculation of conductivity
program foam_diffusion
    use integration, only: integrate
	implicit none
    print*, 'Welcome to Foam aging'
    call integrate
    print*, 'Program Foam aging finished'
end program
