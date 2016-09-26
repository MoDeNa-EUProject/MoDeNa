!> @file
!! Main program. Homogenization approach to heat transfer in polymer foams
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
!! @page deps Dependencies
!! @section dep_foam_cond  Dependencies of Foam conductivity model
!! - LAPACK 3.4.2 or higher
!! - BLAS 3.4.2 or higher
!! - fson
program hahtf
    use tests, only: loadParameters,eqcond,eqcond_por
    use ioutils, only: newunit
    use constants, only: mfi
    implicit none
    write(*,*) 'Welcome in hahtf'
    open (newunit(mfi),file='hahtf.out')
    write(mfi,*) 'Welcome in hahtf'
    call loadParameters
    call eqcond(1)
    ! call eqcond_por
    ! call eqcond_dcell
    ! call eqcond_strut
    write(*,*) 'Program exited normally'
    write(mfi,*) 'Program exited normally'
    close(mfi)
end program
