!> @file
!! Main program. Homogenization approach to heat transfer in polymer foams
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
program hahtf
    use tests
    use ioutils
    use constants, only: mfi
    implicit none
    write(*,*) 'Welcome in hahtf'
    open (newunit(mfi),file='hahtf.out')
    write(mfi,*) 'Welcome in hahtf'
    call loadParameters
    call eqcond(1)
!    call eqcond_por
!    call eqcond_dcell
!    call eqcond_strut
    write(*,*) 'Program exited normally'
    write(mfi,*) 'Program exited normally'
    close(mfi)
end program
