!> @file      foamConductivity/src/src/main.f90
!! @ingroup   src_mod_foamConductivity
!! @author    Pavel Ferkl
!! @brief     Main executable program.
!! @details
!! Simulation of conductive-radiative heat transfer in polymer closed-cell foam.
!! Output is the equivalent foam conductivity.


!> Calls the appropriate subroutine, depending on what we want to simulate.
!!
!! Changes to this program are for testing purposes only.
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
