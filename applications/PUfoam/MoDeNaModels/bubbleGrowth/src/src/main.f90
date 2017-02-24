!> @file      bubbleGrowth/src/src/main.f90
!! @ingroup   src_mod_bubbleGrowth
!! @author    Pavel Ferkl
!! @brief     Main executable program.
!! @details
!! Simulation of growth of a single bubble in a shell of PU reaction mixture.
!! Summary of the model is provided in \cite Ferkl2016.


!> Calls the appropriate subroutine, depending on what we want to simulate.
!!
!! Changes to this program are for testing purposes only.
program singlebubblegrowth
    use tests, only:onegrowth,secondgrowth,shooting_method,shooting_method_test
    use foaming_globals_m
    implicit none
    firstrun=.true.
    call onegrowth
    ! call secondgrowth
    ! call shooting_method_test
end program
