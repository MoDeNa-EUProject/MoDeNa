!> @file
!! simulation of bubble growth;
!! summary in [Ferkl et al., 2016](dx.doi.org/10.1016/j.ces.2016.03.040)
!! @author    Pavel Ferkl
!! @ingroup   bblgr
!! @page deps Dependencies
!! @section dep_bblgr  Dependencies of Bubble growth model
!! - LAPACK 3.4.2 or higher
!! - BLAS 3.4.2 or higher
!! - bspline-fortran
!! - fson
!! - sundials
program singlebubblegrowth
    use tests, only:onegrowth,secondgrowth,shooting_method,shooting_method_test
    use foaming_globals_m
    implicit none
    firstrun=.true.
    call onegrowth
    ! call secondgrowth
    ! call shooting_method_test
end program
