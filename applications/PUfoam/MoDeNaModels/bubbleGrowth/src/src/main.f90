!> @file
!! simulation of bubble growth;
!! inspired by [Feng and Bertelo (2004)](http://dx.doi.org/10.1122/1.1645518)
!! and adapted to polyurethanes
!! @author    Pavel Ferkl
!! @ingroup   bblgr
!! @page deps Dependencies
!! @section dep_bblgr  Dependencies of Bubble growth model
!! - LAPACK 3.4.2 or higher
!! - BLAS 3.4.2 or higher
program singlebubblegrowth
    use tests
    use foaming_globals_m
    implicit none
    firstrun=.true.
    call onegrowth
    firstrun=.false.
    allocate(bub_rad(size(init_bub_rad(:,1)),2))
    bub_rad=init_bub_rad
    bub_inx=1
    call onegrowth
end program
