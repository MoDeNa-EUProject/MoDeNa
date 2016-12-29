!> @file      bubbleGrowth/src/src/ioutils.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_bubbleGrowth
!! @brief     Tools for file input/output.
!! @details
!! Defines several useful funcions for fortran file i/o.
module ioutils
    implicit none
    private
    public newunit,str
contains
!********************************BEGINNING*************************************
!> Returns lowest i/o unit number not in use.
!!
!! Can be used direcly in the open statement.
integer function newunit(unit) result(n)
    integer, intent(out), optional :: unit !< unit number
    logical inuse
    integer, parameter :: &
        nmin=123,&   ! avoid lower numbers which are sometimes reserved
        nmax=999  ! may be system-dependent
    do n = nmin, nmax
        inquire(unit=n, opened=inuse)
        if (.not. inuse) then
            if (present(unit)) unit=n
            return
        end if
    end do
    write(*,*) "newunit ERROR: available unit not found."
    stop
end function newunit
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Converts integer to string.
!!
!! Uses write function for the conversion.
function str(k)
    character(len=20) :: str
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!***********************************END****************************************
end module ioutils
