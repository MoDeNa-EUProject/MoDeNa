!> @file
!! i/o utilities
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module ioutils
    implicit none
    private
    public newunit,str
contains
!********************************BEGINNING*************************************
!> returns lowest i/o unit number not in use
integer function newunit(unit) result(n)
    integer, intent(out), optional :: unit
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
!> converts integer to string
function str(k)
    character(len=20) :: str
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!***********************************END****************************************
end module ioutils
