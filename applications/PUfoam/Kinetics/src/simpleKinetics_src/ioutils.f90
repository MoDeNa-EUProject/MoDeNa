!******************************************************BEGINNING***************************************************************
!stores parameters and commonly used variables
module ioutils
!*****************************************************DECLARATION**************************************************************
    implicit none
    private
    public newunit,str
!*********************************************************BODY*****************************************************************
contains
!****************************BEGINNING*******************************
! returns lowest i/o unit number not in use
integer function newunit(unit) result(n)
!***************************DECLARATION******************************
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
!******************************BODY**********************************
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
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
! converts integer to string
character(len=20) function str(k)
!***************************DECLARATION******************************
    integer, intent(in) :: k
!******************************BODY**********************************
    write (str, *) k
    str = adjustl(str)
end function str
!********************************************************************
!*******************************END**********************************
end module ioutils
!**********************************************************END*****************************************************************
