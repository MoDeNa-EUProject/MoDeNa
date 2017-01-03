!> @file      bubbleGrowth/src/src/model_sundials.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_bubbleGrowth
!! @brief     Interface for sundials.
!! @details
!! Special interface needed for integration by sundials.


!********************************beginning*************************************
!> Model interface for sundials.
!!
!! Must be outside of the module. The name cannot be changed.
!! The name of the model subroutine is hardcoded in here.
subroutine  fcvfun(t, y, ydot, ipar, rpar, ier)
    use iso_c_binding
    use constants, only:dp
    use integration, only:neq
    use model, only: odesystem
    implicit none
    integer :: ier !< error flag
    integer(c_long) :: ipar(1) !< integer parameters
    real(dp) :: t !< time
    real(dp) :: y(neq) !< integrated variables
    real(dp) :: ydot(neq) !< derivatives of integrated variables
    real(dp) :: rpar(1) !< real parameters
    call odesystem(neq, t, y, ydot)
end subroutine fcvfun
!***********************************end****************************************
