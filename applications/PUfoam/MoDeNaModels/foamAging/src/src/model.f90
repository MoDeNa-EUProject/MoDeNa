!> @file      foamAging/src/src/model.f90
!! @ingroup   src_mod_foamAging
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @brief     Model subroutines.
!! @details
!! Implementation of the physical equations.
module model
    use globals, only:dp
    implicit none
    integer :: &
        ngas,&    !< number of gases
        nfv     !< number of finite volumes
    integer, dimension(:), allocatable :: &
        mor   !< morphology: 1=cell,2=wall,3=sheet
    real(dp), dimension(:), allocatable :: &
        dz,&    !< size of finite volume
        dif,&   !< diffusivity
        sol,&   !< solubility
        bc      !< boundary condition
    private
    public modelPU,ngas,nfv,dz,dif,sol,bc,mor
contains
!********************************BEGINNING*************************************
!> Model supplied to the integrator.
!!
!! Model to describes diffusion of gases through polymer. Wall
!! is discretized into finite elements, pressure in cell is in equilibrium to
!! the pressure on the wall, solubility given by Henry's law cpol = H*cgas.
!! Model call is defined for ODEPACK.
!! @param [in] time time
subroutine modelPU(neq, time, ystate, yprime)
    use constants
    integer, intent(in) :: neq !< number of equations
    real(dp), intent(in) :: time !< time
    real(dp), intent(in) :: ystate(neq) !< integrated variables
    real(dp), intent(out) :: yprime(neq) !< derivatives of integrated variables
    integer :: i, j, k
    real(dp) :: fluxw,fluxe
    k=1
    j=1
    do i=1,ngas
        fluxw=-2*dif(k)*sol(k)*(ystate(k)-bc(i))/dz(j)
        fluxe=-2*dif(k)*dif(k+ngas)*sol(k)*sol(k+ngas)*&
            (ystate(k+ngas)-ystate(k))/&
            (dif(k+ngas)*dz(j)*sol(k)+dif(k)*dz(j+1)*sol(k+ngas))
        yprime(k)=(fluxw-fluxe)/dz(j)
        k=k+1
    enddo
    do j=2,nfv-1
        do i=1,ngas
            fluxw=-2*dif(k-ngas)*dif(k)*sol(k)*sol(k-ngas)*&
                (ystate(k)-ystate(k-ngas))/&
                (dif(k)*dz(j-1)*sol(k-ngas)+dif(k-ngas)*dz(j)*sol(k))
            fluxe=-2*dif(k)*dif(k+ngas)*sol(k)*sol(k+ngas)*&
                (ystate(k+ngas)-ystate(k))/&
                (dif(k+ngas)*dz(j)*sol(k)+dif(k)*dz(j+1)*sol(k+ngas))
            yprime(k)=(fluxw-fluxe)/dz(j)
            k=k+1
        enddo
    enddo
    j=nfv
    do i=1,ngas
        fluxw=-2*dif(k-ngas)*dif(k)*sol(k)*sol(k-ngas)*&
            (ystate(k)-ystate(k-ngas))/&
            (dif(k)*dz(j-1)*sol(k-ngas)+dif(k-ngas)*dz(j)*sol(k))
        fluxe=0
        yprime(k)=(fluxw-fluxe)/dz(j)
        k=k+1
    enddo
end subroutine modelPU
!***********************************END****************************************
end module model
