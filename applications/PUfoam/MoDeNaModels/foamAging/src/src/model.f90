!> @file
!! Model subroutines for diffusion of gases in the foam.
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module model
    use globals, only:dp
    implicit none
    integer :: &
        ngas,&    ! number of gases
        nfv     ! number of finite volumes
    integer, dimension(:), allocatable :: &
        mor   ! morphology: 1=cell,2=wall,3=sheet
    real(dp), dimension(:), allocatable :: &
        dz,&    ! size of finite volume
        dif,&   ! diffusivity
        sol,&   ! solubility
        bc      ! boundary condition
    private
    public modelPU,ngas,nfv,dz,dif,sol,bc,mor
contains
!********************************BEGINNING*************************************
!> model supplied to the integrator
subroutine modelPU(neq, time, ystate, yprime)	! ODEPACK call
! 12/07/23 - change of model to describe diffusion on H2 and N2 throug PS, wall
!	is discretized into finite elements, pressure in cell is in equilibrium to
!	the pressure on the wall, solubility given by Henry's law cpol = H*cgas
!	- not ordinary definition of H but the definition of H fits
!	!!!this model is for filling the foam
! 2015/03/12 change in model with respect to ODEPACK calling
! michal.vonka@seznam.cz
! rewritten by pavel.ferkl@vscht.cz
    use constants
    integer :: neq
    real(dp) :: time
    real(dp) :: ystate(neq)
    real(dp) :: yprime(neq)
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
