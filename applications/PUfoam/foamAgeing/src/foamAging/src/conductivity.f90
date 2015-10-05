!> @file
!! subroutines for calculation of equivalent conductivity of the foam
!! using Modena calls
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module conductivity
    use constants
    use fmodena
    use physicalProperties
    implicit none
    private
    public equcond,mixtureConductivity
contains
!********************************BEGINNING*************************************
!> determine equivalent conductivity of the foam
subroutine equcond(keq,ystate,neq,eps,fstrut,temp)
    real(dp), intent(out) :: keq
    real(dp), dimension(:), intent(in) :: ystate
    integer, intent(in) :: neq
    real(dp), intent(in) :: temp,eps,fstrut
    real(dp) :: dcell,ccd,cair,ccyp,xcd,xair,xcyp,kgas
    real(dp), dimension(3) :: kg,yg
    integer :: i,ncell,onecell,nFV
    dcell   = ystate(nEQ + 1 )
    ncell = int(ystate(nEQ + 11))
    onecell = int(ystate(nEQ + 12)) != dble(ipar(4))
    nFV  = onecell*ncell
    !calculate average concentrations
    ccd=0
    cair=0
    ccyp=0
    do i=1,ncell
        cair=cair+ystate(i*onecell)
        ccd=ccd+ystate(nFV+i*onecell)
        ccyp=ccyp+ystate(2*nFV+i*onecell)
    enddo
    cair=cair/ncell
    ccd=ccd/ncell
    ccyp=ccyp/ncell
    xair=cair/(cair+ccd+ccyp)
    xcd=ccd/(cair+ccd+ccyp)
    xcyp=ccyp/(cair+ccd+ccyp)
    kg(1)=cdConductivity(temp)
    kg(2)=airConductivity(temp)
    kg(3)=cypConductivity(temp)
    yg=(/xcd,xair,xcyp/)
    call mixtureConductivity(kgas,kg,yg,Tc,pc,Mg,temp,1._dp)
    ! write(*,*) yg
    ! write(*,*) kg
    ! write(*,*) kgas
    kgas=0
    do i=1,size(kg)
        kgas=kgas+yg(i)*kg(i)
    enddo
    ! write(*,*) yg
    ! write(*,*) kg
    ! write(*,*) kgas
    call modena_inputs_set(kfoamInputs, kfoamEpspos, eps)
    call modena_inputs_set(kfoamInputs, kfoamDcellpos, dcell)
    call modena_inputs_set(kfoamInputs, kfoamFstrutpos, fstrut)
    call modena_inputs_set(kfoamInputs, kfoamKgaspos, kgas)
    call modena_inputs_set(kfoamInputs, kfoamTemppos, temp)
    ret = modena_model_call (kfoamModena, kfoamInputs, kfoamOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    keq = modena_outputs_get(kfoamOutputs, 0_c_size_t); !fetch results
end subroutine equcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine thermal conductivity of a mixture
!! Wassiljewa model, parameters calculated according to Mason and Saxena
!! [link](http://dx.doi.org/10.1016/j.fluid.2007.07.059)
subroutine mixtureConductivity(kmix,k,yin,Tc,pc,M,T,eps)
    real(dp), intent(out) :: kmix   !thermal conductivity of the mixture
    real(dp), dimension(:), intent(in) :: k !thermal conductivities
    real(dp), dimension(:), intent(in) :: yin !molar fractions
    real(dp), dimension(:), intent(in) :: Tc !critical temperatures
    real(dp), dimension(:), intent(in) :: pc !critical pressures
    real(dp), dimension(:), intent(in) :: M !molar masses
    real(dp), intent(in) :: T !temperature
    real(dp), intent(in) :: eps !parameter close to one
    integer :: n,i,j
    real(dp) :: x
    real(dp), dimension(:), allocatable :: gam,y
    real(dp), dimension(:,:), allocatable :: ktr,A
    n=size(k)
    allocate(y(n),gam(n),ktr(n,n),A(n,n))
    y=yin
    if (minval(y)<0) stop 'Input molar fractions to mixtureConductivity &
        cannot be negative.'
    y=y/sum(y)
    gam=210*(Tc*M**3/pc**4)**(1/6._dp)
    do i=1,n
        do j=1,n
            ktr(i,j)=gam(j)*(exp(0.0464_dp*T/Tc(i))-exp(-0.2412_dp*T/Tc(i)))/&
                gam(i)/(exp(0.0464_dp*T/Tc(j))-exp(-0.2412_dp*T/Tc(j)))
            A(i,j)=eps*(1+sqrt(ktr(i,j))*(M(i)/M(j))**0.25_dp)**2/&
                sqrt(8*(1+M(i)/M(j)))
        enddo
    enddo
    kmix=0
    do i=1,n
        x=0
        do j=1,n
            x=x+y(j)*A(i,j)
        enddo
        kmix=kmix+y(i)*k(i)/x
    enddo
end subroutine mixtureConductivity
!***********************************END****************************************
end module conductivity
