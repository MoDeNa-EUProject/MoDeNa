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
    integer :: gasModel=3
    private
    public equcond
contains
!********************************BEGINNING*************************************
!> determine equivalent conductivity of the foam
subroutine equcond(keq,ystate,neq,eps,fstrut,temp)
    real(dp), intent(out) :: keq
    real(dp), dimension(:), intent(in) :: ystate
    integer, intent(in) :: neq
    real(dp), intent(in) :: temp,eps,fstrut
    real(dp) :: dcell,ccd,cair,ccyp,xcd,xair,xcyp,kgas
    real(dp), dimension(4) :: kg,yg,cpg
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
    kg(2)=nitrConductivity(temp)
    kg(3)=oxyConductivity(temp)
    kg(4)=cypConductivity(temp)
    yg=(/xcd,0.79_dp*xair,0.21_dp*xair,xcyp/)
    ! write(*,*) yg
    ! write(*,*) 'x argPos:', kfoamXCO2pos, kfoamXCyPpos, kfoamXO2pos, kfoamXN2pos
    ! yg=(/0.0_dp,0.79_dp*0.5_dp,0.21_dp*0.5_dp,0.5_dp/)
    cpg(1)=cdHeatCapacity(temp)
    cpg(2)=nitrHeatCapacity(temp)
    cpg(3)=oxyHeatCapacity(temp)
    cpg(4)=cypHeatCapacity(temp)
    select case(gasModel) ! determine conductivity of gas mixture
    case (1)
        kgas = weightedAverage(kg,yg)
    case (2)
        kgas = extWassiljewa(kg,yg,Tc,pc,Mg,temp,1._dp)
    case (3)
        kgas = lindsayBromley(kg,yg,Tb,cpg,Mg,temp)
    case (4)
        kgas = pandeyPrajapati(kg,yg,Tb,Mg,temp)
    case default
        stop 'unknown gas conductivity model'
    end select
    ! write(*,*) yg
    ! write(*,*) kg
    ! write(*,*) "kgas: ", kgas
    ! stop
    ! write(*,*) 'eps: ',eps
    ! write(*,*) 'dcell: ',dcell
    ! write(*,*) 'fstrut: ',fstrut
    ! write(*,*) 'temp: ',temp
    ! write(*,*) 'xcd: ',xcd
    ! write(*,*) 'xair: ',xair
    ! write(*,*) 'xcyp: ',xcyp
    call modena_inputs_set(kfoamInputs, kfoamEpspos, eps)
    call modena_inputs_set(kfoamInputs, kfoamDcellpos, dcell)
    call modena_inputs_set(kfoamInputs, kfoamFstrutpos, fstrut)
    ! call modena_inputs_set(kfoamInputs, kfoamKgaspos, kgas)
    call modena_inputs_set(kfoamInputs, kfoamTemppos, temp)
    call modena_inputs_set(kfoamInputs, kfoamXCO2pos, xcd)
    call modena_inputs_set(kfoamInputs, kfoamXCyPpos, xCyP)
    call modena_inputs_set(kfoamInputs, kfoamXO2pos, xAir*0.21_dp)
    call modena_inputs_set(kfoamInputs, kfoamXN2pos, xAir*0.79_dp)
    ret = modena_model_call (kfoamModena, kfoamInputs, kfoamOutputs)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    keq = modena_outputs_get(kfoamOutputs, 0_c_size_t); !fetch results
    ! write(*,*) keq
    ! stop
end subroutine equcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine thermal conductivity of a mixture
!! simple weighted average
real(dp) function weightedAverage(k,yin) result(kmix)
    real(dp), dimension(:), intent(in) :: k !thermal conductivities
    real(dp), dimension(:), intent(in) :: yin !molar fractions
    integer :: n,i,j
    real(dp), dimension(:), allocatable :: y
    n=size(k)
    allocate(y(n))
    y=yin
    if (minval(y)<0) stop 'Input molar fractions to weightedAverage &
        cannot be negative.'
    y=y/sum(y)
    kmix=0
    do i=1,size(k)
        kmix=kmix+yin(i)*k(i)
    enddo
end function weightedAverage
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine thermal conductivity of a mixture
!! extended Wassiljewa model, parameters calculated according to Mason and Saxena
!! [link](http://dx.doi.org/10.1016/j.fluid.2007.07.059)
real(dp) function extWassiljewa(k,yin,Tc,pc,M,T,eps) result(kmix)
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
    if (minval(y)<0) stop 'Input molar fractions to extWassiljewa &
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
end function extWassiljewa
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine thermal conductivity of a mixture
!! Lindsay-Bromley model
!! see [link](http://dx.doi.org/10.1021/ie50488a017)
real(dp) function lindsayBromley(k,yin,Tb,cp,M,T) result(kmix)
    real(dp), dimension(:), intent(in) :: &
        k,& !thermal conductivities
        yin,& !molar fractions
        Tb,& !boiling point temperatures
        cp,& !thermal capacities at constant pressure
        M !molar masses
    real(dp), intent(in) :: &
        T !temperature
    integer :: n,i,j
    real(dp) :: x
    real(dp), dimension(:), allocatable :: &
        y,& !molar fractions
        S,& !Sutherland constants
        gam,& !heat capacity ratio
        cv !thermal capacities at constant volume
    real(dp), dimension(:,:), allocatable :: A
    n=size(k)
    allocate(y(n),cv(n),S(n),gam(n),A(n,n))
    y=yin
    if (minval(y)<0) stop 'Input molar fractions to lindsayBromley &
        cannot be negative.'
    y=y/sum(y)
    do i=1,n
        S(i)=1.5_dp*Tb(i)
    enddo
    cv=cp-Rg    !ideal gas assumed
    gam=cp/cv
    do i=1,n
        do j=1,n
            x=k(i)/k(j)*cp(j)/cp(i)*(9-5/gam(j))/(9-5/gam(i))
            A(i,j)=0.25_dp*(1+sqrt(x*(M(j)/M(i))**0.75_dp*(T+S(i))/(T+S(j))))**2*&
                (T+sqrt(S(i)*S(j)))/(T+S(j))
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
end function lindsayBromley
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine thermal conductivity of a mixture
!! Pandey-Prajapati model
!! see [link](http://www.new1.dli.ernet.in/data1/upload/insa/INSA_1/20005baf_372.pdf)
real(dp) function pandeyPrajapati(k,yin,Tb,M,T) result(kmix)
    real(dp), dimension(:), intent(in) :: &
        k,& !thermal conductivities
        yin,& !molar fractions
        Tb,& !boiling point temperatures
        M !molar masses
    real(dp), intent(in) :: &
        T !temperature
    integer :: n,i,j
    real(dp) :: x
    real(dp), dimension(:), allocatable :: &
        y,& !molar fractions
        S !Sutherland constants
    real(dp), dimension(:,:), allocatable :: A
    n=size(k)
    allocate(y(n),S(n),A(n,n))
    y=yin
    if (minval(y)<0) stop 'Input molar fractions to pandeyPrajapati &
        cannot be negative.'
    y=y/sum(y)
    do i=1,n
        S(i)=1.5_dp*Tb(i)
    enddo
    do i=1,n
        do j=1,n
            A(i,j)=0.25_dp*(1+sqrt(k(i)/k(j)*(M(i)/M(j))**0.25_dp*&
                (T+S(i))/(T+S(j))))**2*(T+sqrt(S(i)*S(j)))/(T+S(j))
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
end function pandeyPrajapati
!***********************************END****************************************
end module conductivity
