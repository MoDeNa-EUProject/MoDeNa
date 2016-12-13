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
    integer :: gasModel=2
    private
    public equcond
contains
!********************************BEGINNING*************************************
!> determine equivalent conductivity of the foam
subroutine equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp)
    real(dp), intent(out) :: keq
    real(dp), dimension(:), intent(in) :: ystate
    integer, intent(in) :: ngas,nfv,mor(:)
    real(dp), intent(in) :: temp,eps,dcell,fstrut
    real(dp) :: ccd,cair,ccyp,co2,cn2,xcd,xair,xcyp,xo2,xn2,kgas
    real(dp), dimension(4) :: kg,yg,cpg
    integer :: i,j
    !calculate average concentrations
    co2=0
    cn2=0
    ccd=0
    ccyp=0
    j=1
    do i=1,nfv
        if (mor(i)==1) then
            co2=co2+ystate(ngas*(i-1)+1)
            cn2=cn2+ystate(ngas*(i-1)+2)
            ccd=ccd+ystate(ngas*(i-1)+3)
            ccyp=ccyp+ystate(ngas*(i-1)+4)
            j=j+1
        endif
    enddo
    co2=co2/j
    cn2=cn2/j
    ccd=ccd/j
    ccyp=ccyp/j
    xo2=co2/(co2+cn2+ccd+ccyp)
    xn2=cn2/(co2+cn2+ccd+ccyp)
    xcd=ccd/(co2+cn2+ccd+ccyp)
    xcyp=ccyp/(co2+cn2+ccd+ccyp)
    kg(1)=cdConductivity(temp)
    kg(2)=nitrConductivity(temp)
    kg(3)=oxyConductivity(temp)
    kg(4)=cypConductivity(temp)
    yg=(/xcd,xn2,xo2,xcyp/)
    do i=1,4
        if (abs(yg(i))<1e-6_dp) then
            yg(i)=0
        endif
    enddo
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
    call modena_inputs_set(kfoamInputs, kfoamEpspos, eps)
    call modena_inputs_set(kfoamInputs, kfoamDcellpos, dcell)
    call modena_inputs_set(kfoamInputs, kfoamFstrutpos, fstrut)
    ! call modena_inputs_set(kfoamInputs, kfoamKgaspos, kgas)
    call modena_inputs_set(kfoamInputs, kfoamTemppos, temp)
    call modena_inputs_set(kfoamInputs, kfoamXCO2pos, yg(1))
    call modena_inputs_set(kfoamInputs, kfoamXCyPpos, yg(4))
    call modena_inputs_set(kfoamInputs, kfoamXO2pos, yg(3))
    call modena_inputs_set(kfoamInputs, kfoamXN2pos, yg(2))
    ret = modena_model_call (kfoamModena, kfoamInputs, kfoamOutputs)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    keq = modena_outputs_get(kfoamOutputs, 0_c_size_t); !fetch results
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
    if (minval(y)<0) then
        print*,  'Input molar fractions to weightedAverage cannot be negative.'
        print*, y
        stop
    endif
    y=y/sum(y)
    kmix=0
    do i=1,size(k)
        kmix=kmix+yin(i)*k(i)
    enddo
end function weightedAverage
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine thermal conductivity of a mixture
!! extended Wassiljewa model (Dohrn)
!! parameters calculated according to Mason and Saxena
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
    if (minval(y)<0) then
        print*,  'Input molar fractions to extWassiljewa cannot be negative.'
        print*, y
        stop
    endif
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
    if (minval(y)<0) then
        print*,  'Input molar fractions to lindsayBromley cannot be negative.'
        print*, y
        stop
    endif
    y=y/sum(y)
    do i=1,n
        S(i)=1.5_dp*Tb(i)
    enddo
    cv=cp-Rg    !ideal gas assumed
    gam=cp/cv
    do i=1,n
        do j=1,n
            x=k(i)/k(j)*cp(j)/cp(i)*(9-5/gam(j))/(9-5/gam(i))
            A(i,j)=&
                0.25_dp*(1+sqrt(x*(M(j)/M(i))**0.75_dp*(T+S(i))/(T+S(j))))**2*&
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
    if (minval(y)<0) then
        print*,  'Input molar fractions to pandeyPrajapati cannot be negative.'
        print*, y
        stop
    endif
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
