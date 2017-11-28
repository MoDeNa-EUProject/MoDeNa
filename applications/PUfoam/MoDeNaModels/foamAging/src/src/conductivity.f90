!> @file      foamAging/src/src/conductivity.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_foamAging
!! @brief     Calculation of thermal conductivity.
!! @details
!! Subroutines for prediction of thermal conductivity of foam and gas phase.
module conductivity
    use constants
    use fmodena
    use physicalProperties
    use globals, only: Mg
    implicit none
    private
    public equcond
contains
!********************************BEGINNING*************************************
!> Determines equivalent conductivity of the foam.
!!
!! Calls the Modena model. Also tests different mixing rules for gas mixture
!! conductivity.
subroutine equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp)
    real(dp), intent(out) :: keq !< equivalent conductivity of foam
    real(dp), dimension(:), intent(in) :: ystate !< integrated variables
    integer, intent(in) :: ngas !< number of gases
    integer, intent(in) :: nfv !< number of grid points
    integer, intent(in) :: mor(:) !< phase information
    real(dp), intent(in) :: temp !< temperature
    real(dp), intent(in) :: eps !< porosity
    real(dp), intent(in) :: dcell !< cell size
    real(dp), intent(in) :: fstrut !< strut content
    real(dp) :: kgas
    real(dp), dimension(ngas) :: kg,yg,cc
    integer :: i,j,k
    integer :: gasModel=2
    !calculate average concentrations
    cc=0
    j=0
    do i=1,nfv
        if (mor(i)==1) then
            do k=1,ngas
                cc(k)=cc(k)+ystate(ngas*(i-1)+k)
            enddo
            j=j+1
        endif
    enddo
    cc=cc/j
    yg=cc/sum(cc)
    do i=1,ngas
        kg(i)=gasConductivity(temp,i)
        if (yg(i)<1e-6_dp) then
            yg(i)=0
        endif
    enddo
    call modena_inputs_set(kfoamInputs, kfoamEpspos, eps)
    call modena_inputs_set(kfoamInputs, kfoamDcellpos, dcell)
    call modena_inputs_set(kfoamInputs, kfoamFstrutpos, fstrut)
    call modena_inputs_set(kfoamInputs, kfoamTemppos, temp)
    do i=1,ngas
        call modena_inputs_set(kfoamInputs, kfoamXg(i), yg(i))
    enddo
    ret = modena_model_call (kfoamModena, kfoamInputs, kfoamOutputs)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    keq = modena_outputs_get(kfoamOutputs, 0_c_size_t); !fetch results
end subroutine equcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Determine thermal conductivity of a mixture.
!!
!! Simple weighted average.
real(dp) function weightedAverage(k,yin) result(kmix)
    real(dp), dimension(:), intent(in) :: k !thermal conductivities
    real(dp), dimension(:), intent(in) :: yin !molar fractions
    integer :: n,i,j
    real(dp), dimension(:), allocatable :: y
    n=size(k)
    allocate(y(n))
    y=yin
    if (minval(y)<0) then
        print*, 'Input molar fractions to weightedAverage cannot be negative.'
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
!> Determine thermal conductivity of a mixture.
!!
!! Extended Wassiljewa model (Dohrn).
!! Parameters calculated according to Mason and Saxena
!! [link](http://dx.doi.org/10.1016/j.fluid.2007.07.059).
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
        print*, 'Input molar fractions to extWassiljewa cannot be negative.'
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
!> Determine thermal conductivity of a mixture
!!
!! Lindsay-Bromley model,
!! see [link](http://dx.doi.org/10.1021/ie50488a017).
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
        print*, 'Input molar fractions to lindsayBromley cannot be negative.'
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
!> Determine thermal conductivity of a mixture.
!!
!! Pandey-Prajapati model,
!! see [link](http://www.new1.dli.ernet.in/data1/upload/insa/INSA_1/20005baf_372.pdf).
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
        print*, 'Input molar fractions to pandeyPrajapati cannot be negative.'
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
