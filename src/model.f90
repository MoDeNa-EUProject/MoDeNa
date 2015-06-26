!******************************************************BEGINNING***************************************************************
!contains model subroutines
!TODO: material balance (to check results)
!TODO: dimensional results
!TODO: complete flow boundary conditions in several models
!TODO: add condition for film breakage when minimum film thickness is below certain threshold
module model
!*****************************************************DECLARATION**************************************************************
    use constants
    use ioutils
    implicit none

    integer :: maxts=5000
    integer :: its=100
    real(dp) :: timestep=1e-3_dp
    real(dp) :: hi=5e-6_dp
    real(dp) :: Rc=1e-4_dp
    real(dp) :: s=1/sqrt(3._dp)
    real(dp) :: Q=0*2e-1_dp
    real(dp) :: mu=1e-1_dp
    real(dp) :: gam=2e-2_dp
    real(dp) :: dstr=0.2_dp

    real(dp) :: dr
    real(dp) :: Rs
    real(dp), dimension(:), allocatable :: R,U

    !time integration variables for lsode
    integer :: IOUT
    real(dp), dimension(:), allocatable :: RWORK,Y
    integer, dimension(:), allocatable :: IWORK
    real(dp) :: JAC,TOUT,RTOL=1e-12_dp,ATOL=1e-12_dp,T=0
    integer :: IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ=400, NNZ, LENRAT, MF=222

    !needed for selection of model subroutine
    abstract interface
        subroutine sub (NEQ, T, Y, YDOT)
            use constants
            INTEGER :: NEQ
            real(dp) ::  T, Y(NEQ), YDOT(NEQ)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => FEX8!null ()
!*********************************************************BODY*****************************************************************
contains
!****************************BEGINNING*******************************
!FVM, equidistant mesh
!cylindrical geometry, sin(alpha)=(dh/dr)/(1+(dh/dr)**2)
subroutine  FEX8 (NEQ, T, Y, YDOT)
!***************************DECLARATION******************************
    integer :: NEQ,i
    real(dp) :: T, Y(NEQ), YDOT(NEQ)
    real(dp) :: z,ze,zw,zee,zww
    real(dp) :: lame,lamw
    real(dp) :: h,he,hw,hee,hww,heee,hwww
    real(dp) :: he1,hw1,he2,hw2,he3,hw3
    real(dp) :: fluxe,fluxw
!******************************BODY**********************************
    dr=(Rc-Y(NEQ))/NEQ
!    dr=Rc/NEQ
    do i=1,NEQ
        if (i==1) then
            z=dr/2
            ze=dr
            zee=3*dr/2
            lame=(ze-z)/(zee-z)
            h=Y(i)
            hee=Y(i+1)
            hw=h
            he=hee*lame+h*(1-lame)
            heee=Y(i+2)*lame+hee*(1-lame)
            he1=(hee-h)/(dr)
            he2=(hee-2*he+h)/(dr**2/4)
            he3=(heee-2*hee+2*h-hw)/(dr**3/4)
            fluxw=0
            fluxe=he**3*((2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp)*gam
        elseif (i==NEQ) then
            zww=z
            zw=ze
            z=zee
            lamw=(zw-zww)/(z-zww)
            hww=Y(i-1)
            h=Y(i)
            hw=h*lamw+hww*(1-lamw)
            he=h+dr*s/2
            hwww=hww*lamw+Y(i-2)*(1-lamw)
            hw1=(h-hww)/(dr)
            hw2=(h-2*hw+hww)/(dr**2/4)
            hw3=(he-2*h+2*hww-hwww)/(dr**3/4)
            fluxw=hw**3*((2*hw1**3/zw+hw1**5/zw+hw1*(1+3*zw**2*hw2**2)/zw-hw2-zw*hw3-hw1**2*(hw2+zw*hw3))/(1+hw1**2)**2.5_dp)*gam
            fluxe=-Q
        else
            zww=z
            zw=ze
            z=zee
            ze=ze+dr
            zee=ze+dr/2
            lamw=(zw-zww)/(z-zww)
            lame=(ze-z)/(zee-z)
            hww=Y(i-1)
            h=Y(i)
            hee=Y(i+1)
            hw=h*lamw+hww*(1-lamw)
            he=hee*lame+h*(1-lame)
            if (i==2) then
                hwww=hww
            else
                hwww=hww*lamw+Y(i-2)*(1-lamw)
            endif
            if (i==NEQ-1) then
                heee=hee+dr*s/2
            else
                heee=Y(i+2)*lame+hee*(1-lame)
            endif
            hw1=(h-hww)/(dr)
            hw2=(h-2*hw+hww)/(dr**2/4)
            hw3=(he-2*h+2*hww-hwww)/(dr**3/4)
            he1=(hee-h)/(dr)
            he2=(hee-2*he+h)/(dr**2/4)
            he3=(heee-2*hee+2*h-hw)/(dr**3/4)
            fluxw=hw**3*((2*hw1**3/zw+hw1**5/zw+hw1*(1+3*zw**2*hw2**2)/zw-hw2-zw*hw3-hw1**2*(hw2+zw*hw3))/(1+hw1**2)**2.5_dp)*gam
            fluxe=he**3*((2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp)*gam
        endif
        YDOT(i)=(fluxe-fluxw)/(z*dr)/3/mu
        if (i==NEQ) YDOT(i)=-fluxw/(3*mu*Rc*dr+3*mu*Rc*Y(i)/sqrt(3.0_dp))
        if (i==NEQ) YDOT(i)=-fluxw/(3*mu*Rc*dr+3*mu*(Rc-Y(i)/sqrt(3.0_dp))*Y(i)/sqrt(3.0_dp))
!        if (i==NEQ) YDOT(i)=-fluxw/(3*mu*Rc*dr+1.5_dp*mu*(Rc+Y(i)/sqrt(3.0_dp))*Y(i)/sqrt(3.0_dp))
    enddo
!    write(*,*) YDOT
!    stop
end subroutine FEX8
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!simulates film drainage
subroutine  drain
!***************************DECLARATION******************************
    integer :: i,j,fi,fi2
    real(dp) :: h,h1,h2,h3,h4,de
!******************************BODY**********************************
    write(*,*) 'Wellcome to wall drainage.'
    allocate(R(NEQ),Y(NEQ),U(-1:NEQ+2))
    NNZ=NEQ**2 !I really don't know, smaller numbers can make problems
    LENRAT=2 !depends on dp
    allocate(RWORK(int(20+(2+1._dp/LENRAT)*NNZ+(11+9._dp/LENRAT)*NEQ)),IWORK(30))
    ITASK = 1
    ISTATE = 1
    IOPT = 1
    LRW = size(RWORK)
    LIW = size(IWORK)
    IWORK(6)=maxts
    TOUT =T+timestep
    ITOL = 1 !don't change, or you must declare ATOL as ATOL(NEQ)
    dr=Rc/NEQ
    Rs=Rc*sqrt((1+s**2)/s**2)*dstr
    R(1)=dr/2
    do i=2,NEQ
        R(i)=R(i-1)+dr !this is correct only for FDM (slight error for FVM)
    enddo
    open (unit=newunit(fi), file = 'radius.out')
    Y=hi
    j=1
    do i=1,NEQ
!        Y(i)=(Rs+hi)-sqrt(Rs**2-R(i)**2)
        if (i>NEQ*(1-dstr)) then
            Y(i)=(Rs+hi)-sqrt(Rs**2-R(j)**2)
            j=j+1
        endif
        write(fi,"(10000es12.4)") R(i)
    enddo
    close(fi)
    open (unit=newunit(fi), file = 'filmthickness.out')
    open (unit=newunit(fi2), file = 'time.out')
    write(*,'(100A12)') 'film: ','strut: ','total: '
    do i=1,its
        write(fi,"(10000es12.4)") Y(1:NEQ)
        write(fi2,"(10000es12.4)") TOUT
        call DLSODES (sub_ptr, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
        call volume_balance
        TOUT = TOUT+timestep
    enddo
    close(fi)
    close(fi2)
    write(*,*) 'Program exited normally.'
end subroutine drain
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!checks whether we are losing some mass or not
subroutine  volume_balance
!***************************DECLARATION******************************
    integer :: i
    real(dp) :: vf,vs,vt
!******************************BODY**********************************
    vf=0
    do i=1,NEQ
!        vf=vf+2*pi*R(i)*Y(i)*dr
        dr=(Rc-Y(NEQ))/NEQ
        vf=vf+2*pi*dr*(0.5_dp+i-1)*Y(i)*dr
    enddo
    vs=2*pi*Rc*Y(NEQ)**2/sqrt(3.0_dp)
    vs=2*pi*(Rc-Y(NEQ)/sqrt(3.0_dp))*Y(NEQ)**2/sqrt(3.0_dp)
!    vs=pi*(Rc+Y(NEQ)/sqrt(3.0_dp))*Y(NEQ)**2/sqrt(3.0_dp)
    vt=vf+vs
    write(*,'(100es12.3)') vf,vs,vt
end subroutine volume_balance
!********************************************************************
!*******************************END**********************************
end module model
!**********************************************************END*****************************************************************
