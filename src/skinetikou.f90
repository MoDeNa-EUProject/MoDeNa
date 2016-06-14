!******************************************************BEGINNING***************************************************************
!contains model subroutines
!TODO: material balance (to check results)
!TODO: dimensional results
!TODO: complete flow boundary conditions in several models
!TODO: add condition for film breakage when minimum film thickness is below certain threshold
module model3
!*****************************************************DECLARATION**************************************************************
    use constants
    use ioutils
    implicit none

    integer :: maxts=10000
    integer :: its=10000
    real(dp) :: timestep=5.e-3_dp
    real(dp) :: hi=20.e-6_dp
    real(dp) :: Rc=5.e-4_dp
    real(dp) :: s=1._dp/sqrt(3._dp)
    real(dp) :: Q= 0._dp!2.e-10_dp
    real(dp) :: mu=1.e-1_dp
    real(dp) :: gam=45.0e-3_dp !2.e-2_dp
    real(dp) :: dstr=0.3_dp!0.2_dp
    real(dp) :: vyska,x,ddr,l
    real(dp) :: v=5.0e-6_dp
    real(dp) :: V1,V2, hpriem

    real(dp) :: MCO2=44._dp
    real(dp) :: OH0=3400._dp    ! mol/m3 initial concentration of polyol
    real(dp) :: W0=300._dp    !  mol/m3  initial concentration of water
    real(dp) :: NCO0=4000._dp    ! mol/m3  initial concentration of isocyanate

    real(dp) :: AOH=1.7348e0_dp    ! m3/ mol/s frequential factor of gelling reaction
    real(dp) :: EOH=4.04e4_dp    ! J/ mol activation energy of gelling reaction
    real(dp) :: AW=1.385e3_dp    ! 1/s frequential factor of blowing reaction
    real(dp) :: EW=3.266e4_dp    ! J/ mol activation energy of blowing reaction



    real(dp) :: Aeta=4.1e-8_dp    !viscosity constant Aeta
    real(dp) :: Eeta=38.3e3_dp    !viscosity constant Eeta
    real(dp) :: Cg=0.85e0_dp    !viscosity constant Cg
    real(dp) :: AA=4.e0_dp    !viscosity constant AA
    real(dp) :: B=-2.e0_dp    !viscosity constant B

    real(dp) :: rW
    real(dp) :: rNCO

    !entalpicka bilance
    real(dp) :: dHW=-8.6e4_dp !J/ mol
    real(dp) :: dHOH=-7.075e4_dp !J/ mol
    real(dp) :: CP=1800._dp !J/kg K
    real(dp) :: Ccd=836.6_dp !J/kg K
    real(dp) :: rhoP=1100._dp !kg/m3

    integer :: XOHe,XWe,Te

    real(dp) :: dr,dr1
    real(dp) :: Rs
    real(dp), dimension(:), allocatable :: R
    integer :: d=500
    integer :: p=3



    !time integration variables for lsode
    integer :: IOUT
    real(dp), dimension(:), allocatable :: RWORK,Y
    integer, dimension(:), allocatable :: IWORK
    real(dp) :: JAC,TOUT,RTOL=1e-12_dp,ATOL=1e-12_dp,T=0
    integer :: IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ=503, NNZ, LENRAT, MF=222

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
    integer :: NEQ,i,j
    real(dp) :: T, Y(NEQ), YDOT(NEQ)
    real(dp) :: z,ze,zw,zee,zww
    real(dp) :: lame,lamw

    real(dp) :: h,he,hw,hee,hww,heee,hwww
    real(dp) :: he1,hw1,he2,hw2,he3,hw3
    real(dp) :: fluxe,fluxw
!******************************BODY**********************************

    do i=1,d
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
            fluxw=0._dp
        fluxe=he**3*((2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp)*gam
        elseif (i==d) then
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
            fluxe=0._dp
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
            if (i==d-1) then
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
         if (i==d) then

            YDOT(i)=-fluxw/3/mu/(Rc*dr+l*Y(i)/sqrt(3._dp)-Y(i)**2/4)

            ! write(*,*) YDOT(i),l
            ! stop

                    else
             YDOT(i)=(fluxe-fluxw)/(z*dr)/3/mu
        endif
        enddo


        YDOT(xOHe) = AOH*exp(-EOH/Rg/Y(Te))*OH0*(1-Y(xOHe))*(rNCO-2*rW*Y(xWe)-Y(xOHe))
        YDOT(xWe) = AW*exp(-EW/Rg/Y(Te))*(1-Y(xWe))
        YDOT(Te) =((-dHOH)*OH0/rhoP*YDOT(xOHe)+(-dHW)*W0/rhoP*YDOT(xWe))/(CP+Ccd*W0*Y(xWe)*MCO2/1000/rhoP)



!    write(*,*) YDOT
!    stop
end subroutine FEX8
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!simulates film drainage
subroutine  drain
!***************************DECLARATION******************************
    integer :: i,j,fi,fi1,fi2,fi3,fi4,fi5,fi6,fi7,fi8,fi10,k
    real(dp) :: h,h1,h2,h3,h4,de
!******************************BODY**********************************
    write(*,*) 'Wellcome to wall drainage.'
    xOHe=d+1
    XWe=d+2
    Te=d+3
    rW=W0/OH0
    rNCO=NCO0/OH0

    allocate(R(NEQ),Y(NEQ))
    NNZ=NEQ**2 !I really don't know, smaller numbers can make problems
    LENRAT=2!2 !depends on dp
    allocate(RWORK(int(20+(2+1._dp/LENRAT)*NNZ+(11+9._dp/LENRAT)*NEQ)),IWORK(30))
    ITASK = 1
    ISTATE = 1
    IOPT = 1
    LRW = size(RWORK)
    LIW = size(IWORK)
    IWORK(6)=maxts
    TOUT =T+timestep
    ITOL = 1 !don't change, or you must declare ATOL as ATOL(NEQ)
    dr=Rc/(d)
    Rs=Rc*sqrt((1+s**2)/s**2)*dstr
    R(1)=dr/2
    do i=2,(d)
        R(i)=R(i-1)+dr !this is correct only for FDM (slight error for FVM)
    enddo
    open (unit=newunit(fi), file = 'radius.out')
    Y=hi
    j=1
    do i=1,d
       !Y(i)=(Rs+hi)-sqrt(Rs**2-R(i)**2)
        if (i>d*(1-dstr)) then
          Y(i)=(Rs+hi)-sqrt(Rs**2-R(j)**2)
            j=j+1
        endif
        write(fi,"(10000es12.4)") R(i)
    enddo
    close(fi)

    !3.60373E-01 8.67306E-01 3.55076E+02 po 250s
     ! 7.90203E-01  9.99960E-01 4.08735E+02 po 400s
     ! 6.68608E-01 7.39291E-01 9.99651E-01 4.02580E+02 po 375 sek

    Y(xOHe)=7.39291E-01_dp
    Y(xWe)=9.99651E-01_dp
    Y(Te)=4.02580E+02_dp
    mu=Aeta*exp(Eeta/Rg/Y(Te))*(Cg/(Cg-Y(xOHe)))**(AA+B*Y(xOHe))
   ! open (unit=newunit(fi), file = 'filmthickness.out')
   ! open (unit=newunit(fi2), file = 'time.out')
    write(*,'(100A12)') 'film: ','strut: ','total: '

    open (unit=50, file = 'lamella')
    open (unit=51, file = 'strut')
    open (unit=52 ,file = 'reakce')


    x=Y(d)*s
    l=x+RC

    open (unit=newunit(fi6), file = 'hpriem.out')


    do i=1,its

         if (i==1) then
                open (unit=newunit(fi1), file = 'filmthickness1.out')
                do j=1,d
                    write(fi1,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi1)
            end if

            if (i==100) then
                open (unit=newunit(fi2), file = 'filmthickness2.out')
                do j=1,d
                    write(fi2,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi2)
            end if

            if (i==500) then
                open (unit=newunit(fi3), file = 'filmthickness3.out')
                do j=1,d
                    write(fi3,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi3)
            end if

            if (i==1000) then
                open (unit=newunit(fi4), file = 'filmthickness4.out')
                do j=1,d
                    write(fi4,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi4)
            end if

            if (i==4000) then
                open (unit=newunit(fi5), file = 'filmthickness5.out')
                do j=1,d
                    write(fi5,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi5)
            end if

            if (i==8000) then
                open (unit=newunit(fi10), file = 'filmthickness6.out')
                do j=1,d
                    write(fi10,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi10)
            end if

            write(52,"(10000es12.4)") mu,Y(XOHe),Y(XWe),Y(Te)

     !   write(fi,"(10000es12.4)") Y(1:NEQ)
     !   write(fi2,"(10000es12.4)") TOUT

    hpriem=0;
    do k=1, (1-dstr)*d
    hpriem=hpriem+Y(k)
    end do
    hpriem=hpriem/((1-dstr)*d)
    write(fi6,"(10000es12.4)") T, hpriem



        call DLSODES (sub_ptr, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
         mu=Aeta*exp(Eeta/Rg/Y(Te))*(Cg/(Cg-Y(xOHe)))**(AA+B*Y(xOHe))

        call volume_balance
        TOUT = TOUT+timestep


    end do
    close(fi)
    close(fi2)
    close(fi3)
    close(fi4)
    close(fi5)
    close(fi6)
    close(52)
    write(*,*) 'Program exited normally.'
end subroutine drain
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!checks whether we are losing some mass or not
subroutine  volume_balance
!***************************DECLARATION******************************
    integer :: i,j
    real(dp) :: vf,vs,vt
!******************************BODY**********************************
        vf=0
        x=Y(d)*s
        Rc=l-x
        dr1=Rc/(d)
        ddr=(dr1-dr)

        x=Y(d)*s
        Rc=l-x
        dr1=Rc/(d)

        if (ddr>0) then

        do j=1,d-1
        Y(j)=Y(j)+(Y(j+1)-Y(j))/dr*ddr
        enddo



        elseif (ddr<0) then

        write(*,*) Y(d)
        Y(1)=Y(1)-(Y(2)-Y(1))/dr*ddr
        do j=2,d-1
        Y(j)=Y(j)-(Y(j)-Y(j-1))/dr*ddr
        enddo

        dr=dr1

        !write(*,*) Y(NEQ)
        !stop

        end if

        R(1)=dr/2
        do j=2,d
        R(j)=R(j-1)+dr !this is correct only for FDM (slight error for FVM)
        enddo


!        V1=pi*((Rc+x)**2-Rc**2)*Y(NEQ)/2
!        l=l+v*timestep
!        x=Y(NEQ)*s
!        Rc=l-x
!        dr1=Rc/NEQ
!
!        do j=1,NEQ
!            Y(j)=Y(j)*pi*((j*dr)**2-((j-1)*dr)**2)/(pi*((j*dr1)**2-((j-1)*dr1)**2))
!        enddo
!        !write(*,*) Y(NEQ)
!        x=Y(NEQ)*s
!        Rc=l-x
!        dr1=Rc/NEQ
!        V2=pi*((Rc+x)**2-Rc**2)*Y(NEQ)/2
!        do j=1,NEQ-1
!            Y(j)=(Y(j)*pi*((j*dr1)**2-((j-1)*dr1)**2)+(V1-V2)/NEQ)/(pi*((j*dr1)**2-((j-1)*dr1)**2))
!        enddo
         !write(*,*) Y(NEQ)
         !stop
!        dr=dr1
!        R(1)=dr/2
!        do j=2,NEQ
!            R(j)=R(j-1)+dr
!        enddo

     vs=0
    do i=1,d*(1-dstr)
    vf = vf + pi*(((i)*dr)**2-((i-1)*dr)**2)*Y(i)
    end do
    do i=d*(1-dstr),d
    vs = vs + pi*(((i)*dr)**2-((i-1)*dr)**2)*Y(i)
    end do
    vs=vs+pi*((Rc+x)**2-Rc**2)*Y(d)/2
     vt=vf+vs

     write(50,"(10000es12.4)") T, vf, vf/vt
     write(51,"(10000es12.4)") T, vs, vs/vt
     write(*,*) vf,vs,vt
end subroutine volume_balance
!********************************************************************
!*******************************END**********************************
end module model3
!**********************************************************END*****************************************************************


