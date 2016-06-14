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
    real(dp) :: hi=2.e-6_dp
    real(dp) :: Rc=5.e-5_dp
    real(dp) :: s=1._dp/sqrt(3._dp)
    real(dp) :: Q= 0._dp!2.e-10_dp
    real(dp) :: mu=2._dp! 1.e-1_dp
    real(dp) :: gam=45.0e-3_dp !2.e-2_dp
    real(dp) :: dstr=0.6_dp!0.2_dp
    real(dp) :: vyska,x,ddr,l
    real(dp) :: v=0.5e-6_dp
    real(dp) :: V1,V2, hpriem

    real(dp) :: dr,dr1
    real(dp) :: Rs
    real(dp), dimension(:), allocatable :: R

    real(dp) :: Rb, vt, eps

    !time integration variables for lsode
    integer :: IOUT
    real(dp), dimension(:), allocatable :: RWORK,Y
    integer, dimension(:), allocatable :: IWORK
    real(dp) :: JAC,TOUT,RTOL=1e-12_dp,ATOL=1e-12_dp,T=0
    integer :: IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ=1000, NNZ, LENRAT, MF=222

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
            fluxw=0._dp
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
         if (i==NEQ) then

            YDOT(i)=-fluxw/3/mu/(Rc*dr+l*Y(i)/sqrt(3._dp)-Y(i)**2/4)

            ! write(*,*) YDOT(i),l
            ! stop

                    else
             YDOT(i)=(fluxe-fluxw)/(z*dr)/3/mu
        endif
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
    integer :: i,j,fi,fi1,fi2,fi3,fi4,fi5,fi6,fi7,fi8,fi10,k
    real(dp) :: h,h1,h2,h3,h4,de
!******************************BODY**********************************
    write(*,*) 'Wellcome to wall drainage.'
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
       !Y(i)=(Rs+hi)-sqrt(Rs**2-R(i)**2)
        if (i>NEQ*(1-dstr)) then
          Y(i)=(Rs+hi)-sqrt(Rs**2-R(j)**2)
            j=j+1
        endif
        write(fi,"(10000es12.4)") R(i)
    enddo
    close(fi)
   ! open (unit=newunit(fi), file = 'filmthickness.out')
   ! open (unit=newunit(fi2), file = 'time.out')
    write(*,'(100A12)') 'film: ','strut: ','total: '

    open (unit=50, file = 'lamella')
    open (unit=51, file = 'strut')


    x=Y(NEQ)*s
    l=x+RC

    open (unit=newunit(fi6), file = 'hpriem.out')

    do i=1,its

         if (i==1) then
                open (unit=newunit(fi1), file = 'filmthickness1.out')
                do j=1,NEQ
                    write(fi1,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi1)
            end if

            if (i==100) then
                open (unit=newunit(fi2), file = 'filmthickness2.out')
                do j=1,NEQ
                    write(fi2,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi2)
            end if

            if (i==500) then
                open (unit=newunit(fi3), file = 'filmthickness3.out')
                do j=1,NEQ
                    write(fi3,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi3)
            end if

            if (i==1000) then
                open (unit=newunit(fi4), file = 'filmthickness4.out')
                do j=1,NEQ
                    write(fi4,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi4)
            end if

            if (i==4000) then
                open (unit=newunit(fi5), file = 'filmthickness5.out')
                do j=1,NEQ
                    write(fi5,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi5)
            end if

            if (i==8000) then
                open (unit=newunit(fi10), file = 'filmthickness6.out')
                do j=1,NEQ
                    write(fi10,"(10000es12.4)") R(j),Y(j)
                end do
                close(fi10)
            end if

     !   write(fi,"(10000es12.4)") Y(1:NEQ)
     !   write(fi2,"(10000es12.4)") TOUT

    hpriem=0;
    do k=1, (1-dstr)*NEQ
    hpriem=hpriem+Y(k)
    end do
    hpriem=hpriem/((1-dstr)*NEQ)
    write(fi6,"(10000es12.4)") T, hpriem



        call DLSODES (sub_ptr, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

        call volume_balance
        TOUT = TOUT+timestep


    end do
    close(fi)
    close(fi2)
    close(fi3)
    close(fi4)
    close(fi5)
    close(fi6)
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
        x=Y(NEQ)*s
        Rc=l-x
        dr1=Rc/NEQ
        ddr=(dr1-dr)

        x=Y(NEQ)*s
        Rc=l-x
        dr1=Rc/NEQ

        if (ddr>0) then

        do j=1,NEQ-1
        Y(j)=Y(j)+(Y(j+1)-Y(j))/dr*ddr
        enddo



        elseif (ddr<0) then

      !  write(*,*) Y(NEQ)
        Y(1)=Y(1)-(Y(2)-Y(1))/dr*ddr
        do j=2,NEQ-1
        Y(j)=Y(j)-(Y(j)-Y(j-1))/dr*ddr
        enddo

        dr=dr1

        !write(*,*) Y(NEQ)
        !stop

        end if

        R(1)=dr/2
        do j=2,NEQ
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
!        do j=1,NEQ
!            Y(j)=(Y(j)*pi*((j*dr1)**2-((j-1)*dr1)**2)+(V1-V2)/NEQ)/(pi*((j*dr1)**2-((j-1)*dr1)**2))
!        enddo
!
!        x=Y(NEQ)*s
!        Rc=l-x
!        dr1=Rc/NEQ
!        ! write(*,*) Y(NEQ)
!        ! stop
!        dr=dr1
!        R(1)=dr/2
!        do j=2,NEQ
!            R(j)=R(j-1)+dr
!        enddo

     vs=0
    do i=1,NEQ*(1-dstr)
    vf = vf + pi*(((i)*dr)**2-((i-1)*dr)**2)*Y(i)
    end do
    do i=NEQ*(1-dstr),NEQ
    vs = vs + pi*(((i)*dr)**2-((i-1)*dr)**2)*Y(i)
    end do
    vs=vs+pi*((Rc+x)**2-Rc**2)*Y(NEQ)/2
     vt=vf+vs

     write(50,"(10000es12.4)") T, vf, vf/vt
     write(51,"(10000es12.4)") T, vs, vs/vt
   !  write(*,*) vf,vs,vt
end subroutine volume_balance
!********************************************************************
!*******************************END**********************************
end module model3
!**********************************************************END*****************************************************************

