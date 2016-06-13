!contains model subroutines
!TODO: material balance (to check results)
!TODO: dimensional results
!TODO: complete flow boundary conditions in several models
!TODO: add condition for film breakage
module model
    use constants
    use ioutils
    implicit none

    integer :: maxts=5000
    integer :: its=100
    real(dp) :: timestep=1e-3_dp
    real(dp) :: hi=5e-6_dp
    real(dp) :: rc=1e-4_dp
    real(dp) :: s=1/sqrt(3._dp)
    real(dp) :: q=0*2e-1_dp
    real(dp) :: mu=1e-1_dp
    real(dp) :: gam=25e-3_dp
    real(dp) :: dstr=0.2_dp

    real(dp) :: dr
    real(dp) :: rs
    real(dp), dimension(:), allocatable :: r,u

    !time integration variables for lsode
    integer :: iout
    real(dp), dimension(:), allocatable :: rwork,y
    integer, dimension(:), allocatable :: iwork
    real(dp) :: jac,tout,rtol=1e-8_dp,atol=0,t=0
    integer :: iopt, istate, itask, itol, liw, lrw, neq=400, nnz, lenrat, mf=222

    !needed for selection of model subroutine
    abstract interface
        subroutine sub (neq, t, y, ydot)
            use constants
            integer :: neq
            real(dp) ::  t, y(neq), ydot(neq)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => fex8!null ()
contains
!********************************BEGINNING*************************************
!fvm, equidistant mesh
!cylindrical geometry, sin(alpha)=(dh/dr)/(1+(dh/dr)**2)
subroutine  fex8 (neq, t, y, ydot)
    integer :: neq,i
    real(dp) :: t, y(neq), ydot(neq)
    real(dp) :: z,ze,zw,zee,zww
    real(dp) :: lame,lamw
    real(dp) :: h,he,hw,hee,hww,heee,hwww
    real(dp) :: he1,hw1,he2,hw2,he3,hw3
    real(dp) :: fluxe,fluxw
    dr=rc/neq
    ! dr=(rc-y(neq))/neq
    do i=1,neq
        if (i==1) then
            z=dr/2
            ze=dr
            zee=3*dr/2
            lame=(ze-z)/(zee-z)
            h=y(i)
            hee=y(i+1)
            hw=h
            he=hee*lame+h*(1-lame)
            heee=y(i+2)*lame+hee*(1-lame)
            he1=(hee-h)/(dr)
            he2=(hee-2*he+h)/(dr**2/4)
            he3=(heee-2*hee+2*h-hw)/(dr**3/4)
            fluxw=0
            fluxe=he**3*((2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-&
                he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp)*gam
        elseif (i==neq) then
            zww=z
            zw=ze
            z=zee
            lamw=(zw-zww)/(z-zww)
            hww=y(i-1)
            h=y(i)
            hw=h*lamw+hww*(1-lamw)
            he=h+dr*s/2
            hwww=hww*lamw+y(i-2)*(1-lamw)
            hw1=(h-hww)/(dr)
            hw2=(h-2*hw+hww)/(dr**2/4)
            hw3=(he-2*h+2*hww-hwww)/(dr**3/4)
            fluxw=hw**3*((2*hw1**3/zw+hw1**5/zw+hw1*(1+3*zw**2*hw2**2)/zw-&
                hw2-zw*hw3-hw1**2*(hw2+zw*hw3))/(1+hw1**2)**2.5_dp)*gam
            fluxe=-q
        else
            zww=z
            zw=ze
            z=zee
            ze=ze+dr
            zee=ze+dr/2
            lamw=(zw-zww)/(z-zww)
            lame=(ze-z)/(zee-z)
            hww=y(i-1)
            h=y(i)
            hee=y(i+1)
            hw=h*lamw+hww*(1-lamw)
            he=hee*lame+h*(1-lame)
            if (i==2) then
                hwww=hww
            else
                hwww=hww*lamw+y(i-2)*(1-lamw)
            endif
            if (i==neq-1) then
                heee=hee+dr*s/2
            else
                heee=y(i+2)*lame+hee*(1-lame)
            endif
            hw1=(h-hww)/(dr)
            hw2=(h-2*hw+hww)/(dr**2/4)
            hw3=(he-2*h+2*hww-hwww)/(dr**3/4)
            he1=(hee-h)/(dr)
            he2=(hee-2*he+h)/(dr**2/4)
            he3=(heee-2*hee+2*h-hw)/(dr**3/4)
            fluxw=hw**3*((2*hw1**3/zw+hw1**5/zw+hw1*(1+3*zw**2*hw2**2)/zw-&
                hw2-zw*hw3-hw1**2*(hw2+zw*hw3))/(1+hw1**2)**2.5_dp)*gam
            fluxe=he**3*((2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-&
                he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp)*gam
        endif
        ydot(i)=(fluxe-fluxw)/(z*dr)/3/mu
        ! if (i==neq) ydot(i)=-fluxw/(3*mu*rc*dr+3*mu*rc*y(i)/sqrt(3.0_dp))
        ! if (i==neq) ydot(i)=-fluxw/(3*mu*rc*dr+&
        !     3*mu*(rc-y(i)/sqrt(3.0_dp))*y(i)/sqrt(3.0_dp))
        ! if (i==neq) ydot(i)=-fluxw/(3*mu*rc*dr+&
        !     1.5_dp*mu*(rc+y(i)/sqrt(3.0_dp))*y(i)/sqrt(3.0_dp))
    enddo
end subroutine fex8
!***********************************END****************************************


!********************************BEGINNING*************************************
!simulates film drainage
subroutine  drain
    integer :: i,j,fi,fi2
    real(dp) :: h,h1,h2,h3,h4,de
    write(*,*) 'wellcome to wall drainage.'
    allocate(r(neq),y(neq),u(-1:neq+2))
    nnz=neq**2 !i really don't know, smaller numbers can make problems
    lenrat=2 !depends on dp
    allocate(rwork(int(20+(2+1._dp/lenrat)*nnz+(11+9._dp/lenrat)*neq)),&
        iwork(30))
    itask = 1
    istate = 1
    iopt = 1
    rwork(5:10)=0
    iwork(5:10)=0
    lrw = size(rwork)
    liw = size(iwork)
    iwork(6)=maxts
    tout =t+timestep
    itol = 1 !don't change, or you must declare atol as atol(neq)
    dr=rc/neq
    rs=rc*sqrt((1+s**2)/s**2)*dstr
    r(1)=dr/2
    do i=2,neq
        r(i)=r(i-1)+dr
    enddo
    open (unit=newunit(fi), file = 'radius.out')
    y=hi
    j=1
    do i=1,neq
        if (i>neq*(1-dstr)) then
            y(i)=(rs+hi)-sqrt(rs**2-r(j)**2)
            j=j+1
        endif
        write(fi,"(10000es12.4)") r(i)
    enddo
    close(fi)
    open (unit=newunit(fi), file = 'filmthickness.out')
    open (unit=newunit(fi2), file = 'time.out')
    write(*,'(1x,100a12)') 'film: ','strut: ','total: '
    call volume_balance
    do i=1,its
        write(fi,"(10000es12.4)") y(1:neq)
        write(fi2,"(10000es12.4)") tout
        call dlsodes (sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        call volume_balance
        tout = tout+timestep
    enddo
    close(fi)
    close(fi2)
    write(*,*) 'program exited normally.'
end subroutine drain
!***********************************END****************************************


!********************************BEGINNING*************************************
!checks whether we are losing some mass or not
subroutine  volume_balance
    integer :: i
    real(dp) :: vf,vs,vt
    vf=0
    do i=1,neq
        ! dr=(rc-y(neq))/neq
        vf=vf+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
    enddo
    vs=2*pi*rc*y(neq)**2/sqrt(3.0_dp)
    ! vs=2*pi*(rc-y(neq)/sqrt(3.0_dp))*y(neq)**2/sqrt(3.0_dp)
    ! vs=pi*(rc+y(neq)/sqrt(3.0_dp))*y(neq)**2/sqrt(3.0_dp)
    vs=pi*y(neq)/3*((rc+y(neq)/sqrt(3.0_dp))**2+(rc+y(neq)/sqrt(3.0_dp))*rc+&
        rc**2)-pi*y(neq)*rc**2
    vt=vf+vs
    write(*,'(100es12.3)') vf,vs,vt
end subroutine volume_balance
!***********************************END****************************************
end module model
