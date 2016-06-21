!contains model subroutines
!TODO: add condition for film breakage
module model
    use constants
    use ioutils
    implicit none

    logical :: cylindrical=.true.
    integer :: maxts=5000
    integer :: its=100
    real(dp) :: timestep=1e0_dp
    real(dp) :: hi=5e-6_dp
    real(dp) :: rd=100e-6_dp
    real(dp) :: s=1/sqrt(3._dp)
    real(dp) :: q=0*2e-1_dp
    real(dp) :: mu=1e1_dp
    real(dp) :: gam=25e-3_dp
    real(dp) :: dstr=0.2_dp
    real(dp) :: ndp=4.0_dp
    real(dp) :: mdp=3.0_dp
    real(dp) :: cdp=0.05_dp
    real(dp) :: hdp=1.0e-7_dp
    real(dp) :: bdp=3.0e3_dp*0
    real(dp) :: gr=1e-6

    real(dp) :: dr
    real(dp) :: rs,rc,rc0
    real(dp) :: vsold,vfold
    real(dp), dimension(:), allocatable :: r,u

    !time integration variables for lsode
    integer :: iout
    real(dp), dimension(:), allocatable :: rwork,y
    integer, dimension(:), allocatable :: iwork
    real(dp) :: jac,tout,rtol=1e-8_dp,atol=0,t=0
    integer :: iopt, istate, itask, itol, liw, lrw, neq=200, nnz, lenrat, mf=222

    !needed for selection of model subroutine
    abstract interface
        subroutine sub (neq, t, y, ydot)
            use constants
            integer :: neq
            real(dp) ::  t, y(neq), ydot(neq)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => null ()
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
    real(dp) :: vf,vs,vt
    real(dp) :: dispr,dph
    rc=rc0+gr*t
    rd=rc-y(neq)/sqrt(3.0_dp)
    dr=rd/neq
    call volume_balance(vf,vs,vt)
    q=(vs-vsold)/timestep-(vf-vfold)/timestep
    q=-(vf-vfold)/timestep
    ! q=(vs-vsold)/timestep
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
            call dispress(he,dispr,dph)
            fluxe=gam*he**3*(2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-&
                he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp+&
                dph*he1*he**3*ze
            ! fluxe=gam*he**3*(-he1/ze+he2+he3*ze)
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
            call dispress(he,dispr,dph)
            fluxw=gam*hw**3*(2*hw1**3/zw+hw1**5/zw+hw1*(1+3*zw**2*hw2**2)/zw-&
                hw2-zw*hw3-hw1**2*(hw2+zw*hw3))/(1+hw1**2)**2.5_dp+&
                dph*hw1*hw**3*zw
            ! fluxw=gam*hw**3*(-hw1/zw+hw2+hw3*zw)
            fluxe=-q*3*mu/2/pi !-gr*he*3*mu*zw
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
            call dispress(he,dispr,dph)
            fluxw=gam*hw**3*(2*hw1**3/zw+hw1**5/zw+hw1*(1+3*zw**2*hw2**2)/zw-&
                hw2-zw*hw3-hw1**2*(hw2+zw*hw3))/(1+hw1**2)**2.5_dp+&
                dph*hw1*hw**3*zw
            ! fluxw=gam*hw**3*(-hw1/zw+hw2+hw3*zw)
            call dispress(he,dispr,dph)
            fluxe=gam*he**3*(2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-&
                he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp+&
                dph*he1*he**3*ze
            ! fluxe=gam*he**3*(-he1/ze+he2+he3*ze)
        endif
        ydot(i)=(fluxe-fluxw)/(z*dr)/3/mu
    enddo
end subroutine fex8
!***********************************END****************************************



!********************************BEGINNING*************************************
!fvm, equidistant mesh
!cartesian coordiantes
subroutine  fex1 (neq, t, y, ydot)
    integer :: neq,i
    real(dp) :: t, y(neq), ydot(neq)
    real(dp) :: z,ze,zw,zee,zww
    real(dp) :: lame,lamw
    real(dp) :: h,he,hw,hee,hww,heee,hwww
    real(dp) :: he1,hw1,he2,hw2,he3,hw3
    real(dp) :: fluxe,fluxw
    real(dp) :: vf,vs,vt
    rc=rc0+gr*t
    rd=rc-y(neq)/sqrt(3.0_dp)
    ! dr=rd/neq
    call volume_balance(vf,vs,vt)
    q=(vs-vsold)/timestep
    q=0
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
            fluxe=he**3*he3*gam
            ! fluxe=gam*he**3/(he3/(1+he1**2)**1.5_dp+&
            !     (-he1*he2**2-he3*he1**2)/(1+he1**2)**2.5_dp-&
            !     5*he1**3*he2**2/(1+he1**2)**3.5_dp)
        elseif (i==neq) then
            zww=z
            zw=ze
            z=zee
            lamw=(zw-zww)/(z-zww)
            hww=y(i-1)
            h=y(i)
            hw=h*lamw+hww*(1-lamw)
            he=h+dr*s/2
            ! he=4*h/3-hw/3+s*dr/3
            hwww=hww*lamw+y(i-2)*(1-lamw)
            hw1=(h-hww)/(dr)
            hw2=(h-2*hw+hww)/(dr**2/4)
            hw3=(he-2*h+2*hww-hwww)/(dr**3/4)
            fluxw=hw**3*hw3*gam
            ! fluxw=gam*hw**3/(hw3/(1+hw1**2)**1.5_dp+&
            !     (-hw1*hw2**2-hw3*hw1**2)/(1+hw1**2)**2.5_dp-&
            !     5*hw1**3*hw2**2/(1+hw1**2)**3.5_dp)
            fluxe=-q*3*mu/2/pi! -gr*he*3*mu*zw/4
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
                ! heee=4*hee/3-h/3+s*dr/3
            else
                heee=y(i+2)*lame+hee*(1-lame)
            endif
            hw1=(h-hww)/(dr)
            hw2=(h-2*hw+hww)/(dr**2/4)
            hw3=(he-2*h+2*hww-hwww)/(dr**3/4)
            he1=(hee-h)/(dr)
            he2=(hee-2*he+h)/(dr**2/4)
            he3=(heee-2*hee+2*h-hw)/(dr**3/4)
            fluxw=hw**3*hw3*gam
            ! fluxw=gam*hw**3/(hw3/(1+hw1**2)**1.5_dp+&
            !     (-hw1*hw2**2-hw3*hw1**2)/(1+hw1**2)**2.5_dp-&
            !     5*hw1**3*hw2**2/(1+hw1**2)**3.5_dp)
            fluxe=he**3*he3*gam
            ! fluxe=gam*he**3/(he3/(1+he1**2)**1.5_dp+&
            !     (-he1*he2**2-he3*he1**2)/(1+he1**2)**2.5_dp-&
            !     5*he1**3*he2**2/(1+he1**2)**3.5_dp)
        endif
        ydot(i)=(fluxe-fluxw)/(dr)/3/mu
    enddo
end subroutine fex1
!***********************************END****************************************


!********************************BEGINNING*************************************
!simulates film drainage
subroutine  drain
    integer :: i,fi,fi2
    real(dp) :: h,h1,h2,h3,h4,de
    real(dp) :: vf,vs,vt
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
    if (cylindrical) then
        sub_ptr => fex8
    else
        sub_ptr => fex1
    endif
    dr=rd/neq
    rs=rd*sqrt((1+s**2)/s**2)*dstr
    do i=1,neq
        r(i)=dr*(0.5_dp+i-1)
    enddo
    y=hi
    do i=1,neq
        if (i>neq*(1-dstr)) then
            y(i)=(rs+hi)-sqrt(rs**2-(r(i)-(1-dstr)*rd)**2)
        endif
        ! y(i)=s/Rd/2*r(i)**2+hi
        ! y(i)=s/Rd**2/3*r(i)**3+hi
    enddo
    rc0=rd+y(neq)/sqrt(3.0_dp)
    rc=rc0
    open (unit=newunit(fi), file = 'filmthickness.csv')
    open (unit=newunit(fi2), file = 'results_1d.csv')
    call volume_balance(vf,vs,vt)
    write(*,'(1x,100a12)') 'time:','dr:','film: ','strut: ','total: '
    write(*,'(100es12.3)') t,dr,vf,vs,vt
    write(fi,"(10000es12.4)") y(1:neq)
    write(unit=fi2, fmt='(10000a12)') '#time','dr','np','vf','vs','vt'
    write(unit=fi2, fmt='(10000es12.4)') t,dr,dble(neq),vf,vs,vt
    vsold=vs
    vfold=vf
    do i=1,its
        write(fi,"(10000es12.4)") y(1:neq)
        call dlsodes (sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
            istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        rc=rc0+gr*t
        rd=rc-y(neq)/sqrt(3.0_dp)
        dr=rd/neq
        call volume_balance(vf,vs,vt)
        write(*,'(100es12.3)') t,dr,vf,vs,vt
        write(unit=fi2, fmt='(10000es12.4)') t,dr,dble(neq),vf,vs,vt
        vsold=vs
        vfold=vf
        tout = tout+timestep
    enddo
    close(fi)
    close(fi2)
    write(*,*) 'program exited normally.'
end subroutine drain
!***********************************END****************************************


!********************************BEGINNING*************************************
!checks whether we are losing some mass or not
pure subroutine  volume_balance(vf,vs,vt)
    integer :: i
    real(dp), intent(out) :: vf,vs,vt
    vf=0
    if (cylindrical) then
        do i=1,neq
            vf=vf+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
        enddo
        vs=pi*y(neq)/3*(rd**2+rc*rd+rc**2)-pi*y(neq)*rd**2
    else
        do i=1,neq
            vf=vf+y(i)*dr
        enddo
        vs=y(neq)**2/2/sqrt(3.0_dp)
    endif
    vt=vf+vs
end subroutine volume_balance
!***********************************END****************************************


!********************************BEGINNING*************************************
!disjoining pressure
subroutine dispress(h,dispr,dph)
    real(dp), intent(in) :: h
    real(dp), intent(out) :: dispr !disjoining pressure
    real(dp), intent(out) :: dph !derivative of disjoining pressure
    dispr=bdp*((hdp/h)**(ndp-1)-(hdp/h)**(mdp-1))*(hdp-cdp)
    dph=(bdp*(hdp/h)**mdp*(hdp*mdp+cdp*(h-h*mdp))-&
        bdp*(hdp/h)**ndp*(hdp*ndp+cdp*(h-h*ndp)))/(h*hdp)
end subroutine dispress
!***********************************END****************************************
end module model
