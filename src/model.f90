!contains main model
module model
    use constants, only: dp
    implicit none
    private
    real(dp), parameter :: &
        s=1/sqrt(3._dp) !film thickness derivative at outer domain boundary
    real(dp) :: q !flux into domain from strut
    public odesystem,set_initial_conditions,q
contains
!********************************BEGINNING*************************************
!fvm, equidistant mesh
!cylindrical geometry, sin(alpha)=(dh/dr)/(1+(dh/dr)**2)
subroutine odesystem(neq, t, y, ydot)
    use constants, only: pi
    use globals
    use phys_prop, only: dispress,Rb_der,visc
    integer :: neq,i
    real(dp) :: t, y(neq), ydot(neq)
    real(dp) :: z,ze,zw,zee,zww
    real(dp) :: lame,lamw
    real(dp) :: h,he,hw,hee,hww,heee,hwww
    real(dp) :: he1,hw1,he2,hw2,he3,hw3
    real(dp) :: fluxe,fluxw
    real(dp) :: vf,vs,vt
    real(dp) :: dispr,dph
    mu=visc(t)
    rc=rc0+gr*t
    rd=rc-y(neq)/sqrt(3.0_dp)
    dr=rd/neq
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
            fluxe=-q*3*mu/2/pi
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
            call dispress(he,dispr,dph)
            fluxe=gam*he**3*(2*he1**3/ze+he1**5/ze+he1*(1+3*ze**2*he2**2)/ze-&
                he2-ze*he3-he1**2*(he2+ze*he3))/(1+he1**2)**2.5_dp+&
                dph*he1*he**3*ze
        endif
        ydot(i)=(fluxe-fluxw)/(z*dr)/3/mu
    enddo
end subroutine odesystem
!***********************************END****************************************


!********************************BEGINNING*************************************
!disjoining pressure
subroutine set_initial_conditions(y)
    use globals
    use in_out, only: load_bubble_growth
    use phys_prop, only: bblgr_res,Rb_spline_ini,visc_spline_ini
    real(dp), dimension(:), intent(out) :: y
    integer :: i,neq
    real(dp) :: tmp
    real(dp), dimension(:), allocatable :: r
    neq=size(y)
    dr=rd/neq
    rs=rd*sqrt((1+s**2)/s**2)*dstr
    allocate(r(neq))
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
    call load_bubble_growth(bblgr_res)
    call Rb_spline_ini
    call visc_spline_ini
end subroutine set_initial_conditions
!***********************************END****************************************
end module model
