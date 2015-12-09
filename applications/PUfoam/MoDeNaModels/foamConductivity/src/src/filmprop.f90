!> @file
!! subroutines for evaluation of film radiative properties
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module filmprop
    use constants
    use solidprop
    implicit none
    private
    public scattwall,trextwall,absorwall,scattcoeffintwall,&
        trextcoeffintwall,abscoeffintwall,filmconst
contains
!********************************BEGINNING*************************************
!> determine film reflectance, transmittance and absorbance
subroutine filmconst(lambda,theta,h,Rwin,Twin,Awin)
    real(dp), intent(in) :: &
        lambda,&  !wavelength
        theta,&  !incident angle
        h  !film thickness
    real(dp), intent(out) :: &
        Rwin,&    !reflectance
        Twin,&    !transmittance
        Awin    !absorptance
    integer :: &
        method=1 !1 is recommended, 3 and 4 should work; 2 is for no interference;
        ! 1,3,4 give similar results
    real(dp) :: &
        theta2,&    !refracted agle
        rho12,rho23,&    !reflectivities
        tau12,tau23,&    !transmissivities
        r12,r23,&    !reflection coefficients
        t12,t23,&    !transmission coefficients
        dzeta,&    !interference parameter
        p,q
    complex(dp) :: &
        cn,beta,rc,tc,rc2,tc2
    call optconst(lambda,n2,k2)
    select case(method)
    case(1)
        !according to Modest
        kappa2=4*pi*k2/lambda
    !    theta2=asin(sin(theta)/n2)
    !    rho12=(((cos(theta2)-n2*cos(theta))/(cos(theta2)+n2*cos(theta)))**2+&
            ! ((cos(theta)-n2*cos(theta2))/(cos(theta)+n2*cos(theta2)))**2)/2
    !    rho23=(((n2*cos(theta)-cos(theta2))/(n2*cos(theta)+cos(theta2)))**2+&
            ! ((n2*cos(theta2)-cos(theta))/(n2*cos(theta2)+cos(theta)))**2)/2
    !    write(*,*) rho23
    !    stop
        p=sqrt((sqrt((n2**2-k2**2-n1**2*sin(theta)**2)**2+4*n2**2*k2**2)+&
            (n2**2-k2**2-n1**2*sin(theta)**2))/2)
        q=sqrt((sqrt((n2**2-k2**2-n1**2*sin(theta)**2)**2+4*n2**2*k2**2)-&
            (n2**2-k2**2-n1**2*sin(theta)**2))/2)
        theta2=atan(sin(theta)*n1/p)
        rho12=((n1*cos(theta)-p)**2+q**2)/((n1*cos(theta)+p)**2+q**2)
        rho12=rho12*(1+((p-n1*sin(theta)*tan(theta))**2+q**2)/&
            ((p+n1*sin(theta)*tan(theta))**2+q**2))/2
        rho23=rho12
    !    write(*,*) rho12
    !    stop
        r12=-sqrt(rho12)
        r23=sqrt(rho23)
        tau12=1-rho12
        tau23=1-rho23
        dzeta=4*pi*n2*h/lambda
        Rwin=(r12**2+2*r12*r23*exp(-kappa2*h)*cos(dzeta)+r23**2*exp(-2*kappa2*h))/&
            (1+2*r12*r23*exp(-kappa2*h)*cos(dzeta)+r12**2*r23**2*exp(-2*kappa2*h))
        Twin=(tau12*tau23*exp(-kappa2*h))/&
            (1+2*r12*r23*exp(-kappa2*h)*cos(dzeta)+r12**2*r23**2*exp(-2*kappa2*h))
        Awin=1-Rwin-Twin
    case(2)
        !no interference
        kappa2=4*pi*k2/lambda
        theta2=asin(sin(theta)/n2)
        rho12=(((cos(theta2)-n2*cos(theta))/(cos(theta2)+n2*cos(theta)))**2+&
            ((cos(theta)-n2*cos(theta2))/(cos(theta)+n2*cos(theta2)))**2)/2
        rho23=(((n2*cos(theta)-cos(theta2))/(n2*cos(theta)+cos(theta2)))**2+&
            ((n2*cos(theta2)-cos(theta))/(n2*cos(theta2)+cos(theta)))**2)/2
        tau12=exp(-kappa2*h)
        Rwin=(rho12+(1-2*rho12)*rho23*tau12**2)/(1-rho12*rho23*tau12**2)
        Twin=(1-rho12)*(1-rho23)*tau12/(1-rho12*rho23*tau12**2)
        Awin=1-Rwin-Twin
    case(3)
        !Coquard (according to Ochsner), adapted
        cn=cmplx(n2,-k2,kind=dp)
        theta2=asin(sin(theta)/n2)
        rho12=(((cos(theta2)-n2*cos(theta))/(cos(theta2)+n2*cos(theta)))**2+&
            ((cos(theta)-n2*cos(theta2))/(cos(theta)+n2*cos(theta2)))**2)/2
        rho23=(((n2*cos(theta)-cos(theta2))/(n2*cos(theta)+cos(theta2)))**2+&
            ((n2*cos(theta2)-cos(theta))/(n2*cos(theta2)+cos(theta)))**2)/2
        r12=-sqrt(rho12)
        r23=sqrt(rho23)
        tau12=1-r12**2
        tau23=1-r23**2
        t12=sqrt(tau12/n2)
        t23=sqrt(tau23*n2)
        beta=2*pi*cn*h*cos(theta2)/lambda
        rc=(r12+r23*exp(-iu*2*beta))/(1+r12*r23*exp(-iu*2*beta))
        tc=(t12*t23*exp(-iu*beta))/(1+r12*r23*exp(-iu*2*beta))
        Rwin=abs(rc)**2
        Twin=abs(tc)**2
        Awin=1-Rwin-Twin
    case(4)
        !Coquard (according to Dombrovsky)
        cn=cmplx(n2,-k2,kind=dp)
        theta2=asin(sin(theta)/n2)
        r12=(cos(theta2)-n2*cos(theta))/(cos(theta2)+n2*cos(theta))
    !    r12=(n2*cos(theta)-cos(theta2))/(n2*cos(theta)+cos(theta2)) !this gives
        ! the same results
        r23=-r12
        t12=1-r12**2 !actually this is t12*t23
        beta=2*pi*cn*h*cos(theta2)/lambda
        rc=(r12+r23*exp(-iu*2*beta))/(1+r12*r23*exp(-iu*2*beta))
        tc=(t12*exp(-iu*beta))/(1+r12*r23*exp(-iu*2*beta))
        r12=(cos(theta)-n2*cos(theta2))/(cos(theta)+n2*cos(theta2))
    !    r12=(n2*cos(theta2)-cos(theta))/(n2*cos(theta2)+cos(theta))
        r23=-r12
        t12=1-r12**2 !actually this is t12*t23
        rc2=(r12+r23*exp(-iu*2*beta))/(1+r12*r23*exp(-iu*2*beta))
        tc2=(t12*exp(-iu*beta))/(1+r12*r23*exp(-iu*2*beta))
        Rwin=(abs(rc)**2+abs(rc2)**2)/2
        Twin=(abs(tc)**2+abs(tc2)**2)/2
        Awin=1-Rwin-Twin
    case default
        stop 'unknown method for calculation of slab reflectivity'
    end select
end subroutine filmconst
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates scattering coefficient as a function of wall thickness
real(dp) function scattwall ( dw )
    use quadpack
    real(dp), intent(in) :: dw
    real(dp) :: wtweight,dwold
    !qags variables
    real(dp) :: a !start point of integration
    real(dp) :: abserr
    real(dp) :: b !end point of integration
    real(dp), parameter :: epsabs = 0.0e0_dp
    real(dp), parameter :: epsrel = 1e-3_dp
    integer :: ier
    integer :: neval
    real(dp) :: res
    !end of qags variables
    a=0
    b=pi/2
    wtweight=1/(sqrt(2*pi)*wsdev*dw)*exp(-(log(dw)-log(dwall))**2/(2*wsdev**2))
    dwold=dwall
    dwall=dw
    call qags(scattcoeffintwall, a, b, epsabs, epsrel, res, abserr, neval, ier)
    if (ier /= 0) then
        write(*,*) 'qags returned',ier
        write(*,*) 'wall thickness',dw
        write(*,*) 'scattering coefficient of walls not calculated'
        stop
    endif
    scattwall=res*(1-por)/dwall*wtweight
    dwall=dwold
end function scattwall
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates scattering coefficient as a function of wall thickness
real(dp) function trextwall ( dw )
    use quadpack
    real(dp), intent(in) :: dw
    real(dp) :: wtweight,dwold
    !qags variables
    real(dp) :: a !start point of integration
    real(dp) :: abserr
    real(dp) :: b !end point of integration
    real(dp), parameter :: epsabs = 0.0e0_dp
    real(dp), parameter :: epsrel = 0.001e0_dp
    integer :: ier
    integer :: neval
    real(dp) :: res
    !end of qags variables
    a=0
    b=pi/2
    wtweight=1/(sqrt(2*pi)*wsdev*dw)*exp(-(log(dw)-log(dwall))**2/(2*wsdev**2))
    dwold=dwall
    dwall=dw
    call qags(trextcoeffintwall, a, b, epsabs, epsrel, res, abserr, neval, ier)
    if (ier /= 0) then
        write(*,*) 'qags returned',ier
        write(*,*) 'wall thickness',dw
        write(*,*) 'transport extinction coefficient of walls not calculated'
        stop
    endif
    trextwall=res*(1-por)/dwall*wtweight
    dwall=dwold
end function trextwall
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates absorption coefficient as a function of wall thickness
real(dp) function absorwall ( dw )
    use quadpack
    real(dp), intent(in) :: dw
    real(dp) :: wtweight,dwold
    !qags variables
    real(dp) :: a !start point of integration
    real(dp) :: abserr
    real(dp) :: b !end point of integration
    real(dp), parameter :: epsabs = 0.0e0_dp
    real(dp), parameter :: epsrel = 0.001e0_dp
    integer :: ier
    integer :: neval
    real(dp) :: res
    !end of qags variables
    a=0
    b=pi/2
    wtweight=1/(sqrt(2*pi)*wsdev*dw)*exp(-(log(dw)-log(dwall))**2/(2*wsdev**2))
    dwold=dwall
    dwall=dw
    call qags ( abscoeffintwall, a, b, epsabs, epsrel, res, abserr, neval, ier )
    if (ier /= 0) then
        write(*,*) 'qags returned',ier
        write(*,*) 'wall thickness',dw
        write(*,*) 'absorption coefficient of walls not calculated'
        stop
    endif
    absorwall=res*(1-por)/dwall*wtweight
    dwall=dwold
end function absorwall
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates scattering coefficient of wall integrand
real(dp) function scattcoeffintwall ( theta )
    real(dp), intent(in) :: theta
    real(dp) :: Rwin,Twin,Awin
    call filmconst(lambda,theta,dwall,Rwin,Twin,Awin)
    scattcoeffintwall=Rwin*sin(theta)*cos(theta)
end function scattcoeffintwall
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates transport extinction coefficient of wall integrand
real(dp) function trextcoeffintwall ( theta )
    real(dp), intent(in) :: theta
    real(dp) :: Rwin,Twin,Awin
    call filmconst(lambda,theta,dwall,Rwin,Twin,Awin)
    trextcoeffintwall=(1-Twin+Rwin*cos(2*theta))*sin(theta)*cos(theta)
end function trextcoeffintwall
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates absorption coefficient of wall integrand
real(dp) function abscoeffintwall ( theta )
    real(dp), intent(in) :: theta
    real(dp) :: Rwin,Twin,Awin
    call filmconst(lambda,theta,dwall,Rwin,Twin,Awin)
    abscoeffintwall=Awin*sin(theta)*cos(theta)
end function abscoeffintwall
!***********************************END****************************************
end module filmprop
