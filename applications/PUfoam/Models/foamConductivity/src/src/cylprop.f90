!> @file
!! subroutines for evaluation of strut radiative properties
!! struts are modeled as cylinders
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module cylprop
    use constants
    use solidprop
    use specfun
    implicit none
    private
    public trextcoeffintstrut,extcoeffintstrut,scattcoeffintstrut,cylconst
contains
!********************************BEGINNING*************************************
!> determine efficiency factors of scattering, extinction and asymmetry factor
!!of infinitely long cylinder
!!according to Dombrovsky: Thermal Radiation in Disperse Systems: An Engineering
!!Approach (2010)
subroutine cylconst(lambda,theta,radius,Qs,Qt,gcyl)
    real(dp), intent(in) :: &
        lambda,&  !wavelength
        theta,&  !incident angle
        radius  !radius of cylinder
    real(dp), intent(out) :: &
        Qs,&    !efficiency factor of scattering
        Qt,&    !efficiency factor of extinction
        gcyl    !anisotropy factor
    integer :: &
        i,k,kmax,nwaw=1001
    real(dp) :: &
        x,&   !diffraction parameter
        tol=1e-10_dp,&   !tolerance for calculation of Bessel functions
        phi,& !scattered angle
        s1,s2,Qte,Qth,Qse,Qsh!,mue,muh
    complex(dp) :: &
        m,&    !complex index of refraction
        u,v,s,Dk,Rk,Ak1,Akm,Bk1,Bkm,Deltak
    complex(dp) :: &
        T1,T2,T3,T4
    real(dp), dimension(:), allocatable :: &
        adsr !angular distribution of scattered radiation
    complex(dp), dimension(:), allocatable :: &
        ae,ah,be,bh,ju,dju,yu,dyu,jv,djv,yv,dyv,h1v,dh1v,h2v,dh2v
    call optconst(lambda,n2,k2)
!    n2=1.6_dp
!    k2=0.01_dp
    m=cmplx(n2,-k2,kind=dp)
    x=2*pi*radius/lambda
    u=x*sqrt(m**2-sin(theta)**2)
    v=x*cos(theta)
    s=1/u**2-1/v**2
    kmax=nint(x+4*x**(1/3.0_dp)+2)+2 !+2 for safety (also there is most likely
    ! a mistake in Mackowski code on this place)
    if (kmax>50) kmax=50
!    write(*,*) kmax
    allocate(ae(0:kmax),ah(0:kmax),be(0:kmax),bh(0:kmax),ju(0:kmax),dju(0:kmax),&
        yu(0:kmax),dyu(0:kmax),jv(0:kmax),djv(0:kmax),yv(0:kmax),dyv(0:kmax),&
        h1v(0:kmax),dh1v(0:kmax),h2v(0:kmax),dh2v(0:kmax))
    call bessel(kmax,u,tol,ju,dju,yu,dyu)
    call bessel(kmax,v,tol,jv,djv,yv,dyv)
    call hankel(kmax,v,tol,h1v,dh1v,h2v,dh2v)
    do k=0,kmax
        Dk=k*s*sin(theta)
        Rk=jv(k)/h2v(k)
        Ak1=dh2v(k)/(v*h2v(k))-dju(k)/(u*ju(k))
        Akm=dh2v(k)/(v*h2v(k))-m**2*dju(k)/(u*ju(k))
        Bk1=djv(k)/(v*jv(k))-dju(k)/(u*ju(k))
        Bkm=djv(k)/(v*jv(k))-m**2*dju(k)/(u*ju(k))
        Deltak=Akm*Ak1-Dk**2
        ae(k)=iu*Dk*Rk*(Bk1-Ak1)/Deltak
        ah(k)=Rk*(Akm*Bk1-Dk**2)/Deltak
        be(k)=Rk*(Ak1*Bkm-Dk**2)/Deltak
        bh(k)=-ae(k)
    enddo
    Qte=real(0.5_dp*be(0),kind=dp)
    Qth=real(0.5_dp*ah(0),kind=dp)
    Qse=real(0.5_dp*be(0)*conjg(be(0)),kind=dp)
    Qsh=real(0.5_dp*ah(0)*conjg(ah(0)),kind=dp)
    do k=1,kmax
        Qte=Qte+real(be(k),kind=dp)
        Qth=Qth+real(ah(k),kind=dp)
        Qse=Qse+real(ae(k)*conjg(ae(k))+be(k)*conjg(be(k)),kind=dp)
        Qsh=Qsh+real(ah(k)*conjg(ah(k))+bh(k)*conjg(bh(k)),kind=dp)
    enddo
    Qte=4/x*Qte
    Qth=4/x*Qth
    Qse=4/x*Qse
    Qsh=4/x*Qsh
    Qs=(Qse+Qsh)/2
    Qt=(Qte+Qth)/2
!    write(*,*) Qte,Qse,Qte-Qse
!    write(*,*) Qth,Qsh,Qth-Qsh
    !see Kaemmerlen, 2010: 10.1016/j.jqsrt.2009.11.018
!    Qs=Qs*(0.0015_dp*(lambda/(2*radius))**3-0.0219_dp*(lambda/(2*radius))**2+&
        ! 0.0688_dp*(lambda/(2*radius))+0.7639_dp)
!    Qt=Qt*(0.0015_dp*(lambda/(2*radius))**3-0.0219_dp*(lambda/(2*radius))**2+&
        ! 0.0688_dp*(lambda/(2*radius))+0.7639_dp)

    ! asymmetry factor is not needed, we need gcyl according to Placido,
    ! calculated as Dombrovsky
!    mue=Qse*sin(theta)**2*x/(4*cos(theta)**2)
!    muh=Qsh*sin(theta)**2*x/(4*cos(theta)**2)
!    do k=0,kmax-1
!        mue=mue+real(ae(k)*conjg(ae(k+1))+be(k)*conjg(be(k+1)),kind=dp)
!        muh=muh+real(ah(k)*conjg(ah(k+1))+bh(k)*conjg(bh(k+1)),kind=dp)
!    enddo
!    mue=mue*4*cos(theta)**2/x
!    muh=muh*4*cos(theta)**2/x
!    mu=(mue+muh)/(Qse+Qsh) !mue is actually mue*Qse

    allocate(adsr(nwaw))
    do i=1,nwaw
        phi=(i-1)*pi/(nwaw)
        T1=be(0)
        T2=ah(0)
        T3=0
        T4=0
        do k=1,kmax
            T1=T1+2*be(k)*cos(k*phi)
            T2=T2+2*ah(k)*cos(k*phi)
            T3=T3-2*iu*ae(k)*sin(k*phi)
            T4=-T3
        enddo
        adsr(i)=(abs(T1)**2+abs(T2)**2)/(pi*x*Qs)
        adsr(i)=(abs(T1)**2+abs(T2)**2+abs(T3)**2+abs(T4)**2)/2
    enddo
    s1=adsr(1)*cos(0e0_dp)+adsr(nwaw)*cos(pi)
    s2=adsr(1)+adsr(nwaw)
    do i=2,nwaw-1
        phi=(i-1)*pi/(nwaw-1)
        s1=s1+2*adsr(i)*cos(phi)
        s2=s2+2*adsr(i)
    enddo
    s1=s1*pi/(2*(nwaw-1))
    s2=s2*pi/(2*(nwaw-1))
    gcyl=s1/s2
!    write(*,*) theta,s2*lambda/pi**2/radius,Qs
!    write(*,*) s1/s2,mu
!    phi=pi/2
!    T1=be(0)
!    T2=ah(0)
!    T3=0
!    T4=0
!    do k=1,kmax
!        T1=T1+2*be(k)*cos(k*phi)
!        T2=T2+2*ah(k)*cos(k*phi)
!        T3=T3-2*iu*ae(k)*sin(k*phi)
!        T4=-T3
!    enddo
!    write(*,*) abs(T1)**2/Qs,abs(T2)**2/Qs,abs(T3)**2/Qs

!    Qs=abs(be(0))**2+abs(ah(0))**2 !according to Modest, doesn't work
!    do i=1,kmax
!        Qs=Qs+abs(be(i))**2+abs(ah(i))**2+abs(bh(i))**2+abs(ae(i))**2
!    enddo
!    Qs=Qs/x
!    write(*,*) Qs
!    write(*,*) abs(T1)**2/Qs
    deallocate(adsr)
end subroutine cylconst
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates transport extinction coefficient of strut integrand
real(dp) function trextcoeffintstrut ( theta )
    real(dp), intent(in) :: theta
    real(dp) :: Qs,Qt,gcyl
    call cylconst(lambda,theta,dstrut/2,Qs,Qt,gcyl)
    trextcoeffintstrut=(Qt-Qs*sin(theta)**2-gcyl*Qs*cos(theta)**2)*cos(theta)
end function trextcoeffintstrut
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates extinction coefficient of strut integrand
real(dp) function extcoeffintstrut ( theta )
    real(dp), intent(in) :: theta
    real(dp) :: Qs,Qt,gcyl
    call cylconst(lambda,theta,dstrut/2,Qs,Qt,gcyl)
    extcoeffintstrut=Qt*cos(theta)
end function extcoeffintstrut
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates scattering coefficient of strut integrand
real(dp) function scattcoeffintstrut ( theta )
    real(dp), intent(in) :: theta
    real(dp) :: Qs,Qt,gcyl
    call cylconst(lambda,theta,dstrut/2,Qs,Qt,gcyl)
    scattcoeffintstrut=Qs*cos(theta)
end function scattcoeffintstrut
!***********************************END****************************************
end module cylprop
