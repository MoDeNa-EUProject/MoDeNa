!> @file
!! subroutines for evaluation of effective radiative properties of foam
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module foamprop
    use constants
    use ioutils
    use gasprop
    use filmprop
    use cylprop
    implicit none
    private
    public effrad,fbep
    real(dp) :: unin=1.57_dp
    real(dp), dimension(:), allocatable :: &
        lambdaf,&
        !foam radiative properties - wall contribution
        kappafwall,sigmafwall,betafwall,omegafwall,betatrfwall,&
        !foam radiative properties - strut contribution
        kappafstrut,sigmafstrut,betafstrut,omegafstrut,betatrfstrut,&
        !foam radiative properties - gas contribution
        kappafgas,&
        !foam radiative properties - overall
        kappaf,sigmaf,betaf,omegaf,betatrf
contains
!********************************BEGINNING*************************************
!> determine effective radiative properties of foam
subroutine effrad(spectra)
    use quadpack
    character(len=*), intent(in) :: spectra
    integer :: i,j,fi,nwawel=100
!    real(dp) :: theta  !incident angle
!    real(dp) :: Rwin    !reflectance
!    real(dp) :: Twin    !transmittance
!    real(dp) :: Awin    !absorptance
!    real(dp) :: rn    !random number
!    real(dp) :: absp    !absorption parameter
!    real(dp) :: scap    !scattering parameter
!    real(dp) :: npart    !number of particlesper unit volume
    real(dp) :: dwmin,dwmax,lambdamin,lambdamax
    !qags variables
    real(dp) :: a !start point of integration
    real(dp) :: abserr
    real(dp) :: b !end point of integration
    real(dp), parameter :: epsabs = 0.0e0_dp
    real(dp) :: epsrel = 1e-3_dp
    integer :: ier
    integer :: neval
    real(dp) :: res
    !end of qags variables
    write(*,*) 'Radiative properties of the foam:'
    write(mfi,*) 'Radiative properties of the foam:'
!    call cylconst(2*pi,pi/4,1.0e-1_dp,Qs,Qt,mu)
!    stop
    !monte carlo method, does not work properly
!    swin=dcell**2
!    absp=0
!    scap=0
!    npart=3/dcell**3
!    lambda=10e-6
!    call srand(int(abs(sin(dble(time()))*1e6)))
!    do i=1,nrays
!        theta=acos(rand())
!!        theta=rand()*pi/2
!        call filmconst(lambda,theta,dwall,Rwin,Twin,Awin)
!        absp=absp+Swin*Awin
!        scap=scap+Swin*Rwin
!!        absp=absp+cos(theta)*sin(theta)*Swin*Awin*pi/2
!!        scap=scap+cos(theta)*sin(theta)*Swin*Rwin*pi/2
!!        rn=rand()
!!        if (rn>1-Awin) then
!!            absp=absp+Swin
!!        elseif(rn<Rwin) then
!!            scap=scap+Swin
!!        endif
!    enddo
!    kappaf=npart*absp/nrays
!    sigmaf=npart*scap/nrays
!    betaf=kappaf+sigmaf
!    omegaf=sigmaf/betaf
!    write(*,*) kappaf
!    write(*,*) sigmaf
!    write(*,*) betaf
!    write(*,*) omegaf
!    a=0
!    b=pi/2
!    call qags ( reflectance, a, b, epsabs, epsrel, res, abserr, neval, ier )
!    write(*,*)
!    write(*,*) sigmaf
!    sigmaf=3*res/dcell
!    write(*,*) sigmaf

    allocate(lambdaf(nwawel),kappafwall(nwawel),sigmafwall(nwawel),&
        betafwall(nwawel),omegafwall(nwawel),betatrfwall(nwawel),&
        kappafstrut(nwawel),sigmafstrut(nwawel),betafstrut(nwawel),&
        omegafstrut(nwawel),betatrfstrut(nwawel),kappafgas(nwawel),&
        kappaf(nwawel),sigmaf(nwawel),betaf(nwawel),omegaf(nwawel),&
        betatrf(nwawel))
    lambdamin=lambdan(1)
    lambdamax=lambdan(size(lambdan))
    do i=1,nwawel
        lambda=lambdamin+(i-1)*(lambdamax-lambdamin)/nwawel
!        write(*,*) lambda
        lambdaf(i)=lambda
        kappafgas(i)=abscoeffgas(lambda)
        if (fs<1-struttol) then
            epsrel=1e-3_dp
            if (wdist) then
                dwmax=dwall*10 !if you choose larger interval, it can fail for
                ! small standard deviation;
                ! this should work for wsdev between 0.01 and 1
                dwmin=dwall/10
                call qag ( scattwall, dwmin, dwmax, epsabs, epsrel, 1, res, &
                    abserr, neval, ier )
                if (ier /= 0) then
                    write(*,*) 'qag returned',ier
                    write(*,*) 'wavelength',lambda
                    write(*,*) 'scattering coefficient of walls not calculated'
                    write(mfi,*) 'qag returned',ier
                    write(mfi,*) 'wavelength',lambda
                    write(mfi,*) 'scattering coefficient of walls not calculated'
                    stop
                endif
                sigmafwall(i)=res
                call qag ( trextwall, dwmin, dwmax, epsabs, epsrel, 1, res, &
                    abserr, neval, ier )
                if (ier /= 0) then
                    write(*,*) 'qag returned',ier
                    write(*,*) 'wavelength',lambda
                    write(*,*) 'transport extinction coefficient of walls &
                        not calculated'
                    write(mfi,*) 'qag returned',ier
                    write(mfi,*) 'wavelength',lambda
                    write(mfi,*) 'transport extinction coefficient of walls &
                        not calculated'
                    stop
                endif
                betatrfwall(i)=res
                call qag ( absorwall, dwmin, dwmax, epsabs, epsrel, 1, res, &
                    abserr, neval, ier )
                if (ier /= 0) then
                    write(*,*) 'qag returned',ier
                    write(*,*) 'wavelength',lambda
                    write(*,*) 'absorption coefficient of walls not calculated'
                    write(mfi,*) 'qag returned',ier
                    write(mfi,*) 'wavelength',lambda
                    write(mfi,*) 'absorption coefficient of walls not calculated'
                    stop
                endif
                kappafwall(i)=res
            else
                a=0
                b=pi/2
                call qags ( scattcoeffintwall, a, b, epsabs, epsrel, res, &
                    abserr, neval, ier )
                if (ier /= 0) then
                    write(*,*) 'qags returned',ier
                    write(*,*) 'wavelength',lambda
                    write(*,*) 'scattering coefficient of walls not calculated'
                    write(mfi,*) 'qags returned',ier
                    write(mfi,*) 'wavelength',lambda
                    write(mfi,*) 'scattering coefficient of walls not calculated'
                    stop
                endif
                sigmafwall(i)=res*(1-por)/dwall
                call qags ( trextcoeffintwall, a, b, epsabs, epsrel, res, &
                    abserr, neval, ier )
                if (ier /= 0) then
                    write(*,*) 'qags returned',ier
                    write(*,*) 'wavelength',lambda
                    write(*,*) 'transport extinction coefficient of walls &
                        not calculated'
                    write(mfi,*) 'qags returned',ier
                    write(mfi,*) 'wavelength',lambda
                    write(mfi,*) 'transport extinction coefficient of walls &
                        not calculated'
                    stop
                endif
                betatrfwall(i)=res*(1-por)/dwall
                call qags ( abscoeffintwall, a, b, epsabs, epsrel, res, &
                    abserr, neval, ier )
                if (ier /= 0) then
                    write(*,*) 'qags returned',ier
                    write(*,*) 'wavelength',lambda
                    write(*,*) 'absorption coefficient of walls not calculated'
                    write(mfi,*) 'qags returned',ier
                    write(mfi,*) 'wavelength',lambda
                    write(mfi,*) 'absorption coefficient of walls not calculated'
                    stop
                endif
                kappafwall(i)=res*(1-por)/dwall
            endif
            betafwall(i)=kappafwall(i)+sigmafwall(i)
            omegafwall(i)=sigmafwall(i)/betafwall(i)
        else
            kappafwall=0
            sigmafwall=0
            betafwall=0
            omegafwall=0
            betatrfwall=0
        endif
        if (fs>struttol) then
            a=0
            b=pi/2
            if (pi*dstrut/lambda>10) then
                epsrel=1e-1_dp
            elseif(pi*dstrut/lambda>1) then
                epsrel=1e-2_dp
            else
                epsrel=1e-3_dp
            endif
            call qags ( trextcoeffintstrut, a, b, epsabs, epsrel, res, abserr, &
                neval, ier )
            if (ier /= 0) then
                write(*,*) 'qags returned',ier,neval
                write(*,*) 'qags returned',res,abserr
                write(*,*) 'wavelength',lambda
                write(*,*) 'transport extinction coefficient of struts &
                    not calculated'
                stop
            endif
            betatrfstrut(i)=res*4/pi/dstrut*(1-por)
            call qags ( extcoeffintstrut, a, b, epsabs, epsrel, res, abserr, &
                neval, ier )
            if (ier /= 0) then
                write(*,*) 'qags returned',ier
                write(*,*) 'wavelength',lambda
                write(*,*) 'extinction coefficient of struts not calculated'
                stop
            endif
            betafstrut(i)=res*4/pi/dstrut*(1-por)
            call qags ( scattcoeffintstrut, a, b, epsabs, epsrel, res, abserr, &
                neval, ier )
            if (ier /= 0) then
                write(*,*) 'qags returned',ier
                write(*,*) 'wavelength',lambda
                write(*,*) 'scattering coefficient of struts not calculated'
                stop
            endif
            sigmafstrut(i)=res*4/pi/dstrut*(1-por)
            kappafstrut(i)=betafstrut(i)-sigmafstrut(i)
            omegafstrut(i)=sigmafstrut(i)/betafstrut(i)
        else
            kappafstrut=0
            sigmafstrut=0
            betafstrut=0
            omegafstrut=0
            betatrfstrut=0
        endif
    enddo
    kappafwall=(1-fs)*kappafwall
    sigmafwall=(1-fs)*sigmafwall
    betafwall=(1-fs)*betafwall
    betatrfwall=(1-fs)*betatrfwall
    kappafstrut=fs*kappafstrut
    sigmafstrut=fs*sigmafstrut
    betafstrut=fs*betafstrut
    betatrfstrut=fs*betatrfstrut
    kappaf=kappafwall+kappafstrut
    sigmaf=sigmafwall+sigmafstrut
    betaf=betafwall+betafstrut
    betatrf=betatrfwall+betatrfstrut
!    kappaf=(1-fs)*kappafwall+fs*kappafstrut
!    sigmaf=(1-fs)*sigmafwall+fs*sigmafstrut
!    betaf=(1-fs)*betafwall+fs*betafstrut
!    betatrf=(1-fs)*betatrfwall+fs*betatrfstrut
    omegaf=sigmaf/betaf
    open(newunit(fi),file=trim(adjustl(spectra)))
    write(fi,'(1000A23)') '#wavelength','abs.coeff.wall','scatt.coeff.wall',&
        'ext.coeff.wall','albedo wall','tr.ext.coeff.wall','abs.coeff.strut',&
        'scatt.coeff.strut','ext.coeff.strut','albedo strut',&
        'tr.ext.coeff.strut','abs.coeff','scatt.coeff','ext.coeff','albedo ',&
        'tr.ext.coeff'
    do i=1,nwawel
        write(fi,'(1000es23.15)') lambdaf(i),kappafwall(i),sigmafwall(i),&
            betafwall(i),omegafwall(i),betatrfwall(i),kappafstrut(i),&
            sigmafstrut(i),betafstrut(i),omegafstrut(i),betatrfstrut(i),&
            kappaf(i),sigmaf(i),betaf(i),omegaf(i),betatrf(i)
    enddo
    close(fi)
    !calculate gray properties
    a=lambdaf(1)
    b=lambdaf(size(lambdaf))
    b=1e-4_dp !100 um is relatively enough (if you use lower value, you should
    ! incorporate fraction of blackbody radiation function)
    effn=por*n1+(1-por)*unin !this is not perfect, use of some average value
    ! of n2 instead of unin would be better
    call qags ( rosextc, a, b, epsabs, epsrel, res, abserr, neval, ier )
    if (ier /= 0) then
        write(*,*) 'qags returned',ier
        write(*,*) 'rosseland extinction coefficient not calculated'
        stop
    endif
    rossextcoeff=(4*effn**2*sigmab*tmean**3)/res
    call qags ( planckextc, a, b, epsabs, epsrel, res, abserr, neval, ier )
    if (ier /= 0) then
        write(*,*) 'qags returned',ier
        write(*,*) 'planck extinction coefficient not calculated'
        stop
    endif
    planckextcoeff=res/(effn**2*sigmab*tmean**4)
    call qags ( planckalbedo, a, b, epsabs, epsrel, res, abserr, neval, ier )
    if (ier /= 0) then
        write(*,*) 'qags returned',ier
        write(*,*) 'scattering albedo not calculated'
        stop
    endif
    albedo=res/(effn**2*sigmab*tmean**4)
    krad=16*sigmab*tmean**3/(3*rossextcoeff)
    !calculate actual (usually non-gray) properties
    if (nbox<1) then
        stop 'choose number of gray boxes 1 for gray approximation, greater &
            than 1 for non-gray simulation'
    else
        allocate(lambdabox(nbox+1),trextcoeffbox(nbox),albedobox(nbox),&
            abscoeffbox(nbox),scattcoeffbox(nbox))
        if (.not. allocated(fbepbox)) allocate(fbepbox(nbox))
        lambdamin=2e-6_dp
        lambdamax=25e-6_dp
        lambdabox(nbox+1)=100e-6_dp
        lambdabox(nbox)=lambdamax
        lambdabox(1)=lambdamin !if nbox equals 1, then we have only one
        ! box: 2-100 um
        do i=2,nbox-1
            lambdabox(i)=lambdamin+(i-1)*(lambdamax-lambdamin)/(nbox-1._dp)
        enddo
        do i=1,nbox
            fbepbox(i)=fbep(effn,lambdabox(i+1),tmean)-fbep(effn,lambdabox(i)&
                ,tmean)
            call qags(planckextc,lambdabox(i),lambdabox(i+1), epsabs, epsrel, &
                res, abserr, neval, ier)
            if (ier /= 0) then
                write(*,*) 'qags returned',ier
                write(*,*) 'gray box',i
                write(*,*) 'planck extinction coefficient not calculated'
                stop
            endif
            trextcoeffbox(i)=res/(effn**2*sigmab*tmean**4*fbepbox(i))
            call qags(planckalbedo,lambdabox(i),lambdabox(i+1), epsabs, &
                epsrel, res, abserr, neval, ier)
            if (ier /= 0) then
                write(*,*) 'qags returned',ier
                write(*,*) 'gray box',i
                write(*,*) 'scattering albedo not calculated'
                stop
            endif
            albedobox(i)=res/(effn**2*sigmab*tmean**4*fbepbox(i))
        enddo
        scattcoeffbox=albedobox*trextcoeffbox
        abscoeffbox=trextcoeffbox-scattcoeffbox
    endif
    write(*,'(2x,A,1x,es9.3)') 'effective index of refraction:',effn
    write(*,'(2x,A,1x,es9.3,1x,A)') 'Planck mean extinction coefficient:',&
        planckextcoeff,'m^-1'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'Rosseland extinction coefficient:',&
        rossextcoeff,'m^-1'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'Rosseland extinction coefficient:',&
        rossextcoeff/(1-por)/rhos,'kg/m^2'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'radiative conductivity:',&
        krad*1e3_dp,'mW/m/K'
    write(mfi,'(2x,A,1x,es9.3)') 'effective index of refraction:',effn
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'Planck mean extinction coefficient:',&
        planckextcoeff,'m^-1'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'Rosseland extinction coefficient:',&
        rossextcoeff,'m^-1'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'Rosseland extinction coefficient:',&
        rossextcoeff/(1-por)/rhos,'kg/m^2'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'radiative conductivity:',&
        krad*1e3_dp,'mW/m/K'
    deallocate(lambdaf)
    deallocate(kappafwall,sigmafwall,betafwall,omegafwall,betatrfwall,&
        kappafstrut,sigmafstrut,betafstrut,omegafstrut,betatrfstrut,&
        kappaf,sigmaf,betaf,omegaf,betatrf,kappafgas,lambdabox,trextcoeffbox,&
        albedobox)
end subroutine effrad
!***********************************END****************************************



!********************************BEGINNING*************************************
!> evaluate integrand for Rosseland extinction coefficient
real(dp) function rosextc(lambda)
    use interpolation
    real(dp), intent(in) :: lambda  !wavelength
    real(dp) :: beta
    real(dp) :: n  !refractive index
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
!    call optconst(lambda,n,k) !could yield recursive call to qags
    n=effn !used Planck law is valid only constant index of refraction
    if (lambda<lambdaf(1)) then
        write(*,*) 'No data for such low wavelength.'
        write(*,*) 'Minimum wavelength is',lambdaf(1)
        stop
    elseif (lambda>lambdaf(size(lambdaf))) then
!        write(*,*) 'No data for such high wavelength.'
!        write(*,*) 'Maximum wavelength is',lambdaf(size(lambdaf))
!        stop
        beta=betatrf(size(betatrf))
    else
        xi(1)=lambda
        call pwl_interp_1d ( size(lambdaf), lambdaf, betatrf, ni, xi, yi )
        beta=yi(1)
    endif
    rosextc=2*c0**3*exp(hPc*c0/(kb*n*tmean*lambda))*hPc**2*pi/&
        ((exp(hPc*c0/(kb*n*tmean*lambda))-1)**2*kb*n**3*tmean**2*lambda**6)/beta
end function rosextc
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluate integrand for Planck mean extinction coefficient
real(dp) function planckextc(lambda)
    use interpolation
    real(dp), intent(in) :: lambda  !wavelength
    real(dp) :: beta
    real(dp) :: n  !refractive index
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
!    call optconst(lambda,n,k) !could yield recursive call to qags
    n=effn !used Planck law is valid only constant index of refraction
    if (lambda<lambdaf(1)) then
        write(*,*) 'No data for such low wavelength.'
        write(*,*) 'Minimum wavelength is',lambdaf(1)
        stop
    elseif (lambda>lambdaf(size(lambdaf))) then
!        write(*,*) 'No data for such high wavelength.'
!        write(*,*) 'Maximum wavelength is',lambdaf(size(lambdaf))
!        stop
        beta=betatrf(size(betatrf))
    else
        xi(1)=lambda
        call pwl_interp_1d ( size(lambdaf), lambdaf, betatrf, ni, xi, yi )
        beta=yi(1)+abscoeffgas(lambda)
    endif
    planckextc=2*pi*hPc*c0**2/&
        (n**2*lambda**5*(Exp(hPc*c0/(n*lambda*kb*tmean))-1))*beta
end function planckextc
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluate integrand for scattering albedo - Planck style
real(dp) function planckalbedo(lambda)
    use interpolation
    real(dp), intent(in) :: lambda  !wavelength
    real(dp) :: omega
    real(dp) :: n  !refractive index
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
!    call optconst(lambda,n,k) !could yield recursive call to qags
    n=effn !used Planck law is valid only constant index of refraction
    if (lambda<lambdaf(1)) then
        write(*,*) 'No data for such low wavelength.'
        write(*,*) 'Minimum wavelength is',lambdaf(1)
        stop
    elseif (lambda>lambdaf(size(lambdaf))) then
!        write(*,*) 'No data for such high wavelength.'
!        write(*,*) 'Maximum wavelength is',lambdaf(size(lambdaf))
!        stop
        omega=omegaf(size(omegaf))
    else
        xi(1)=lambda
        call pwl_interp_1d ( size(lambdaf), lambdaf, omegaf, ni, xi, yi )
        omega=yi(1)
    endif
    planckalbedo=2*pi*hPc*c0**2/&
        (n**2*lambda**5*(Exp(hPc*c0/(n*lambda*kb*tmean))-1))*omega
end function planckalbedo
!***********************************END****************************************


!********************************BEGINNING*************************************
!> fraction of blackbody radiation according to eq. (1-33) in Siegel's and
!! Howell's 4th Thermal Radiation Heat Transfer
real(dp) function fbep(n,lambda,T)
    integer :: i,maxit=100
    real(dp), intent(in) :: n,lambda,T !index of refraction,wavelength,temperature
    real(dp) :: dzeta,old,res,tol=1e-12_dp
    dzeta=C2/(n*lambda*T)
    old=0
    do i=1,maxit
        fbep=old+15/pi**4*(exp(-i*dzeta)/i*&
            (dzeta**3+3*dzeta**2/i+6*dzeta/i**2+6e0_dp/i**3))
        res=abs(fbep-old)
        if (res<tol) exit
        old=fbep
        if (i==maxit) stop 'fraction of blackbody radiation was &
            not calculated (did not converge)'
    enddo
end function fbep
!***********************************END****************************************
end module foamprop
