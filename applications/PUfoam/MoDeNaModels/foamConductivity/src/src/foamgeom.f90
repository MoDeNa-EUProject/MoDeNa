!> @file
!! subroutines for calculation of geometric properties of the foam
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module foamgeom
    use constants
    use Solve_NonLin
    implicit none
    private
    public foam_morpholgy
contains
!********************************BEGINNING*************************************
!> determine all geometric parameters of the foam
subroutine foam_morpholgy
    integer, parameter :: n=2
    integer :: info,i
    real (dp) :: tol=1e-8_dp
    real (dp), dimension(n) :: x,fvec,diag
    write(*,*) 'Foam morphology:'
    write(mfi,*) 'Foam morphology:'
    select case(morph_input)
    case(1) !dwall is input; fs and dstrut are calculated
        x(1)=fs
        x(2)=dstrut
        call hbrd(fcn_dwall,n,x,fvec,epsilon(pi),tol,info,diag)
        if (info /= 1) then
            write(*,*) 'unable to determine foam morphology parameters, &
                hbrd returned',info
            write(mfi,*) 'unable to determine foam morphology parameters, &
                hbrd returned',info
            stop
        endif
        fs=x(1)
        dstrut=x(2)
        if (fs<0 .or. dstrut< 0) then
            write(*,*) 'unable to determine foam morphology &
                parameters, try different initial guess'
            write(mfi,*) 'unable to determine foam morphology &
                parameters, try different initial guess'
            stop
        endif
    case(2) !fs is input; dwall and dstrut are calculated
        if  (fs<struttol) then
            dwall=(1-por)*dcell/3.775_dp
            dstrut=0
        else
            do i=1,10
                x(1)=dwall*i
                x(2)=dstrut*i
                call hbrd(fcn_fs,n,x,fvec,epsilon(pi),tol,info,diag)
                if (info /= 1) then
                    write(*,*) 'unable to determine foam morphology &
                        parameters, hbrd returned',info, 'restarting'
                    write(mfi,*) 'unable to determine foam morphology &
                        parameters, hbrd returned',info, 'restarting'
                else
                    exit
                endif
            enddo
            if (info /= 1) then
                write(*,*) 'unable to determine foam morphology parameters, &
                    hbrd returned',info
                write(mfi,*) 'unable to determine foam morphology parameters, &
                    hbrd returned',info
                stop
            endif
            dwall=x(1)
            dstrut=x(2)
            if (dwall<0 .or. dstrut< 0) then
                write(*,*) 'unable to determine foam &
                    morphology parameters, try different initial guess'
                write(mfi,*) 'unable to determine foam &
                    morphology parameters, try different initial guess'
                stop
            endif
        endif
    case(3) !dstrut is input; dwall and fs are calculated
        x(1)=dwall
        x(2)=fs
        call hbrd(fcn_dstrut,n,x,fvec,epsilon(pi),tol,info,diag)
        if (info /= 1) then
            write(*,*) 'unable to determine foam morphology parameters, &
                hbrd returned',info
            write(mfi,*) 'unable to determine foam morphology parameters, &
                hbrd returned',info
            stop
        endif
        dwall=x(1)
        fs=x(2)
        if (dwall<0 .or. fs< 0) then
            write(*,*) 'unable to determine foam morphology &
                parameters, try different initial guess'
            write(mfi,*) 'unable to determine foam morphology &
                parameters, try different initial guess'
            stop
        endif
    case(4) !fs is input; dwall and dstrut are calculated
        if  (fs<struttol) then
            dwall=(1-por)*dcell/3.775_dp
            dstrut=0
        else
            x(1)=dwall
            x(2)=dstrut
            call hbrd(fcn_fs2,n,x,fvec,epsilon(pi),tol,info,diag)
            if (info /= 1) then
                write(*,*) 'unable to determine foam morphology parameters, &
                    hbrd returned',info
                write(mfi,*) 'unable to determine foam morphology parameters, &
                    hbrd returned',info
                stop
            endif
            dwall=x(1)
            dstrut=x(2)
            if (dwall<0 .or. dstrut< 0) then
                write(*,*) 'unable to determine foam &
                    morphology parameters, try different initial guess'
                write(mfi,*) 'unable to determine foam &
                    morphology parameters, try different initial guess'
                stop
            endif
        endif
    case default
        write(*,*) 'unknown foam morphology input'
        write(mfi,*) 'unknown foam morphology input'
        stop
    end select
    rhof=(1-por)*rhos
    write(*,'(2x,A,1x,e9.3)') 'porosity:', por
    write(*,'(2x,A,1x,e9.3,1x,A)') 'foam density:', rhof, 'kg/m^3'
    write(*,'(2x,A,1x,e9.3,1x,A)') 'cell size:', dcell*1e6, 'um'
    write(*,'(2x,A,1x,e9.3,1x,A)') 'wall thickness:', dwall*1e6, 'um'
    write(*,'(2x,A,1x,e9.3)') 'strut content:', fs
    write(*,'(2x,A,1x,e9.3,1x,A)') 'strut diameter:', dstrut*1e6, 'um'
    write(mfi,'(2x,A,1x,e9.3)') 'porosity:', por
    write(mfi,'(2x,A,1x,e9.3,1x,A)') 'foam density:', (1-por)*rhos, 'kg/m^3'
    write(mfi,'(2x,A,1x,e9.3,1x,A)') 'cell size:', dcell*1e6, 'um'
    write(mfi,'(2x,A,1x,e9.3,1x,A)') 'wall thickness:', dwall*1e6, 'um'
    write(mfi,'(2x,A,1x,e9.3)') 'strut content:', fs
    write(mfi,'(2x,A,1x,e9.3,1x,A)') 'strut diameter:', dstrut*1e6, 'um'
end subroutine foam_morpholgy
!***********************************END****************************************


!********************************BEGINNING*************************************
!> residual function for fs and dstrut
subroutine fcn_dwall(n,x,fvec,iflag)
    integer, intent(in) :: n
    real (dp), intent(in) :: x(n)
    real (dp), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    real(dp) :: Vcell,Vstruts,Vwalls,fs,dstrut,dcelldd
    fs=x(1)
    dstrut=x(2)
    dcelldd=dcell*(pi/6/0.348_dp)**(1/3._dp)
    Vcell=0.348_dp*dcelldd**3
    Vstruts=2.8_dp*dstrut**2*dcelldd-3.93_dp*dstrut**3
    Vwalls=(1.3143_dp*dcelldd**2-7.367_dp*dstrut*dcelldd+10.323_dp*dstrut**2)*&
        dwall
    fvec(1)=fs-Vstruts/(Vstruts+Vwalls)
    fvec(2)=1-por-(Vstruts+Vwalls)/Vcell
end subroutine fcn_dwall
!***********************************END****************************************


!********************************BEGINNING*************************************
!> residual function for dwall and dstrut
subroutine fcn_fs(n,x,fvec,iflag)
    integer, intent(in) :: n
    real (dp), intent(in) :: x(n)
    real (dp), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    real(dp) :: Vcell,Vstruts,Vwalls,dwall,dstrut,dcelldd
    dwall=x(1)
    dstrut=x(2)
    dcelldd=dcell*(pi/6/0.348_dp)**(1/3._dp)
    Vcell=0.348_dp*dcelldd**3
    Vstruts=2.8_dp*dstrut**2*dcelldd-3.93_dp*dstrut**3
    Vwalls=(1.3143_dp*dcelldd**2-7.367_dp*dstrut*dcelldd+10.323_dp*dstrut**2)*&
        dwall
    fvec(1)=fs-Vstruts/(Vstruts+Vwalls)
    fvec(2)=1-por-(Vstruts+Vwalls)/Vcell
end subroutine fcn_fs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> residual function for dwall and fs
subroutine fcn_dstrut(n,x,fvec,iflag)
    integer, intent(in) :: n
    real (dp), intent(in) :: x(n)
    real (dp), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    real(dp) :: Vcell,Vstruts,Vwalls,dwall,fs,dcelldd
    dwall=x(1)
    fs=x(2)
    dcelldd=dcell*(pi/6/0.348_dp)**(1/3._dp)
    Vcell=0.348_dp*dcelldd**3
    Vstruts=2.8_dp*dstrut**2*dcelldd-3.93_dp*dstrut**3
    Vwalls=(1.3143_dp*dcelldd**2-7.367_dp*dstrut*dcelldd+10.323_dp*dstrut**2)*&
        dwall
    fvec(1)=fs-Vstruts/(Vstruts+Vwalls)
    fvec(2)=1-por-(Vstruts+Vwalls)/Vcell
end subroutine fcn_dstrut
!***********************************END****************************************


!********************************BEGINNING*************************************
!> residual function for dwall and dstrut
!! based on Kaemmerlen (10.1016/j.jqsrt.2009.11.018)
subroutine fcn_fs2(n,x,fvec,iflag)
    integer, intent(in) :: n
    real (dp), intent(in) :: x(n)
    real (dp), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
    real(dp) :: Vcell,Vstruts,Vwalls,dwall,dstrut,dcelldd,x1,x2,x3
    dwall=x(1)
    dstrut=x(2)
    dcelldd=dcell*(pi/6/0.348_dp)**(1/3._dp)
    Vcell=0.348_dp*dcelldd**3
    Vstruts=2.8_dp*dstrut**2*dcelldd
    Vwalls=(1.317_dp*dcelldd**2-13.4284_dp*dstrut*dcelldd+34.2375_dp*dstrut**2)*&
        dwall+(4.639_dp*dcell-17.976_dp*dstrut)*dwall**2
    fvec(1)=fs-Vstruts/(Vstruts+Vwalls)
    fvec(2)=1-por-(Vstruts+Vwalls)/Vcell
end subroutine fcn_fs2
!***********************************END****************************************
end module foamgeom
