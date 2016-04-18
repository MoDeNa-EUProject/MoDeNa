!> @file
!! contains model subroutines for bubble growth model
!! non-dimensional subroutines
!! not updated to latest version, still needs work to be used
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module model_nd
    use iso_c_binding
    use constants
    use globals
    use foaming_globals_m
    use fmodena
    use modenastuff
    use model
    use phys_prop, only:physical_properties,Rb,Rderiv
    implicit none
    private
    !time integration variables for lsode
    integer :: iout, iopt, istate, itask, itol, liw, lrw, nnz, lenrat, neq, mf
    real(dp) :: jac,tout,rtol,atol,t
    real(dp), dimension(:), allocatable :: rwork!,y
    integer, dimension(:), allocatable :: iwork
    !mesh variables
    integer :: info
    real(dp),allocatable :: atri(:),btri(:),ctri(:),rtri(:),utri(:)
    ! interpolation variables
    logical :: Rb_initialized
    integer :: Rb_kx=5,Rb_iknot=0,Rb_inbvx
    real(dp), dimension(:), allocatable :: Rb_tx,Rb_coef
    !needed for selection of subroutine for evaluation of derivatives
    abstract interface
        subroutine sub (neq, t, y, ydot)
            use constants
            integer :: neq
            real(dp) ::  t, y(neq), ydot(neq)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => odesystem_nd
    public odesystem_nd
contains
!********************************BEGINNING*************************************
!> model supplied to integrator, FVM, nonequidistant mesh
subroutine  odesystem_nd (neq, t, y, ydot)
    integer :: neq,i,j
    real(dp) :: t,y(neq),ydot(neq),z,zw,ze,zww,zee,lamw,lame,cw,ce,cww,cee,&
        c,dcw,dce,dil,bll
    call dim_var_nd
    call molar_balance(y)
    call nondim_var
    ydot=0
    ydot(xOHeq) = AOH*exp(-EOH/Rg/y(teq)/temp0)*(1-y(xOHeq))*&
        (NCO0-2*W0*y(xWeq)-OH0*y(xOHeq))*R0**2/D0(1) !polyol conversion
    if (kin_model==3) then
        if (y(xOHeq)>0.5_dp .and. y(xOHeq)<0.87_dp) ydot(xOHeq)=ydot(xOHeq)*&
            (-2.027_dp*y(xOHeq)+2.013_dp) !gelling influence on kinetics
        if (y(xOHeq)>0.87_dp) ydot(xOHeq)=ydot(xOHeq)*&
            (3.461_dp*y(xOHeq)-2.761_dp)
    endif
    if (W0>1e-3) then
        ! water conversion
        ! ydot(xWeq) = AW*exp(-EW/Rg/y(teq)/temp0)*(1-y(xWeq))*&
        !     (NCO0-2*W0*y(xWeq)-OH0*y(xOHeq))*R0**2/D0 2nd order
        ydot(xWeq) = AW*exp(-EW/Rg/y(teq)/temp0)*(1-y(xWeq))*R0**2/D0(1) !1st order
    endif
    if (dilution) then
        if (co2_pos==1) then
            bll=mb(2)/Vsh*Mbl(2)/rhop
        else
            bll=mb(1)/Vsh*Mbl(1)/rhop
        endif
        dil=1/(1+rhop/rhobl*bll)
        ydot(xOHeq)=ydot(xOHeq)*dil
        ydot(xWeq)=ydot(xWeq)*dil
    endif
    if (kin_model==4) then
        call kinModel(y)
        do i=1,size(kineq)
            ydot(kineq(i))=kinsource(i)
        enddo
    endif
    !temperature (enthalpy balance)
    ydot(teq) = -dHOH*OH0/(rhop*cp/temp0)*ydot(xOHeq)&
                -dHW*W0/(rhop*cp/temp0)*ydot(xWeq)
    do i=1,ngas
        ydot(teq) = ydot(teq) - dHv(i)*12*pi*Mbl(i)*D(i)*radius**4*R0**4/&
            (rhop*cp*Vsh)*(y(fceq+i-1)-KH(i)*y(fpeq+i-1))/(dz(1)/2)*R0**2/D0(1)
    enddo
    if (kin_model==4) then
        ! ydot(kineq(19))=ydot(teq)
    endif
    if (firstrun) then
        if (inertial_term) then
            ydot(req) = y(req+1)    !radius (momentum balance)
            ydot(req+1) = sum(y(fpeq:lpeq))/y(req)/Rey + &
                pairst*y(teq)/y(req)**4/Rey - pambst/y(req)/Rey - &
                2/y(req)**2/Ca/Rey - 4*y(req+1)/y(req)**2/Rey - &
                3._dp/2*y(req+1)**2/y(req)
        else
            ydot(req) = sum(y(fpeq:lpeq))*y(req)/4 + &
                pairst*y(teq)/y(req)**2/4 - pambst*y(req)/4 - &
                1._dp/2/Ca   !radius (momentum balance)
        endif
        do i=fpeq,lpeq
            ydot(i) = -3*y(i)*ydot(req)/y(req) + y(i)/y(teq)*ydot(teq) + &
                9*Rg*temp0*D(i-fpeq+1)*R0**5/eta/D0(1)*y(req)*y(teq)*&
                (y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i)*eta*D0(1)/R0**2)/&
                (dz(1)/2)    !partial pressure (molar balance)
        enddo
    else
        do i=fpeq,lpeq
            ydot(i) = -3*y(i)*Rderiv(t)/radius*R0**2/D0(1) + &
                y(i)/y(teq)*ydot(teq) + &
                9*Rg*temp0*D(i-fpeq+1)*R0**5/eta/D0(1)*radius/R0*y(teq)*&
                (y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i)*eta*D0(1)/R0**2)/&
                (dz(1)/2)    !partial pressure (molar balance)
        enddo
    endif
    do j=1,ngas
        do i=1,p+1
            if (i==1) then !bubble boundary
                zw=0e0_dp
                z=dz(i)/2
                ze=dz(i)
                zee=ze+dz(i+1)/2
                lame=(ze-z)/(zee-z)
                c=y(fceq+(i-1)*ngas+j-1)
                cee=y(fceq+i*ngas+j-1)
                cw=KH(j)*y(fpeq+j-1)
                ce=cee*lame+c*(1-lame)
                dcw=(c-cw)/(z-zw)
                dce=(cee-c)/(zee-z)
            elseif(i==p+1) then !outer boundary
                zww=z
                zw=ze
                z=zee
                ze=ze+dz(i)
                lamw=(zw-zww)/(z-zww)
                cww=y(fceq+(i-2)*ngas+j-1)
                c=y(fceq+(i-1)*ngas+j-1)
                cw=c*lamw+cww*(1-lamw)
                ce=c
                dcw=(c-cww)/(z-zww)
                dce=0e0_dp
            else
                zww=z
                zw=ze
                z=zee
                ze=ze+dz(i)
                zee=ze+dz(i+1)/2
                lamw=(zw-zww)/(z-zww)
                lame=(ze-z)/(zee-z)
                cww=y(fceq+(i-2)*ngas+j-1)
                c=y(fceq+(i-1)*ngas+j-1)
                cee=y(fceq+i*ngas+j-1)
                cw=c*lamw+cww*(1-lamw)
                ce=cee*lame+c*(1-lame)
                dcw=(c-cww)/(z-zww)
                dce=(cee-c)/(zee-z)
            endif
            !concentration (molar balance)
            ydot(fceq+(i-1)*ngas+j-1) = 9*D(j)*((ze+radius**3)**(4._dp/3)*dce -&
                (zw+radius**3)**(4._dp/3)*dcw)/dz(i)*R0**2/D0(1)
            if (j==co2_pos) ydot(fceq+(i-1)*ngas+j-1) = &
                ydot(fceq+(i-1)*ngas+j-1) + W0*ydot(xWeq) !reaction source
        enddo
    enddo
end subroutine odesystem_nd
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate dimensional variables
subroutine dim_var_nd
	integer :: i
    time=t*R0**2/D0(1)
    if (firstrun) then
        radius=y(req)*R0 ! calculate bubble radius
    else
        radius=Rb(time) ! use calculated bubble radius
    endif
    temp=y(teq)*temp0
    conv=y(xOHeq)
    call physical_properties(temp,conv,radius)
    eqconc=y(fpeq)*KH(1)  !only first gas
    do i=1,ngas
    	pressure(i)=y(fpeq+i-1)*eta*D0(1)/R0**2
    enddo
    do i=1,ngas
        wblpol(i)=mb(i)*Mbl(i)/(rhop*4*pi/3*(S0**3-R0**3))
    enddo
    avconc=mb/Vsh
    porosity=radius**3/(radius**3+S0**3-R0**3)
    rhofoam=(1-porosity)*rhop
    st=(S0**3+radius**3-R0**3)**(1._dp/3)-radius !thickness of the shell
    pair=pair0*R0**3/radius**3*temp/temp0
end subroutine dim_var_nd
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate non-dimensional variables
subroutine nondim_var
	Rey=rhop*D0(1)/eta
    pairst=pair0*R0**2/eta/D0(1)
    pambst=pamb*R0**2/eta/D0(1)
    Ca=eta*D0(1)/sigma/R0
end subroutine nondim_var
!***********************************END****************************************


!********************************BEGINNING*************************************
!> set initial conditions
subroutine set_initial_conditions_nd
    integer :: i,j
    t = tstart*R0**2/D0(1)
    tout = (t+timestep)*R0**2/D0(1)
    allocate(y(neq))
    y=0
    if (firstrun) then
        y(req)=radius/R0   !radius
        if (inertial_term) y(req+1) = 0        !velocity
    endif
    y(teq) = temp/temp0   !temperature
    y(xOHeq) = conv        !xOH
    y(xWeq) = 0        !xW
    if (kin_model==4) then
        y(kineq(1)) = 6.73000e-02_dp
        y(kineq(2)) = 1.92250e+00_dp
        y(kineq(3)) = 2.26920e+00_dp
        y(kineq(4)) = 0.00000e+00_dp
        y(kineq(5)) = 5.46200e-01_dp
        ! y(kineq(5)) = 1.0924e+00_dp
        y(kineq(6)) = 2.19790e+00_dp
        y(kineq(7)) = 1.64000e+00_dp
        y(kineq(8)) = 1.71030e+00_dp
        y(kineq(9)) = 0.00000e+00_dp
        y(kineq(10)) = 0.00000e+00_dp
        y(kineq(11)) = 0.00000e+00_dp
        y(kineq(12)) = 0.00000e+00_dp
        y(kineq(13)) = 0.00000e+00_dp
        y(kineq(14)) = 0.00000e+00_dp
        y(kineq(15)) = 0.00000e+00_dp
        y(kineq(16)) = 4.45849e+00_dp
        y(kineq(17)) = 0.00000e+00_dp
        y(kineq(18)) = 1.00000e+00_dp
        y(kineq(19)) = 60!2.27000e+01_dp
        y(kineq(20)) = 1e0_dp!8.46382e-01_dp
    endif
    do j=1,ngas
        do i=1,p+1
            y(fceq+(i-1)*ngas+j-1) = cbl(j)      !blowing agent concentration
        enddo
    enddo
    do i=1,ngas
        y(fpeq+i-1) = xgas(i+1)*(pamb+2*sigma/R0)*R0**2/eta/D0(1) !pressure
        if (y(fpeq+i-1)<1e-16_dp) y(fpeq+i-1)=1e-16_dp
    enddo
end subroutine set_initial_conditions_nd
!***********************************END****************************************
end module model_nd
