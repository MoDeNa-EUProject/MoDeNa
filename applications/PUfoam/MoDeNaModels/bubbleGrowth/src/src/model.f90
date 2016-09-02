!> @file
!! contains model subroutines for bubble growth model
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module model
    use constants
    use globals
    use foaming_globals_m
    use fmodena
    use modenastuff
    implicit none
    private
    public molar_balance,kinModel,odesystem,dim_var
contains
!********************************BEGINNING*************************************
!> model supplied to integrator, FVM, nonequidistant mesh
subroutine  odesystem (neq, t, y, ydot)
    use phys_prop, only:Rderiv
    integer :: neq,i,j
    real(dp) :: t,y(neq),ydot(neq),z,zw,ze,zww,zee,lamw,lame,cw,ce,cww,cee,&
        c,dcw,dce,dil,bll
    call dim_var(t,y)
    call molar_balance(y)
    ydot=0
    if (kin_model==4) then
        call kinModel(y)
        do i=1,size(kineq)
            ydot(kineq(i))=kinsource(i)
        enddo
        if (W0>1e-8_dp) then
            ydot(xWeq)=-ydot(kineq(5))/W0*1.0e3_dp
        else
            ydot(xWeq)=0
        endif
        ydot(xOHeq)=-(ydot(kineq(2))+ydot(kineq(3)))/OH0*1.0e3_dp
    else
        ydot(xOHeq) = AOH*exp(-EOH/Rg/temp)*(1-y(xOHeq))*&
            (NCO0-2*W0*y(xWeq)-OH0*y(xOHeq)) !polyol conversion
        if (kin_model==3) then
            if (y(xOHeq)>0.5_dp .and. y(xOHeq)<0.87_dp) ydot(xOHeq)=ydot(xOHeq)*&
                (-2.027_dp*y(xOHeq)+2.013_dp) !gelling influence on kinetics
            if (y(xOHeq)>0.87_dp) ydot(xOHeq)=ydot(xOHeq)*&
                (3.461_dp*y(xOHeq)-2.761_dp)
        endif
        if (W0>1e-3) then
            ! water conversion
            ! ydot(xWeq) = AW*exp(-EW/Rg/temp)*(1-y(xWeq))*&
            !     (NCO0-2*W0*y(xWeq)-OH0*y(xOHeq)) 2nd order
            ydot(xWeq) = AW*exp(-EW/Rg/temp)*(1-y(xWeq)) !1st order
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
    endif
    !temperature (enthalpy balance)
    ydot(teq) = -dHOH*OH0/(rhop*cp)*ydot(xOHeq)-dHW*W0/(rhop*cp)*ydot(xWeq)
    do i=1,ngas
        if (geometry=="3D") then
            ydot(teq) = ydot(teq) - dHv(i)*12*pi*Mbl(i)*D(i)*radius**4/&
                (rhop*cp*Vsh)*(y(fceq+i-1)-KH(i)*y(fpeq+i-1))/(dz(1)/2)
        elseif (geometry=="2D") then
            ydot(teq) = ydot(teq) - dHv(i)*4*pi*Mbl(i)*D(i)*radius**2/&
                (rhop*cp*Vsh)*(y(fceq+i-1)-KH(i)*y(fpeq+i-1))/(dz(1)/2)
        endif
    enddo
    if (kin_model==4) then
        ydot(kineq(19))=ydot(teq)
    endif
    if (firstrun) then
        if (inertial_term) then
            if (geometry=="3D") then
                ydot(req) = y(req+1)    !radius (momentum balance)
                ydot(req+1) = (sum(y(fpeq:lpeq)) + pair - pamb - &
                    2*sigma/radius - 4*eta*y(req+1)/radius - &
                    3._dp/2*y(req+1)**2)/(radius*rhop)
            elseif (geometry=="2D") then
                ydot(req) = y(req+1)    !radius (momentum balance)
                ydot(req+1) = (sum(y(fpeq:lpeq)) + pair - pamb - &
                    sigma/radius - 2*eta*y(req+1)/radius - &
                    3._dp/2*y(req+1)**2)/(radius*rhop)
            endif
        else
            if (geometry=="3D") then
                ydot(req) = (sum(y(fpeq:lpeq)) + pair - pamb - &
                    2*sigma/radius)*radius/(4*eta)   !radius (momentum balance)
            elseif (geometry=="2D") then
                ydot(req) = (sum(y(fpeq:lpeq)) + pair - pamb - &
                    sigma/radius)*radius/(2*eta)   !radius (momentum balance)
            endif
        endif
        !partial pressure (molar balance)
        do i=fpeq,lpeq
            if (geometry=="3D") then
                ydot(i) = -3*y(i)*ydot(req)/radius + y(i)/temp*ydot(teq) + &
                    9*Rg*temp*D(i-fpeq+1)*radius*&
                    (y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i))/(dz(1)/2)
            elseif (geometry=="2D") then
                ydot(i) = -2*y(i)*ydot(req)/radius + y(i)/temp*ydot(teq) + &
                    4*Rg*temp*D(i-fpeq+1)*&
                    (y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i))/(dz(1)/2)
            endif
        enddo
    else
        !partial pressure (molar balance)
        do i=fpeq,lpeq
            if (geometry=="3D") then
                ydot(i) = -3*y(i)*Rderiv(time)/radius + y(i)/temp*ydot(teq) + &
                    9*Rg*temp*D(i-fpeq+1)*radius*&
                    (y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i))/(dz(1)/2)
            elseif (geometry=="2D") then
                ydot(i) = -2*y(i)*Rderiv(time)/radius + y(i)/temp*ydot(teq) + &
                    4*Rg*temp*D(i-fpeq+1)*&
                    (y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i))/(dz(1)/2)
            endif
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
            if (geometry=="3D") then
                ydot(fceq+(i-1)*ngas+j-1) = &
                    9*D(j)*((ze+radius**3)**(4._dp/3)*dce -&
                    (zw+radius**3)**(4._dp/3)*dcw)/dz(i)
                !reaction source
                if (j==co2_pos) ydot(fceq+(i-1)*ngas+j-1) = &
                    ydot(fceq+(i-1)*ngas+j-1) + W0*ydot(xWeq)
            elseif (geometry=="2D") then
                ydot(fceq+(i-1)*ngas+j-1) = &
                    4*D(j)*((ze+radius**2)*dce -&
                    (zw+radius**2)*dcw)/dz(i)
                !reaction source
                !WARNING: reaction source exaggerated to provide reasonable
                !growth
                if (j==co2_pos) ydot(fceq+(i-1)*ngas+j-1) = &
                    ydot(fceq+(i-1)*ngas+j-1) + W0*ydot(xWeq)*1e4
            endif
        enddo
    enddo
end subroutine odesystem
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates molar amount of blowing agents in bubble and shell
subroutine molar_balance(y)
    integer :: i,j
    real(dp) :: y(:)
    !numerical integration
    mb=0e0_dp
    !rectangle rule
    do i=1,p+1
        do j=1,ngas
            mb(j)=mb(j)+y(fceq+j-1+(i-1)*ngas)*dz(i)
        enddo
    enddo
    if (geometry=="3D") then
        mb=mb*4*pi/3 !moles in polymer
    elseif (geometry=="2D") then
        mb=mb*pi !moles in polymer
    endif
    do i=1,ngas
        if (geometry=="3D") then
            mb2(i)=pressure(i)*radius**3*4*pi/(3*Rg*temp) !moles in bubble
        elseif (geometry=="2D") then
            mb2(i)=pressure(i)*radius**2*pi/(Rg*temp) !moles in bubble
        endif
    enddo
    mb3=mb+mb2 !total moles
end subroutine molar_balance
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate dimensional variables
subroutine dim_var(t,y)
    use phys_prop, only:physical_properties,Rb
	integer :: i
    real(dp) :: t,y(:)
    time=t
    if (firstrun) then
        radius=y(req) ! calculate bubble radius
    else
        radius=Rb(time) ! use calculated bubble radius
    endif
    temp=y(teq)
    conv=y(xOHeq)
    call physical_properties(temp,conv,radius)
    eqconc=y(fpeq)*KH(1)  !only first gas
    do i=1,ngas
    	pressure(i)=y(fpeq+i-1)
    enddo
    do i=1,ngas
        wblpol(i)=mb(i)*Mbl(i)/(rhop*Vsh)
    enddo
    avconc=mb/Vsh
    if (geometry=="3D") then
        laplace_pres=2*sigma/radius
        pair=pair0*R0**3/radius**3*temp/temp0
        porosity=radius**3/(radius**3+S0**3-R0**3)
        st=(S0**3+radius**3-R0**3)**(1._dp/3)-radius !thickness of the shell
    elseif (geometry=="2D") then
        laplace_pres=sigma/radius
        pair=pair0*R0**2/radius**2*temp/temp0
        porosity=radius**2/(radius**2+S0**2-R0**2)
        st=(S0**2+radius**2-R0**2)**(1._dp/2)-radius !thickness of the shell
    endif
    rhofoam=(1-porosity)*rhop
    bub_pres=sum(pressure)+pair-pamb
end subroutine dim_var
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates kinetic source terms
!! modena models
subroutine kinModel(y)
    real(dp), intent(in) :: y(:)
    integer :: i
    if (kin_model==4) then
        do i=1,size(kineq)
            call modena_inputs_set(kinInputs, kinInputsPos(i+1), y(kineq(i)))
        enddo
    endif
    call modena_inputs_set(kinInputs, kinInputsPos(1), time)
    call modena_inputs_set(kinInputs, kinInputsPos(20), y(teq)-273.15_dp)
    ret = modena_model_call (kinModena, kinInputs, kinOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    if (kin_model==4) then
        kinsource=0
        do i=1,size(kineq)
            kinsource(i) = modena_outputs_get(kinOutputs, kinOutputsPos(i))
        enddo
    endif
end subroutine kinModel
!***********************************END****************************************
end module model
