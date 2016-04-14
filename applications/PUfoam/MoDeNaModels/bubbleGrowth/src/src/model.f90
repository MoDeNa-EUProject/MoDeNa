!> @file
!! contains model subroutines for bubble growth model
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module model
    use iso_c_binding
    use constants
    use foaming_globals_m
    use fmodena
    use modenastuff
    use in_out
    implicit none
    private
    !time integration variables for lsode
    integer :: iout, iopt, istate, itask, itol, liw, lrw, nnz, lenrat, neq, mf
    real(dp) :: jac,tout,rtol,atol,t
    real(dp), dimension(:), allocatable :: rwork!,y
    integer, dimension(:), allocatable :: iwork
    !mesh variables
    integer :: info
    real(dp),allocatable :: atri(:),btri(:),ctri(:),rtri(:),utri(:),dz(:)
    !needed for selection of subroutine for evaluation of derivatives
    abstract interface
        subroutine sub (neq, t, y, ydot)
            use constants
            integer :: neq
            real(dp) ::  t, y(neq), ydot(neq)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => odesystem
    public bblpreproc,bblinteg
contains
!********************************BEGINNING*************************************
!> model supplied to integrator, FVM, nonequidistant mesh
subroutine  odesystem (neq, t, y, ydot)
    integer :: neq,i,j
    real(dp) :: t,y(neq),ydot(neq),z,zw,ze,zww,zee,lamw,lame,cw,ce,cww,cee,&
        c,dcw,dce,dil,bll
    call dim_var
    call molar_balance
    ydot=0
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
    if (kin_model==4) then
        call kinModel
        do i=1,size(kineq)
            ydot(kineq(i))=kinsource(i)
        enddo
    endif
    !temperature (enthalpy balance)
    ydot(teq) = -dHOH*OH0/(rhop*cp)*ydot(xOHeq)-dHW*W0/(rhop*cp)*ydot(xWeq)
    do i=1,ngas
        ydot(teq) = ydot(teq) - dHv(i)*12*pi*Mbl(i)*D(i)*radius**4/&
            (rhop*cp*Vsh)*(y(fceq+i-1)-KH(i)*y(fpeq+i-1))/(dz(1)/2)
    enddo
    if (kin_model==4) then
        ! ydot(kineq(19))=ydot(teq)
    endif
    if (firstrun) then
        if (inertial_term) then
            ydot(req) = y(req+1)    !radius (momentum balance)
            ydot(req+1) = (sum(y(fpeq:lpeq)) + pair - pamb - &
                2*sigma/radius - 4*eta*y(req+1)/radius - &
                3._dp/2*y(req+1)**2)/(radius*rhop)
        else
            ydot(req) = (sum(y(fpeq:lpeq)) + pair - pamb - &
                2*sigma/radius)*radius/(4*eta)   !radius (momentum balance)
        endif
        do i=fpeq,lpeq
            ydot(i) = -3*y(i)*ydot(req)/radius + y(i)/temp*ydot(teq) + &
                9*Rg*temp*D(i-fpeq+1)*radius*(y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i))/&
                (dz(1)/2)    !partial pressure (molar balance)
        enddo
    else
        do i=fpeq,lpeq
            ydot(i) = -3*y(i)*Rderiv(t)/radius + y(i)/temp*ydot(teq) + &
                9*Rg*temp*D(i-fpeq+1)*radius*(y(fceq+i-fpeq)-KH(i-fpeq+1)*y(i))/&
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
                (zw+radius**3)**(4._dp/3)*dcw)/dz(i)
            if (j==co2_pos) ydot(fceq+(i-1)*ngas+j-1) = &
                ydot(fceq+(i-1)*ngas+j-1) + W0*ydot(xWeq) !reaction source
        enddo
    enddo
end subroutine odesystem
!***********************************END****************************************


!********************************BEGINNING*************************************
!> model supplied to integrator, FVM, nonequidistant mesh
subroutine  odesystem_nd (neq, t, y, ydot)
    integer :: neq,i,j
    real(dp) :: t,y(neq),ydot(neq),z,zw,ze,zww,zee,lamw,lame,cw,ce,cww,cee,&
        c,dcw,dce,dil,bll
    call dim_var_nd
    call molar_balance
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
        call kinModel
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
!> calculates values of physical properties
subroutine physical_properties(temp,conv)
    real(dp), intent(in) :: temp,conv
    integer :: i
    if (.not. gelpoint .and. temp<500) then
        select case(visc_model)
        case(1)
        case(2)
            eta=Aeta*exp(Eeta/(Rg*temp))*(Cg/(Cg-conv))**(AA+B*conv)
        case(3)
            !set input vector
            call modena_inputs_set(viscInputs, viscTpos, temp);
            call modena_inputs_set(viscInputs, viscXPos, conv);
            !call model
            ret = modena_model_call(viscModena, viscInputs, viscOutputs)
            if(ret /= 0) then
                call exit(ret)
            endif
            !fetch results
            eta = modena_outputs_get(viscOutputs, 0_c_size_t);
        end select
        if (eta>maxeta .or. isnan(eta)) then
            eta=maxeta
            gelpoint=.true.
            write(*,'(2x,A,es8.2,A)') 'gel point reached at time t = ',TOUT,' s'
            write(*,'(2x,A,es8.2,A)') 'temperature at gel point T = ',temp,' K'
            write(*,'(2x,A,es8.2)') 'conversion at gel point X = ',conv
        endif
    else
        eta=maxeta
    endif
    select case(itens_model)
    case(1)
    case(2)
        call modena_inputs_set(itensInputs, itensTpos, temp)
        ret = modena_model_call(itensModena, itensInputs, itensOutputs)
        if(ret /= 0) then
            call exit(ret)
        endif
        sigma = modena_outputs_get(itensOutputs, 0_c_size_t)*1e-3_dp
    end select
    do i=1,ngas
        select case(diff_model(i))
        case(1)
        case(2)
            call modena_inputs_set(diffInputs(i), diffTpos(i), temp)
            ret = modena_model_call(diffModena(i), diffInputs(i), diffOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            D(i) = modena_outputs_get(diffOutputs(i), 0_c_size_t)
        end select
        select case(sol_model(i))
        case(1)
        case(2)
            ! TODO: implement properly
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            call modena_inputs_set(solInputs(i), solXgasPos(i), 1.0e-4_dp)
            call modena_inputs_set(solInputs(i), solXmdiPos(i), 0.5_dp)
            call modena_inputs_set(solInputs(i), solXpolyolPos(i), 0.5_dp)
            ret = modena_model_call(solModena(i), solInputs(i), solOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            KH(i) = modena_outputs_get(solOutputs(i), 0_c_size_t)
            KH(i)=rhop/Mbl(i)/KH(i)
        case(3)
            KH(i)=-rhop/Mbl(i)/pamb*3.3e-4_dp*(exp((2.09e4_dp-67.5_dp*(temp-&
                35.8_dp*log(pamb/1e5_dp)))/(8.68e4_dp-(temp-35.8_dp*&
                log(pamb/1e5_dp))))-1.01_dp)**(-1)
        case(4)
            KH(i)=rhop/Mbl(i)/pamb*(0.0064_dp+0.0551_dp*exp(-(temp-298)**2/&
                (2*17.8_dp**2)))
        case(5)
            KH(i)=rhop/Mbl(i)/pamb*(0.00001235_dp*temp**2-0.00912_dp*temp+&
                1.686_dp)
        case(6)
            KH(i)=rhop/Mbl(i)/pamb*(1e-7_dp+4.2934_dp*&
                exp(-(temp-203.3556_dp)**2/(2*40.016_dp**2)))
        end select
    enddo
    if (solcorr) KH=KH*exp(2*sigma*Mbl/(rhop*Rg*temp*radius))
    cp=cppol+sum(cbl*Mbl*cpbll)/rhop
end subroutine physical_properties
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates molar amount of blowing agents in bubble and shell
subroutine molar_balance
    integer :: i,j
    !numerical integration
    mb=0e0_dp
    !rectangle rule
    do i=1,p+1
        do j=1,ngas
            mb(j)=mb(j)+y(fceq+j-1+(i-1)*ngas)*dz(i)
        enddo
    enddo
    mb=mb*4*pi/3 !moles in polymer
    do i=1,ngas
    	mb2(i)=pressure(1)*radius**3*4*pi/(3*Rg*temp) !moles in bubble
    enddo
    mb3=mb+mb2 !total moles
end subroutine molar_balance
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate dimensional variables
subroutine dim_var
	integer :: i
    time=t
    if (firstrun) then
        radius=y(req) ! calculate bubble radius
    else
        radius=Rb(time) ! use calculated bubble radius
    endif
    temp=y(teq)
    conv=y(xOHeq)
    call physical_properties(temp,conv)
    eqconc=y(fpeq)*KH(1)  !only first gas
    do i=1,ngas
    	pressure(i)=y(fpeq+i-1)
    enddo
    do i=1,ngas
        wblpol(i)=mb(i)*Mbl(i)/(rhop*4*pi/3*(S0**3-R0**3))
    enddo
    avconc=mb/Vsh
    porosity=radius**3/(radius**3+S0**3-R0**3)
    rhofoam=(1-porosity)*rhop
    st=(S0**3+radius**3-R0**3)**(1._dp/3)-radius !thickness of the shell
    pair=pair0*R0**3/radius**3*temp/temp0
end subroutine dim_var
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
    call physical_properties(temp,conv)
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
!> calculate growth rate
subroutine growth_rate
	integer :: i
    do i=1,ngas
        grrate(i)=(mb2(i)-nold(i))/timestep
    enddo
    do i=1,ngas
        nold(i)=mb2(i)
    enddo
end subroutine growth_rate
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates kinetic source terms
!! modena models
subroutine kinModel
    integer :: i
    if (kin_model==4) then
        do i=1,size(kineq)
            call modena_inputs_set(kinInputs, kinInputsPos(i), y(kineq(i)))
        enddo
    endif
    !TODO: implement temperature
    call modena_inputs_set(kinInputs, kinInputsPos(19), 60.0_dp)
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


!********************************BEGINNING*************************************
!> prepares integration
subroutine bblpreproc
    integer :: i,j
    write(*,*) 'preparing simulation...'
    select case (kin_model)
    case(1)
    case(3)
    case(4)
    case default
        stop 'unknown kinetic model'
    end select
    !determine number of equations and their indexes
    neq=(p+1)*ngas
    if (firstrun) then
        neq = neq+4+ngas
        req=1   !radius index
        fpeq=2  !pressure index
        if (inertial_term) then
            neq=neq+1
            fpeq=fpeq+1
        endif
    else
        neq = neq+3+ngas
        fpeq=1  !pressure index
    endif
    lpeq=fpeq+ngas-1
    teq=lpeq+1   !temperature index
    xOHeq=teq+1 !polyol conversion index
    xWeq=xOHeq+1  !water conversion index
    fceq = xWeq+1    !concentration index
    if (kin_model==4) then
        allocate(kineq(20),kinsource(20))
    endif
    if (kin_model==4) then
        neq=neq+size(kineq)
        do i=1,size(kineq)
            kineq(i)=xWeq+i
        enddo
        fceq=kineq(size(kineq))+1
    endif

    ! determine physical properties
    radius = R0
    temp=temp0
    conv=0.0_dp
    if (sum(xgas) /= 1) then
        write(*,*) 'Sum of initial molar fractions of gases in the bubble is &
            not equal to one. Normalizing...'
        xgas=xgas/sum(xgas)
        write(*,*) 'New initial molar fractions of gases in the bubble'
        write(*,*) xgas
    endif
    call createModenaModels
    select case(rhop_model) !density is kept constant, calculate it only once
    case(1)
    case(2)
        call modena_inputs_set(rhopInputs, rhopTpos, temp)
        call modena_inputs_set(rhopInputs, rhopXOHPos, 0.1_dp)
        ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
        if(ret /= 0) then
            call exit(ret)
        endif
        rhop = modena_outputs_get(rhopOutputs, 0_c_size_t)
        call modena_inputs_set(rhopInputs, rhopTpos, temp+100)
        call modena_inputs_set(rhopInputs, rhopXOHPos, 0.9_dp)
        ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
        if(ret /= 0) then
            call exit(ret)
        endif
        !average density during foaming
        rhop=(rhop + modena_outputs_get(rhopOutputs, 0_c_size_t))/2
    end select
    call physical_properties(temp,conv)
    D0=D
    surface_tension=sigma
    pair0=(pamb+2*sigma/R0)*xgas(1)
    S0=Sn*radius
    Vsh=4*pi/3*(S0**3-R0**3)
    gelpoint=.false.
    timestep=(tend-tstart)/its
    ! write(*,'(2x,A,2x,e12.6)') 'NN',Sn**(-3)/(1-Sn**(-3))/&
    !     exp(log(4._dp/3*pi*R0**3))
    if (firstrun) then
        allocate(etat(its,2),port(its,2),init_bub_rad(its,2))
    else
        ! call load_old_results
    endif

    !calculate spatial grid points
    allocate(atri(p),btri(p),ctri(p),rtri(p),utri(p),dz(p+1))
    atri=-mshco !lower diagonal
    btri=1+mshco    !main diagonal
    ctri=-1 !upper diagonal
    rtri=0  !rhs
    rtri(p)=S0**3-R0**3
    utri=rtri
    call dgtsl(p,atri,btri,ctri,utri,info)
    if ((utri(2)-utri(1))/utri(1)<10*epsilon(utri)) stop 'set smaller mesh &
        coarsening parameter'
    dz(1)=utri(1)+rtri(1)/mshco
    do i=2,p
        dz(i)=utri(i)-utri(i-1)
    enddo
    dz(p+1)=rtri(p)-utri(p)
    deallocate(atri,btri,ctri,rtri,utri)

    call set_init

    !choose and set integrator
    mf=int_meth
    select case(integrator)
    case(1)
        select case(mf)
        case(10)
            allocate(rwork(20+16*neq),iwork(20))
        case(22)
            allocate(rwork(22+9*neq+neq**2),iwork(20+neq))
        case default
            stop 'unknown mf'
        end select
    case(2)
        select case(mf)
        case(10)
            allocate(rwork(20+16*neq),iwork(30))
        case(222)
            nnz=neq**2 !not sure, smaller numbers make problems for low p
            lenrat=2 !depends on dp
            allocate(rwork(int(20+(2+1._dp/lenrat)*nnz+(11+9._dp/lenrat)*neq)),&
                iwork(30))
        case default
            stop 'unknown mf'
        end select
    case default
        stop 'unknown integrator'
    end select
    itask = 1
    istate = 1
    iopt = 1
    rwork(5:10)=0
    iwork(5:10)=0
    lrw = size(rwork)
    liw = size(iwork)
    iwork(6)=maxts
    itol = 2 !don't change, or you must declare atol as atol(neq)
    rtol=rel_tol
    atol=abs_tol
    allocate(atol2(neq))
    atol2=atol
    atol2(req)=atol2(req)!/1.e-1_dp
    do i=fpeq,lpeq
        atol2(i)=atol2(i)!*1.e-2_dp
    enddo
    write(*,*) 'done: simulation prepared'
    write(*,*)
end subroutine bblpreproc
!***********************************END****************************************


!********************************BEGINNING*************************************
!> set initial conditions
subroutine set_init
    integer :: i,j
    t = tstart
    tout = t+timestep
    allocate(y(neq))
    y=0
    if (firstrun) then
        y(req)=radius   !radius
        if (inertial_term) y(req+1) = 0        !velocity
    endif
    y(teq) = temp   !temperature
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
        y(fpeq+i-1) = xgas(i+1)*(pamb+2*sigma/R0) !pressure
        if (y(fpeq+i-1)<1e-16_dp) y(fpeq+i-1)=1e-16_dp
    enddo
end subroutine set_init
!***********************************END****************************************


!********************************BEGINNING*************************************
!> set initial conditions
subroutine set_init_nd
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
end subroutine set_init_nd
!***********************************END****************************************


!********************************BEGINNING*************************************
!> performs integration
subroutine bblinteg
    write(*,*) 'integrating...'
    if (firstrun) then
        call save_integration_header
    endif
    do iout = 1,its
        select case (integrator)
        case(1)
            call dlsode (sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
                istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        case(2)
            call dlsodes (sub_ptr, neq, y, t, tout, itol, rtol, atol2, itask, &
                istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        case default
            stop 'unknown integrator'
        end select
        call dim_var
        call molar_balance
        call growth_rate
        if (firstrun) then
            call save_integration_step(iout)
        endif
        ! write(*,*) tout,kinsource(2)
        write(*,'(2x,A4,F8.3,A3,A13,F10.3,A4,A25,F8.3,A4,A9,EN12.3,A4)') &
            't = ', t, ' s,',&
            'p_b - p_o = ', sum(pressure)+pair-pamb, ' Pa,', &
            'p_b - p_o - p_Laplace = ', sum(pressure)+pair-pamb-2*sigma/radius, ' Pa,',&
            'dR/dt = ', (sum(pressure)+pair-pamb-2*sigma/radius)*radius/4/eta, ' m/s'
        ! write(*,*) tout, radius**3/(radius**3+S0**3-R0**3), radius, eta
        tout = t+timestep
        if (gelpoint) exit
    end do
    if (firstrun) then
        call save_integration_close(iout)
    else
        bub_pres=sum(pressure)+pair-pamb
    endif
    write(*,*) 'done: integration'
    call destroyModenaModels
    deallocate(D,cbl,xgas,KH,fic,Mbl,dHv,mb,mb2,mb3,avconc,pressure,&
        diff_model,sol_model,cpblg,cpbll,RWORK,IWORK,dz,Y,atol2,wblpol,D0)
end subroutine bblinteg
!***********************************END****************************************


!********************************BEGINNING*************************************
!> time derivation of bubble radius as function of time
real(dp) function Rderiv(t)
    real(dp) :: t,dt=1.e-4_dp
    ! Rderiv=(Rb(t+dt)-Rb(t))/dt
    ! Rderiv=(-0.5_dp*Rb(t+2*dt)+2*Rb(t+dt)-1.5_dp*Rb(t))/dt
    Rderiv=(2*Rb(t+3*dt)-9*Rb(t+2*dt)+18*Rb(t+dt)-11*Rb(t))/(6*dt)
endfunction Rderiv
!***********************************END****************************************


!********************************BEGINNING*************************************
!> bubble radius as function of time
real(dp) function Rb(t)
    use interpolation
    real(dp) :: t
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=t
    call pwl_interp_1d ( size(bub_rad(:,1)), bub_rad(:,1), &
        bub_rad(:,bub_inx+1), ni, xi, yi )
    Rb=yi(1)
endfunction Rb
!***********************************END****************************************
end module model
