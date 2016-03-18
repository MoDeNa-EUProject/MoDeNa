!> @file
!! contains model subroutines for bubble growth model
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module model
    use constants
    use in_out
    use iso_c_binding
    use fmodena
    use modenastuff
    implicit none
    integer :: info
    real(dp) :: Pair0,timestep,GR,Rold,Told,pold(2),nold(2),Vsh
    !time integration variables for lsode
    integer :: IOUT, IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NNZ, LENRAT!, MF, NEQ
    real(dp) :: JAC,TOUT!,RTOL,ATOL,T
    real(dp), dimension(:), allocatable :: RWORK!,Y
    integer, dimension(:), allocatable :: IWORK
    !mesh variables
    real(dp),allocatable :: atri(:),btri(:),ctri(:),rtri(:),utri(:),dz(:)
    !needed for selection of subroutine for evaluation of derivatives
    abstract interface
        subroutine sub (NEQ, T, Y, YDOT)
            use constants
            INTEGER :: NEQ
            real(dp) ::  T, Y(NEQ), YDOT(NEQ)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => FEX
contains
!********************************BEGINNING*************************************
!> model supplied to integrator, FVM, nonequidistant mesh
SUBROUTINE  FEX (NEQ, T, Y, YDOT)
    INTEGER :: NEQ,i,j
    real(dp) ::  T, Y(NEQ), YDOT(NEQ),z,zw,ze,zww,zee,lamw,lame,cw,ce,cww,cee,&
        c,dcw,dce,dil,bll
    call molar_balance
    YDOT=0
    YDOT(xOHeq) = AOH*exp(-EOH/Rg/Y(teq))*(1-Y(xOHeq))*&
        (NCO0-2*W0*Y(xWeq)-OH0*Y(xOHeq)) !polyol conversion
    if (kin_model==3) then
        if (Y(xOHeq)>0.5_dp .and. Y(xOHeq)<0.87_dp) YDOT(xOHeq)=YDOT(xOHeq)*&
            (-2.027_dp*Y(xOHeq)+2.013_dp) !gelling influence on kinetics
        if (Y(xOHeq)>0.87_dp) YDOT(xOHeq)=YDOT(xOHeq)*&
            (3.461_dp*Y(xOHeq)-2.761_dp)
    endif
    if (W0>1e-3) then
        ! water conversion
        ! YDOT(xWeq) = AW*exp(-EW/Rg/Y(teq))*(1-Y(xWeq))*&
        !     (NCO0-2*W0*Y(xWeq)-OH0*Y(xOHeq)) 2nd order
        YDOT(xWeq) = AW*exp(-EW/Rg/Y(teq))*(1-Y(xWeq)) !1st order
    endif
    if (dilution) then
        if (co2_pos==1) then
            bll=mb(2)/Vsh*Mbl(2)/rhop
        else
            bll=mb(1)/Vsh*Mbl(1)/rhop
        endif
        dil=1/(1+rhop/rhobl*bll)
        YDOT(xOHeq)=YDOT(xOHeq)*dil
        YDOT(xWeq)=YDOT(xWeq)*dil
    endif
    if (kin_model==2 .or. kin_model==4) then
        call kinModel
        do i=1,size(kineq)
            YDOT(kineq(i))=kinsource(i)
        enddo
        ! YDOT(xOHeq) = -YDOT(kineq(2))/OH0
        ! YDOT(xWeq) = -YDOT(kineq(3))/W0
    endif
    !temperature (enthalpy balance)
    YDOT(teq) = -dHOH*OH0/(rhop*cp)*YDOT(xOHeq)-dHW*W0/(rhop*cp)*YDOT(xWeq)
    do i=1,ngas
        YDOT(teq) = YDOT(teq) - dHv(i)*12*pi*Mbl(i)*D(i)*Y(req)**4/&
            (rhop*cp*Vsh)*(Y(fceq+i-1)-KH(i)*Y(fpeq+i-1))/(dz(1)/2)
    enddo
    if (kin_model==2) then
        YDOT(kineq(12))=YDOT(teq)
    elseif (kin_model==4) then
        ! YDOT(kineq(19))=YDOT(teq)
    endif
    if (inertial_term) then
        YDOT(req) = Y(req+1)    !radius (momentum balance)
        YDOT(req+1) = (sum(Y(fpeq:lpeq)) + Pair0*R0**3/Y(req)**3 - Pamb - &
            2*sigma/Y(req) - 4*eta*Y(req+1)/Y(req) - &
            3._dp/2*Y(req+1)**2)/(Y(req)*rhop)
    else
        YDOT(req) = (sum(Y(fpeq:lpeq)) + Pair0*R0**3/Y(req)**3 - Pamb - &
            2*sigma/Y(req))*Y(req)/(4*eta)   !radius (momentum balance)
    endif
    do i=fpeq,lpeq
        YDOT(i) = -3*Y(i)*YDOT(req)/Y(req) + Y(i)/Y(teq)*YDOT(teq) + &
            9*Rg*Y(teq)*D(i-fpeq+1)*Y(req)*(Y(fceq+i-fpeq)-KH(i-fpeq+1)*Y(i))/&
            (dz(1)/2)    !partial pressure (molar balance)
    enddo
    do j=1,ngas
        do i=1,p+1
            if (i==1) then !bubble boundary
                zw=0e0_dp
                z=dz(i)/2
                ze=dz(i)
                zee=ze+dz(i+1)/2
                lame=(ze-z)/(zee-z)
                c=Y(fceq+(i-1)*ngas+j-1)
                cee=Y(fceq+i*ngas+j-1)
                cw=KH(j)*Y(fpeq+j-1)
                ce=cee*lame+c*(1-lame)
                dcw=(c-cw)/(z-zw)
                dce=(cee-c)/(zee-z)
            elseif(i==p+1) then !outer boundary
                zww=z
                zw=ze
                z=zee
                ze=ze+dz(i)
                lamw=(zw-zww)/(z-zww)
                cww=Y(fceq+(i-2)*ngas+j-1)
                c=Y(fceq+(i-1)*ngas+j-1)
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
                cww=Y(fceq+(i-2)*ngas+j-1)
                c=Y(fceq+(i-1)*ngas+j-1)
                cee=Y(fceq+i*ngas+j-1)
                cw=c*lamw+cww*(1-lamw)
                ce=cee*lame+c*(1-lame)
                dcw=(c-cww)/(z-zww)
                dce=(cee-c)/(zee-z)
            endif
            !concentration (molar balance)
            YDOT(fceq+(i-1)*ngas+j-1) = 9*D(j)*((ze+Y(req)**3)**(4._dp/3)*dce -&
                (zw+Y(req)**3)**(4._dp/3)*dcw)/dz(i)
            if (j==co2_pos) YDOT(fceq+(i-1)*ngas+j-1) = &
                YDOT(fceq+(i-1)*ngas+j-1) + W0*YDOT(xWeq) !reaction source
        enddo
    enddo
END subroutine FEX
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates values of physical properties
subroutine physical_properties(Y)
    real(dp), dimension(:), intent(in) :: Y
    integer :: i
    if (.not. gelpoint .and. Y(teq)<500) then
        select case(visc_model)
        case(1)
        case(2)
            eta=Aeta*exp(Eeta/(Rg*Y(teq)))*(Cg/(Cg-Y(xOHeq)))**(AA+B*Y(xOHeq))
        case(3)
            !set input vector
            call modena_inputs_set(viscInputs, viscTpos, Y(teq));
            call modena_inputs_set(viscInputs, viscXPos, Y(xOHeq));
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
            write(*,'(2x,A,es8.2,A)') 'temperature at gel point T = ',Y(teq),' K'
            write(*,'(2x,A,es8.2)') 'conversion at gel point X = ',Y(xOHeq)
        endif
    else
        eta=maxeta
    endif
    select case(itens_model)
    case(1)
    case(2)
        call modena_inputs_set(itensInputs, itensTpos, Y(teq))
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
            call modena_inputs_set(diffInputs(i), diffTpos(i), Y(teq))
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
            call modena_inputs_set(solInputs(i), solTpos(i), Y(teq))
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
            KH(i)=-rhop/Mbl(i)/Pamb*3.3e-4_dp*(exp((2.09e4_dp-67.5_dp*(Y(teq)-&
                35.8_dp*log(Pamb/1e5_dp)))/(8.68e4_dp-(Y(teq)-35.8_dp*&
                log(Pamb/1e5_dp))))-1.01_dp)**(-1)
        case(4)
            KH(i)=rhop/Mbl(i)/Pamb*(0.0064_dp+0.0551_dp*exp(-(Y(teq)-298)**2/&
                (2*17.8_dp**2)))
        case(5)
            KH(i)=rhop/Mbl(i)/Pamb*(0.00001235_dp*Y(teq)**2-0.00912_dp*Y(teq)+&
                1.686_dp)
        case(6)
            KH(i)=rhop/Mbl(i)/Pamb*(1e-7_dp+4.2934_dp*&
                exp(-(Y(teq)-203.3556_dp)**2/(2*40.016_dp**2)))
        end select
    enddo
    if (solcorr) KH=KH*exp(2*sigma*Mbl/(rhop*Rg*Y(teq)*Y(req)))
    cp=cppol+sum(cbl*Mbl*cpbll)/rhop
end subroutine physical_properties
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates molar amount in bubble and shell and thickness of the shell
subroutine molar_balance
    integer :: i,j
    call physical_properties(Y)
    !numerical integration
    mb=0e0_dp
    !rectangle rule
    do i=1,p+1
        do j=1,ngas
            mb(j)=mb(j)+Y(fceq+j-1+(i-1)*ngas)*dz(i)
        enddo
    enddo
    mb=mb*4*pi/3 !moles in polymer
    do i=1,ngas
    	mb2(i)=Y(fpeq+i-1)*Y(req)**3*4*pi/(3*Rg*Y(teq)) !moles in bubble
    enddo
    mb3=mb+mb2 !total moles
    st=(S0**3+Y(req)**3-R0**3)**(1._dp/3)-Y(req) !thickness of the shell
end subroutine molar_balance
!***********************************END****************************************


!********************************BEGINNING*************************************
!> restores dimensional variables
subroutine restoreDV
	integer :: i
    time=TOUT
    radius=Y(req)
    eqconc=Y(fpeq)*KH(1)  !only first gas
    do i=1,ngas
    	pressure(i)=Y(fpeq+i-1)
        grrate(i)=(mb2(i)-nold(i))/timestep
    enddo
    i=1
    Rold=Y(req)
    Told=Y(teq)
    do i=1,ngas
        pold(i)=Y(fpeq+i-1)
        nold(i)=mb2(i)
    enddo
    avconc=mb/Vsh
end subroutine restoreDV
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates kinetic source terms
!! modena models
subroutine kinModel
    integer :: i
    if (kin_model==2) then
        call modena_inputs_set(kinInputs, kinNCOPos, Y(kineq(1)));
    	call modena_inputs_set(kinInputs, kinOHPos, Y(kineq(2)));
    	call modena_inputs_set(kinInputs, kinH2OPos, Y(kineq(3)));
        call modena_inputs_set(kinInputs, kinCO2Pos, Y(kineq(4)));
        call modena_inputs_set(kinInputs, kinPentanePos, Y(kineq(5)));
        call modena_inputs_set(kinInputs, kinPolymerPos, Y(kineq(6)));
        call modena_inputs_set(kinInputs, kinPolymerBlowPos, Y(kineq(7)));
        call modena_inputs_set(kinInputs, kinUreaPos, Y(kineq(8)));
        call modena_inputs_set(kinInputs, kinR1Pos, Y(kineq(9)));
        call modena_inputs_set(kinInputs, kinRmassPos, Y(kineq(10)));
        call modena_inputs_set(kinInputs, kinRvolPos, Y(kineq(11)));
        call modena_inputs_set(kinInputs, kinRtempPos, Y(kineq(12)));
    elseif (kin_model==4) then
        do i=1,size(kineq)
            call modena_inputs_set(kinInputs, kinInputsPos(i), Y(kineq(i)))
        enddo
    endif
    call modena_inputs_set(kinInputs, kinInputsPos(19), 60.0_dp)
    ret = modena_model_call (kinModena, kinInputs, kinOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    if (kin_model==2) then
        kinsource(kineq(1)) = modena_outputs_get(kinOutputs, kinSourceNCOPos);
        kinsource(kineq(2)) = modena_outputs_get(kinOutputs, kinSourceOHPos);
        kinsource(kineq(3)) = modena_outputs_get(kinOutputs, kinSourceH2OPos);
        kinsource(kineq(4)) = modena_outputs_get(kinOutputs, kinSourceCO2Pos);
        kinsource(kineq(5)) = modena_outputs_get(kinOutputs, kinSourcePentanePos);
        kinsource(kineq(6)) = modena_outputs_get(kinOutputs, kinSourcePolymerPos);
        kinsource(kineq(7)) = modena_outputs_get(kinOutputs, &
            kinSourcePolymerBlowPos);
        kinsource(kineq(8)) = modena_outputs_get(kinOutputs, kinSourceUreaPos);
        kinsource(kineq(9)) = modena_outputs_get(kinOutputs, kinSourceR1Pos);
        kinsource(kineq(10)) = modena_outputs_get(kinOutputs, kinSourceRmassPos);
        kinsource(kineq(11)) = modena_outputs_get(kinOutputs, kinSourceRvolPos);
        kinsource(kineq(12)) = modena_outputs_get(kinOutputs, kinSourceRtempPos);
    elseif (kin_model==4) then
        kinsource=0
        do i=1,size(kineq)
            kinsource(i) = modena_outputs_get(kinOutputs, kinOutputsPos(i))
            ! write(*,*) i,kinsource(i)
        enddo
        ! stop
    endif
end subroutine kinModel
!***********************************END****************************************


!********************************BEGINNING*************************************
!> prepares integration
subroutine bblpreproc
    integer :: i,j
    write(*,*) 'preparing simulation...'
    !determine number of equations and their indexes
    NEQ=(p+1)*ngas
    NEQ = NEQ+4+ngas
    req=1   !radius index
    fpeq=2  !pressure index
    if (inertial_term) then
        NEQ=NEQ+1
        fpeq=fpeq+1
    endif
    lpeq=fpeq+ngas-1
    teq=lpeq+1   !temperature index
    xOHeq=teq+1 !polyol conversion index
    xWeq=xOHeq+1  !water conversion index
    fceq = xWeq+1    !concentration index
    if (kin_model==2) then
        allocate(kineq(12),kinsource(12))
    elseif (kin_model==4) then
        allocate(kineq(20),kinsource(20))
    endif
    if (kin_model==2 .or. kin_model==4) then
        NEQ=NEQ+size(kineq)
        do i=1,size(kineq)
            kineq(i)=xWeq+i
        enddo
        fceq=kineq(size(kineq))+1
    endif

    !set initial values
    allocate(Y(NEQ))
    Y=0
    Y(req) = R0        !radius
    if (inertial_term) Y(req+1) = 0        !velocity
    Y(teq) = Temp0   !temperature
    Y(xOHeq) = 0        !xOH
    Y(xWeq) = 0        !xW
    if (kin_model==2) then
        Y(kineq(1))=NCO0*1e3_dp
        Y(kineq(2))=OH0*1e3_dp
        Y(kineq(3))=W0*1e3_dp
    endif
    if (kin_model==4) then
        Y(kineq(1)) = 6.73000e-02_dp
        Y(kineq(2)) = 1.92250e+00_dp
        Y(kineq(3)) = 2.26920e+00_dp
        Y(kineq(4)) = 0.00000e+00_dp
        Y(kineq(5)) = 5.46200e-01_dp
        ! Y(kineq(5)) = 1.0924e+00_dp
        Y(kineq(6)) = 2.19790e+00_dp
        Y(kineq(7)) = 1.64000e+00_dp
        Y(kineq(8)) = 1.71030e+00_dp
        Y(kineq(9)) = 0.00000e+00_dp
        Y(kineq(10)) = 0.00000e+00_dp
        Y(kineq(11)) = 0.00000e+00_dp
        Y(kineq(12)) = 0.00000e+00_dp
        Y(kineq(13)) = 0.00000e+00_dp
        Y(kineq(14)) = 0.00000e+00_dp
        Y(kineq(15)) = 0.00000e+00_dp
        Y(kineq(16)) = 4.45849e+00_dp
        Y(kineq(17)) = 0.00000e+00_dp
        Y(kineq(18)) = 1.00000e+00_dp
        Y(kineq(19)) = 60!2.27000e+01_dp
        Y(kineq(20)) = 1e0_dp!8.46382e-01_dp
    endif
    do j=1,ngas
        do i=1,p+1
            Y(fceq+(i-1)*ngas+j-1) = cbl(j)      !blowing agent concentration
        enddo
    enddo
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
        call modena_inputs_set(rhopInputs, rhopTpos, Y(teq))
        call modena_inputs_set(rhopInputs, rhopXOHPos, 0.1_dp)
        ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
        if(ret /= 0) then
            call exit(ret)
        endif
        rhop = modena_outputs_get(rhopOutputs, 0_c_size_t)
        call modena_inputs_set(rhopInputs, rhopTpos, Y(teq)+100)
        call modena_inputs_set(rhopInputs, rhopXOHPos, 0.9_dp)
        ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
        if(ret /= 0) then
            call exit(ret)
        endif
        !average density during foaming
        rhop=(rhop + modena_outputs_get(rhopOutputs, 0_c_size_t))/2
    end select
    call physical_properties(Y)
    Pair0=(Pamb+2*sigma/R0)*xgas(1)
    do i=1,ngas
        Y(fpeq+i-1) = xgas(i+1)*(Pamb+2*sigma/R0) !pressure
        if (Y(fpeq+i-1)<1e-16_dp) Y(fpeq+i-1)=1e-16_dp
    enddo
    Rold=Y(req)
    Told=Y(teq)
    pold(1)=Y(fpeq)
    pold(2)=Y(fpeq+1)
    S0=Sn*Y(req)
    Vsh=4*pi/3*(S0**3-R0**3)
    gelpoint=.false.
    timestep=(TEND-T)/its
    write(*,'(2x,A,2x,e12.6)') 'NN',Sn**(-3)/(1-Sn**(-3))/&
        exp(log(4._dp/3*pi*R0**3))

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

    !choose and set integrator
    select case(integrator)
    case(1)
        select case(MF)
        case(10)
            allocate(RWORK(20+16*NEQ),IWORK(20))
        case(22)
            allocate(RWORK(22+9*NEQ+NEQ**2),IWORK(20+NEQ))
        case default
            stop 'unknown MF'
        end select
    case(2)
        select case(MF)
        case(10)
            allocate(RWORK(20+16*NEQ),IWORK(30))
        case(222)
            NNZ=NEQ**2 !Not sure, smaller numbers make problems for low p
            LENRAT=2 !depends on dp
            allocate(RWORK(int(20+(2+1._dp/LENRAT)*NNZ+(11+9._dp/LENRAT)*NEQ)),&
                IWORK(30))
        case default
            stop 'unknown MF'
        end select
    case default
        stop 'unknown integrator'
    end select
    ITASK = 1
    ISTATE = 1
    IOPT = 1
    RWORK(5:10)=0
    IWORK(5:10)=0
    LRW = size(RWORK)
    LIW = size(IWORK)
    IWORK(6)=maxts
    TOUT =T+timestep
    ITOL = 1 !don't change, or you must declare ATOL as ATOL(NEQ)
    write(*,*) 'done: simulation prepared'
    write(*,*)
end subroutine bblpreproc
!***********************************END****************************************


!********************************BEGINNING*************************************
!> performs integration
subroutine bblinteg(outputs_1d,outputs_GR,outputs_GR_c,outputs_GR_p,concloc)
    character(*),intent(in) :: outputs_1d,outputs_GR,outputs_GR_c,outputs_GR_p,&
        concloc !file names
    write(*,*) 'integrating...'
    call save_integration_header(outputs_1d,outputs_GR,outputs_GR_c,&
        outputs_GR_p,concloc)
    DO IOUT = 1,its
        select case (integrator)
        case(1)
            CALL DLSODE (sub_ptr, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
                ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
        case(2)
            call DLSODES (sub_ptr, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
                ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
        case default
            stop 'unknown integrator'
        end select
        call molar_balance
        call restoreDV
        call save_integration_step
        ! write(*,*) tout,kinsource(2)
        TOUT = TOUT+timestep
        if (eta==maxeta) exit
    END DO
    call save_integration_close
    write(*,*) 'done: integration'
    call destroyModenaModels
    call exit(0)
end subroutine bblinteg
!***********************************END****************************************
end module model
