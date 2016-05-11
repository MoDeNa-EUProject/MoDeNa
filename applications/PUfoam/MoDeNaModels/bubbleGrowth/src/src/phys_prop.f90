!> @file
!! calculates material properties of the system
!! results are stored in global variables
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module phys_prop
    use constants
    use globals
    use foaming_globals_m
    use fmodena
    use modenastuff
    implicit none
    private
    ! interpolation variables
    logical :: Rb_initialized
    integer :: Rb_kx=5,Rb_iknot=0,Rb_inbvx
    real(dp), dimension(:), allocatable :: Rb_tx,Rb_coef
    public set_initial_physical_properties,physical_properties,Rb,Rderiv
contains
!********************************BEGINNING*************************************
!> determine physical properties
subroutine set_initial_physical_properties
    time=tstart
    if (.not. firstrun) R0=Rb(time)
    radius=R0
    temp=temp0
    conv=0.0_dp
    ! initial bubble contains only air
    xgas=0
    xgas(1)=1
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
    call physical_properties(temp,conv,radius)
    D0=D
    surface_tension=sigma
    pair0=(pamb+2*sigma/radius)*xgas(1)
    Sn=(1._dp/(nb0*exp(log(4._dp/3*pi*R0**3)))+1)**(1._dp/3)
    ! write(*,'(2x,A,2x,e12.6)') 'NN',Sn**(-3)/(1-Sn**(-3))/&
    !     exp(log(4._dp/3*pi*R0**3))
    S0=Sn*radius
    Vsh=4*pi/3*(S0**3-R0**3)
    gelpoint=.false.
    timestep=(tend-tstart)/its
    if (firstrun) then
        allocate(etat(0:its,2),port(0:its,2),init_bub_rad(0:its,2))
    endif
end subroutine set_initial_physical_properties
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates values of physical properties
subroutine physical_properties(temp,conv,radius)
    real(dp), intent(in) :: temp,conv,radius
    integer :: i
    real(dp) :: Aeta,Eeta,Cg,AA,B !viscosity model constants
    if (.not. gelpoint .and. temp<500) then
        select case(visc_model)
        case(1)
        case(2)
            Aeta=4.1e-8_dp
            Eeta=38.3e3_dp
            Cg=0.85_dp
            AA=4.e0_dp
            B=-2.e0_dp
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
            write(*,'(2x,A,es8.2,A)') 'gel point reached at time t = ',time,' s'
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
            ret = modena_model_call(diffModena(i),diffInputs(i),diffOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            D(i) = modena_outputs_get(diffOutputs(i), 0_c_size_t)
        end select
        select case(sol_model(i))
        case(1)
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            ret = modena_model_call(solModena(i), solInputs(i), solOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            KH(i) = modena_outputs_get(solOutputs(i), 0_c_size_t)
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
            ! KH(i)=rhop/Mbl(i)/KH(i)
        case(3) !n-pentane, Gupta, 10.1002/pen.11405
            ! KH(i)=rhop/Mbl(i)/pamb*3.3e-4_dp*(exp((2.09e4_dp-67.5_dp*(temp-&
            !     35.8_dp*log(pamb/1e5_dp)))/(8.68e4_dp-(temp-35.8_dp*&
            !     log(pamb/1e5_dp))))-1.01_dp)**(-1)
            ! KH(i)=-rhop/Mbl(i)/pamb*3.3e-4_dp/(exp((2.09e4_dp-67.5_dp*temp)/&
            !     (8.69e4_dp-temp))-1.01_dp)
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            ret = modena_model_call(solModena(i), solInputs(i), solOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            KH(i) = modena_outputs_get(solOutputs(i), 0_c_size_t)
        case(4) !n-pentane, Winkler Ph.D.
            ! KH(i)=rhop/Mbl(i)/pamb*(0.0064_dp+0.0551_dp*exp(-(temp-298)**2/&
            !     (2*17.8_dp**2)))
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            ret = modena_model_call(solModena(i), solInputs(i), solOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            KH(i) = modena_outputs_get(solOutputs(i), 0_c_size_t)
        case(6) !R11, Baser, 10.1002/pen.760340804
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            ret = modena_model_call(solModena(i), solInputs(i), solOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            KH(i) = modena_outputs_get(solOutputs(i), 0_c_size_t)
            ! print*, KH(i),rhop,Mbl(i),pamb
            ! KH(i)=rhop/Mbl(i)/pamb*(1e-7_dp+4.2934_dp*&
            !     exp(-(temp-203.3556_dp)**2/(2*40.016_dp**2)))
            ! print*, KH(i)
            ! stop
        end select
    enddo
    if (solcorr) KH=KH*exp(2*sigma*Mbl/(rhop*Rg*temp*radius))
    cp=cppol+sum(cbl*Mbl*cpbll)/rhop
end subroutine physical_properties
!***********************************END****************************************


!********************************BEGINNING*************************************
!> time derivation of bubble radius as function of time
real(dp) function Rderiv(t)
    use bspline_module
    real(dp) :: t
    integer :: nx,idx,iflag
    nx=size(bub_rad(:,1))
    if (.not. Rb_initialized) then
        allocate(Rb_tx(nx+Rb_kx),Rb_coef(nx))
        call db1ink(bub_rad(:,1),nx,bub_rad(:,bub_inx+1),&
            Rb_kx,Rb_iknot,Rb_tx,Rb_coef,iflag)
        Rb_inbvx=1
        Rb_initialized=.true.
    endif
    idx=1
    call db1val(t,idx,Rb_tx,nx,Rb_kx,Rb_coef,Rderiv,iflag,Rb_inbvx)
endfunction Rderiv
!***********************************END****************************************


!********************************BEGINNING*************************************
!> bubble radius as function of time
real(dp) function Rb(t)
    use bspline_module
    real(dp) :: t
    integer :: nx,idx,iflag
    nx=size(bub_rad(:,1))
    if (.not. Rb_initialized) then
        allocate(Rb_tx(nx+Rb_kx),Rb_coef(nx))
        call db1ink(bub_rad(:,1),nx,bub_rad(:,bub_inx+1),&
            Rb_kx,Rb_iknot,Rb_tx,Rb_coef,iflag)
        Rb_inbvx=1
        Rb_initialized=.true.
    endif
    idx=0
    call db1val(t,idx,Rb_tx,nx,Rb_kx,Rb_coef,Rb,iflag,Rb_inbvx)
endfunction Rb
!***********************************END****************************************
end module phys_prop
