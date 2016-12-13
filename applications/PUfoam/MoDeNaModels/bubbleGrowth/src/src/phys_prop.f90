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
    integer :: Rb_kx=4,Rb_iknot=0,Rb_inbvx
    real(dp), dimension(:), allocatable :: Rb_tx,Rb_coef
    public set_initial_physical_properties,physical_properties,Rb,Rderiv,&
        Rb_initialized
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
    case(3)
        call modena_inputs_set(rhopInputs, rhopTpos, temp)
        ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
        if(ret /= 0) then
            call exit(ret)
        endif
        rhop = modena_outputs_get(rhopOutputs, 0_c_size_t)
        call modena_inputs_set(rhopInputs, rhopTpos, temp+100)
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
    if (geometry=="3D") then
        pair0=(pamb+2*sigma/radius)*xgas(1)
        Sn=(1._dp/(nb0*4._dp/3*pi*R0**3)+1)**(1._dp/3)
        S0=Sn*radius
        Vsh=4*pi/3*(S0**3-R0**3)
    elseif (geometry=="2D") then
        pair0=(pamb+sigma/radius)*xgas(1)
        Sn=(1._dp/(nb0*pi*R0**2)+1)**(1._dp/2)
        S0=Sn*radius
        Vsh=pi*(S0**2-R0**2)
    endif
    gelpoint=.false.
    timestep=(tend-tstart)/its
    if (firstrun) then
        if (allocated(etat)) deallocate(etat)
        if (allocated(port)) deallocate(port)
        if (allocated(init_bub_rad)) deallocate(init_bub_rad)
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
    if (temp>500) then
        print*, 'temperature > 500', temp
        return
    endif
    if (.not. gelpoint) then
        select case(visc_model)
        case(1)
        case(2)
            Aeta=4.1e-8_dp
            Eeta=38.3e3_dp
            AA=4.0_dp
            B=-2.0_dp
            eta=Aeta*exp(Eeta/(Rg*temp))*&
                (gelpointconv/(gelpointconv-conv))**(AA+B*conv)
            if (eta>1e10) eta=1.0e10_dp
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
        if (conv>gelpointconv*0.9_dp) then
            gelpoint=.true.
        endif
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
        case(2)
            call modena_inputs_set(diffInputs(i), diffTpos(i), temp)
            ret = modena_model_call(diffModena(i),diffInputs(i),diffOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            D(i) = modena_outputs_get(diffOutputs(i), 0_c_size_t)
        end select
        select case(sol_model(i))
        case(1) !constant Baser 10.1002/pen.760340805
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            ret = modena_model_call(solModena(i), solInputs(i), solOutputs(i))
            if(ret /= 0) then
                call exit(ret)
            endif
            KH(i) = modena_outputs_get(solOutputs(i), 0_c_size_t)
        case(2) !pcsaft
            ! TODO: implement properly
            call modena_inputs_set(solInputs(i), solTpos(i), temp)
            ! call modena_inputs_set(solInputs(i), solXgasPos(i), 1.0e-4_dp)
            ! call modena_inputs_set(solInputs(i), solXmdiPos(i), 0.5_dp)
            ! call modena_inputs_set(solInputs(i), solXpolyolPos(i), 0.5_dp)
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
        case(4) !pentane, Winkler Ph.D.
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
        case(7) !pentane, Winkler Ph.D.
            KH(i)=rhop/Mbl(i)/pamb*(0.0064_dp+0.0551_dp*exp(-(temp-298)**2/&
                (2*17.8_dp**2)))
        case(8) !constant Baser 10.1002/pen.760340805
        end select
    enddo
    if (solcorr) then
        if (geometry=="3D") then
            KH=KH*exp(2*sigma*Mbl/(rhop*Rg*temp*radius))
        elseif (geometry=="2D") then
            KH=KH*exp(sigma*Mbl/(rhop*Rg*temp*radius))
        endif
    endif
    cp=cppol+sum(cbl*Mbl*cpbll)/rhop
end subroutine physical_properties
!***********************************END****************************************


!********************************BEGINNING*************************************
!> time derivation of bubble radius as function of time
real(dp) function Rderiv(t)
    use bspline_module
    real(dp) :: t
    real(dp) :: dt
    integer :: nx,idx,iflag
    ! dt=timestep/100
    ! Rderiv=(Rb(t+dt)-Rb(t))/dt
    nx=size(bub_rad(:,1))
    if (.not. Rb_initialized) then
        call Rb_spline_ini
    endif
    idx=1
    call db1val(t,idx,Rb_tx,nx,Rb_kx,Rb_coef,Rderiv,iflag,Rb_inbvx)
    if (iflag /= 0) then
        print*, 'evaluation of bubble radius derivative from spline failed',&
            iflag
        stop
    endif
endfunction Rderiv
!***********************************END****************************************


!********************************BEGINNING*************************************
!> bubble radius as function of time
real(dp) function Rb(t)
    use bspline_module
    use interpolation
    real(dp) :: t
    integer :: nx,idx,iflag
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    nx=size(bub_rad(:,1))
    ! xi(1)=t
    ! call pwl_interp_1d ( nx, bub_rad(:,1), bub_rad(:,bub_inx+1), ni, xi, yi )
    ! Rb=yi(1)
    if (.not. Rb_initialized) then
        call Rb_spline_ini
    endif
    idx=0
    call db1val(t,idx,Rb_tx,nx,Rb_kx,Rb_coef,Rb,iflag,Rb_inbvx)
    if (iflag /= 0) then
        print*, 'evaluation of bubble radius from spline failed',iflag
        stop
    endif
endfunction Rb
!***********************************END****************************************


!********************************BEGINNING*************************************
!> initialization of spline
subroutine Rb_spline_ini
    use bspline_module
    integer :: nx,iflag
    nx=size(bub_rad(:,1))
    if (allocated(Rb_tx)) deallocate(Rb_tx)
    if (allocated(Rb_coef)) deallocate(Rb_coef)
    allocate(Rb_tx(nx+Rb_kx),Rb_coef(nx))
    Rb_tx=0
    call db1ink(bub_rad(:,1),nx,bub_rad(:,bub_inx+1),&
        Rb_kx,Rb_iknot,Rb_tx,Rb_coef,iflag)
    Rb_inbvx=1
    Rb_initialized=.true.
    if (iflag /= 0) then
        print*, 'initialization of spline failed'
        stop
    endif
end subroutine Rb_spline_ini
!***********************************END****************************************
end module phys_prop
