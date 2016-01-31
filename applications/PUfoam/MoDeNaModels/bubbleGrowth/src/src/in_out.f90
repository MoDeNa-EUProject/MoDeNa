!> @file
!! handles input and output
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module in_out
    use constants
    use ioutils, only:newunit,str
    implicit none
    logical :: inertial_term,solcorr,gelpoint,dilution
    integer :: fi1,fi2,fi3,fi4,&
        integrator,p,maxts,its,visc_model,rhop_model,itens_model,ngas,co2_pos,&
        kin_model,MF,NEQ,&
        fceq,& !first concentration equation (index)
        fpeq,lpeq,& !first and last pressure equation (index)
        req,& !radius equation (index)
        teq,& !temperature equation (index)
        xOHeq,xWeq !conversion equations (indexes)
    real(dp) :: mshco,& !mesh coarsening parameter
        Temp0,R0,Sn,OH0,W0,NCO0,AOH,EOH,AW,EW,dHOH,dHW,&
        time,radius,eqconc,grrate(2),st,S0,&
        T,TEND,RTOL,ATOL,&
        eta,maxeta,Aeta,Eeta,Cg,AA,B,&
        Pamb,sigma,rhop,cp,cppol,rhobl
    integer, dimension(:), allocatable :: diff_model,sol_model,fic,&
        kineq !kinetics state variable equations (indexes)
    real(dp), dimension(:), allocatable :: Y,cbl,xgas,&
        kinsource,& !kinetic source term
        D,KH,Mbl,dHv,cpblg,cpbll,&
        mb,mb2,mb3,avconc,pressure
contains
!********************************BEGINNING*************************************
!> reads input values from a file
subroutine read_inputs(inputs)
    character(len=80),intent(in) :: inputs
    integer :: fi
    write(*,*) 'loading input file ',TRIM(inputs)
    open(newunit(fi),file=inputs)
        read(fi,*) integrator !integrator. 1=dlsode,2=dlsodes
        read(fi,*) MF   !10=nonstiff,22=stiff,automatic Jacobian(dlsode),
            ! 222=stiff,automatic Jacobian(dlsodes)
        read(fi,*)
        read(fi,*) inertial_term    !include inertial term in equations (t/f)
        read(fi,*) solcorr  !use solubility correction on bubble radius (t/f)
        read(fi,*) mshco    !mesh coarsening parameter
        read(fi,*)
        read(fi,*) p    !number of internal nodes
        read(fi,*) T    !initial time
        read(fi,*) TEND    !final time
        read(fi,*) its    !number of outer integration time steps (how many
            ! times are values written)
        read(fi,*) maxts    !maximum inner time steps between t and t+h
            ! (default 500, recommended 50000)
        read(fi,*) RTOL    !relative tolerance
        read(fi,*) ATOL    !absolute tolerance
        read(fi,*)
        read(fi,*) ngas     !number of dissolved gases
        allocate(D(ngas),cbl(ngas),xgas(ngas+1),KH(ngas),fic(ngas),Mbl(ngas),&
            dHv(ngas),mb(ngas),mb2(ngas),mb3(ngas),avconc(ngas),pressure(ngas),&
            diff_model(ngas),sol_model(ngas),cpblg(ngas),cpbll(ngas))
        read(fi,*) co2_pos     !carbon dioxide position
        read(fi,*) Pamb    !ambient pressure
        read(fi,*) Mbl    !blowing agent molar mass (for each dissolved gas)
        read(fi,*) cppol    !heat capacity of polymer
        read(fi,*) cpbll    !heat capacity of blowing agent in liquid phase
            ! (for each)
        read(fi,*) cpblg    !heat capacity of blowing agent in gas phase
            ! (for each)
        read(fi,*) dHv    !evaporation heat of blowing agent (for each gas)
        read(fi,*) rhobl    !density of liquid physical blowing agent
        read(fi,*)
        read(fi,*) Temp0    !initial temperature
        read(fi,*) R0    !initial radius
        read(fi,*) Sn    !how many times is initial shell larger than initial
            ! bubble radius
        read(fi,*) OH0    !initial concentration of polyol (don't set to zero -
            ! division by zero; if you don't want reaction, set water to zero)
        read(fi,*) W0    !initial concentration of water (if you set this to
            ! zero, water conversion results are meanigless)
        read(fi,*) NCO0    !initial concentration of isocyanate
        read(fi,*) cbl    !initial concentration of disolved blowing agent
            ! (for each dissolved gas)
        read(fi,*) xgas    !initial molar fraction of gases in the bubble (for
            ! air and each dissolved gas)
        read(fi,*)
        read(fi,*) kin_model   !reaction kinetics model. 1=Baser,2=modena simple
            ! kinetics, 3=Baser with R(x), 4=modena RF-1-private
        read(fi,*) dilution   !use dilution effect
        read(fi,*) AOH    !frequential factor of gelling reaction
        read(fi,*) EOH    !activation energy of gelling reaction
        read(fi,*) AW    !frequential factor of blowing reaction
        read(fi,*) EW    !activation energy of blowing reaction
        read(fi,*) dHOH    !gelling reaction enthalpy
        read(fi,*) dHW    !blowing reaction enthalpy
        read(fi,*)
        read(fi,*) rhop_model   !polymer density model. 1=constant,2=modena
        read(fi,*) rhop    !polymer density
        read(fi,*)
        read(fi,*) itens_model  !interfacial tension model. 1=constant,2=modena
        read(fi,*) sigma    !interfacial tension
        read(fi,*)
        read(fi,*) diff_model  !diffusivity model (for each dissolved gas).
            ! 1=constant,2=modena
        read(fi,*) D    !diffusion coefficients (for each dissolved gas)
        read(fi,*)
        read(fi,*) sol_model   !solubility model (for each dissolved gas).
            ! 1=constant,2=modena
        read(fi,*) KH    !Henry constants (for each dissolved gas)
        read(fi,*)
        read(fi,*) visc_model    !viscosity model. 1=constant,2=Castro and
            ! Macosko,3=modena
        read(fi,*) eta    !viscosity (if constant viscosity is used)
        read(fi,*) maxeta    !maximum viscosity
        read(fi,*) Aeta    !viscosity constant Aeta
        read(fi,*) Eeta    !viscosity constant Eeta
        read(fi,*) Cg    !viscosity constant Cg
        read(fi,*) AA    !viscosity constant AA
        read(fi,*) B    !viscosity constant B
    close(fi)
    write(*,*) 'done: inputs loaded'
    write(*,*)
end subroutine read_inputs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves parameters of surrogate model
subroutine save_surrogate_parameters(spar)
    character(len=80),intent(in) :: spar
    integer :: fi,i
    write(*,*) 'saving parameters of surrogate model to ',TRIM(spar)
    open(newunit(fi),file=spar)
        write(fi,*) ngas    !number of dissolved gases
        write(fi,*) sigma    !interfacial tension
        do i=1,ngas
        	write(fi,*) KH(i)    !Henry constants (for each dissolved gas)
        enddo
    close(fi)
    write(*,*) 'done: parameters of surrogate model saved'
    write(*,*)
end subroutine save_surrogate_parameters
!***********************************END****************************************


!********************************BEGINNING*************************************
!> opens output files and writes a header
subroutine save_integration_header(outputs_1d,outputs_GR,outputs_GR_c,&
    outputs_GR_p,concloc)
    character(*),intent(in) :: outputs_1d,outputs_GR,outputs_GR_c,outputs_GR_p,&
        concloc !file names
    integer :: i
    open (unit=newunit(fi1), file = outputs_1d)
    write(fi1,'(1000A23)') '#time', 'radius','pressure','conversion of polyol',&
        'conversion of water', 'eq. concentration', 'first concentration', &
        'viscosity', 'moles in polymer', 'moles in bubble', 'total moles', &
        'shell thickness', 'temperature', 'foam density', 'weight fraction'
    open (unit=newunit(fi2), file = outputs_GR)
    write(fi2,'(1000A23)') '#GrowthRate1', 'GrowthRate2', 'temperature', &
        'bubbleRadius', 'KH1','KH2','c1','c2','p1','p2'
    open (unit=newunit(fi3), file = '../results/kinetics.out')
    write(fi3,'(1000A23)') "time","Catalyst_1","CE_A0","CE_A1","CE_B","CE_B2",&
        "CE_I0","CE_I1","CE_I2","CE_PBA","CE_Breac","CE_Areac0","CE_Areac1",&
        "CE_Ireac0","CE_Ireac1","CE_Ireac2","Bulk","R_1","R_1_mass","R_1_temp",&
        "R_1_vol"
    open (unit=newunit(fi4), file = outputs_GR_c)
end subroutine save_integration_header
!***********************************END****************************************


!********************************BEGINNING*************************************
!> writes an integration step to output file
subroutine save_integration_step
    integer :: i
    write(fi1,"(1000es23.15)") time,radius, pressure(1), Y(xOHeq), Y(xWeq), &
        eqconc,Y(fceq),eta,mb(1),mb2(1),mb3(1),st,Y(teq),(1-radius**3/&
        (radius**3+S0**3-R0**3))*rhop,mb(2)*Mbl(2)/(rhop*4*pi/3*(S0**3-R0**3)),&
        mb(1)*Mbl(1)/(rhop*4*pi/3*(S0**3-R0**3))
    write(fi2,"(1000es23.15)") grrate, Y(teq), radius, KH, avconc, pressure
    ! write(fi3,"(1000es23.15)") time,Y(kineq(1)),Y(kineq(2)),Y(kineq(3)),&
    !     Y(kineq(4)),Y(kineq(5)),Y(kineq(6)),Y(kineq(7)),Y(kineq(8)),&
    !     Y(kineq(9)),Y(kineq(10)),Y(kineq(11)),Y(kineq(12)),Y(kineq(13)),&
    !     Y(kineq(14)),Y(kineq(15)),Y(kineq(16)),Y(kineq(17)),Y(kineq(18)),&
    !     Y(kineq(19)),Y(kineq(20))
    write(fi4,"(1000es23.15)") (Y(fceq+i+1),i=0,ngas*p,ngas)
end subroutine save_integration_step
!***********************************END****************************************


!********************************BEGINNING*************************************
!> closes output files
subroutine save_integration_close
!***************************DECLARATION******************************
    integer :: i
!******************************BODY**********************************
    close(fi1)
    close(fi2)
    close(fi3)
    close(fi4)
end subroutine save_integration_close
!***********************************END****************************************
end module in_out
