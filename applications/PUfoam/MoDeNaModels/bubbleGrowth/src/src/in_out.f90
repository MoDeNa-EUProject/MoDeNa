!> @file
!! handles input and output
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module in_out
    use foaming_globals_m
    use constants
    use ioutils, only:newunit,str
    implicit none
    character(len=99) :: &
        fileplacein,& ! location of input files
        fileplaceout,& ! location of output files
        inputs,& ! input file
        outputs_1d,& ! output file with scalar variables
        outputs_GR,& ! output file for the surrogate model fitting
        outputs_c,& ! output file with concentration profiles
        outputs_kin,& ! output file with variables of detailed kinetic model
        outputs_drdt ! test output file for box with multiple growing bubbles
    logical :: inertial_term,solcorr,gelpoint,dilution
    integer :: fi1,fi2,fi3,fi4,fi5,&
        integrator,p,maxts,its,visc_model,rhop_model,itens_model,ngas,co2_pos,&
        kin_model,int_meth,&
        fceq,& !first concentration equation (index)
        fpeq,lpeq,& !first and last pressure equation (index)
        req,& !radius equation (index)
        teq,& !temperature equation (index)
        xOHeq,xWeq !conversion equations (indexes)
    real(dp) :: mshco,& !mesh coarsening parameter
        temp0,R0,Sn,OH0,W0,NCO0,AOH,EOH,AW,EW,dHOH,dHW,&
        time,radius,eqconc,grrate(2),st,S0,&
        rel_tol,abs_tol,&
        eta,maxeta,Aeta,Eeta,Cg,AA,B,&
        pamb,sigma,rhop,cp,cppol,rhobl,porosity,rhofoam,&
        pair0,pair,timestep,gr,nold(2),vsh,temp,conv,&
        Rey,pairst,pambst,Ca
    integer, dimension(:), allocatable :: diff_model,sol_model,fic,&
        kineq !kinetics state variable equations (indexes)
    real(dp), dimension(:), allocatable :: y,cbl,xgas,&
        kinsource,& !kinetic source term
        D,D0,KH,Mbl,dHv,cpblg,cpbll,&
        mb,mb2,mb3,avconc,pressure,times,dRdt,Rt,pt,ATOL2,wblpol
contains
!********************************BEGINNING*************************************
!> set paths to all files
subroutine set_paths
    fileplacein='../'
    fileplaceout='../results/'
    inputs='inputs.in'
    outputs_1d='outputs_1d.out'
    outputs_GR='outputs_GR.out'
    outputs_c='outputs_c.out'
    outputs_kin='kinetics.out'
    outputs_drdt='dRdt.out'
    inputs=TRIM(ADJUSTL(fileplacein))//TRIM(ADJUSTL(inputs))
    outputs_1d=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_1d))
    outputs_GR=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR))
    outputs_c=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_c))
    outputs_kin=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_kin))
    outputs_drdt=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_drdt))
end subroutine set_paths
!***********************************END****************************************


!********************************BEGINNING*************************************
!> reads input values from a file
subroutine read_inputs
    integer :: fi
    write(*,*) 'loading input file ',TRIM(inputs)
    open(newunit(fi),file=inputs)
        read(fi,*) integrator !integrator. 1=dlsode,2=dlsodes
        read(fi,*) int_meth   !10=nonstiff,22=stiff,automatic Jacobian(dlsode),
            ! 222=stiff,automatic Jacobian(dlsodes)
        read(fi,*)
        read(fi,*) inertial_term    !include inertial term in equations (t/f)
        read(fi,*) solcorr  !use solubility correction on bubble radius (t/f)
        read(fi,*) mshco    !mesh coarsening parameter
        read(fi,*)
        read(fi,*) p    !number of internal nodes
        read(fi,*) tstart    !initial time
        if (firstrun) then
            read(fi,*) tend    !final time
        else
            read(fi,*) ! final time is already set
        endif
        read(fi,*) its    !number of outer integration time steps (how many
            ! times are values written)
        read(fi,*) maxts    !maximum inner time steps between t and t+h
            ! (default 500)
        read(fi,*) rel_tol    !relative tolerance
        read(fi,*) abs_tol    !absolute tolerance
        read(fi,*)
        read(fi,*) ngas     !number of dissolved gases
        allocate(D(ngas),cbl(ngas),xgas(ngas+1),KH(ngas),fic(ngas),Mbl(ngas),&
            dHv(ngas),mb(ngas),mb2(ngas),mb3(ngas),avconc(ngas),pressure(ngas),&
            diff_model(ngas),sol_model(ngas),cpblg(ngas),cpbll(ngas),&
            wblpol(ngas),D0(ngas))
        read(fi,*) co2_pos     !carbon dioxide position
        read(fi,*) pamb    !ambient pressure
        read(fi,*) Mbl    !blowing agent molar mass (for each dissolved gas)
        read(fi,*) cppol    !heat capacity of polymer
        read(fi,*) cpbll    !heat capacity of blowing agent in liquid phase
            ! (for each)
        read(fi,*) cpblg    !heat capacity of blowing agent in gas phase
            ! (for each)
        read(fi,*) dHv    !evaporation heat of blowing agent (for each gas)
        read(fi,*) rhobl    !density of liquid physical blowing agent
        read(fi,*)
        read(fi,*) temp0    !initial temperature
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
        read(fi,*) kin_model   !reaction kinetics model. 1=Baser,
            ! 3=Baser with R(x), 4=modena RF-1-private
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
!> opens output files and writes a header
subroutine save_integration_header
    integer :: i
    open (unit=newunit(fi1), file = outputs_1d)
    write(fi1,'(1000A23)') '#time', 'radius','pressure1', 'pressure2',&
        'conversion_of_polyol',&
        'conversion_of_water', 'eq.concentration', 'first_concentration', &
        'viscosity', 'moles_in_polymer', 'moles_in_bubble', 'total_moles', &
        'shell_thickness', 'temperature', 'foam_density', 'weight_fraction1', &
        'weight_fraction2','porosity'
    open (unit=newunit(fi2), file = outputs_GR)
    write(fi2,'(1000A23)') '#GrowthRate1', 'GrowthRate2', 'temperature', &
        'bubbleRadius', 'KH1','KH2','c1','c2','p1','p2'
    if (kin_model==4) then
        open (unit=newunit(fi3), file = outputs_kin)
        write(fi3,'(1000A23)') "time","Catalyst_1","CE_A0","CE_A1","CE_B",&
            "CE_B2","CE_I0","CE_I1","CE_I2","CE_PBA","CE_Breac","CE_Areac0",&
            "CE_Areac1","CE_Ireac0","CE_Ireac1","CE_Ireac2","Bulk","R_1",&
            "R_1_mass","R_1_temp","R_1_vol"
    endif
    open (unit=newunit(fi4), file = outputs_c)
    open (unit=newunit(fi5), file = outputs_drdt)
end subroutine save_integration_header
!***********************************END****************************************


!********************************BEGINNING*************************************
!> writes an integration step to output file
subroutine save_integration_step(iout)
    integer :: i,iout
    real(dp) :: rder
    rder=0
    write(fi1,"(1000es23.15)") time,radius,pressure,Y(xOHeq),Y(xWeq),&
        eqconc,Y(fceq),eta,mb(1),mb2(1),mb3(1),st,temp,rhofoam,&
        wblpol,porosity
    write(fi2,"(1000es23.15)") grrate, temp, radius, KH, avconc, pressure
    if (kin_model==4) then
        write(fi3,"(1000es23.15)") time,Y(kineq(1)),Y(kineq(2)),Y(kineq(3)),&
            Y(kineq(4)),Y(kineq(5)),Y(kineq(6)),Y(kineq(7)),Y(kineq(8)),&
            Y(kineq(9)),Y(kineq(10)),Y(kineq(11)),Y(kineq(12)),Y(kineq(13)),&
            Y(kineq(14)),Y(kineq(15)),Y(kineq(16)),Y(kineq(17)),Y(kineq(18)),&
            Y(kineq(19)),Y(kineq(20))
    endif
    write(fi4,"(1000es23.15)") (Y(fceq+i+1),i=0,ngas*p,ngas)
    write(fi5,"(1000es23.15)") time,radius,rder,pressure(1)
    ! save arrays, which are preserved for future use
    etat(iout,1)=time
    etat(iout,2)=eta
    port(iout,1)=time
    port(iout,2)=porosity
    init_bub_rad(iout,1)=time
    init_bub_rad(iout,2)=radius
end subroutine save_integration_step
!***********************************END****************************************


!********************************BEGINNING*************************************
!> closes output files
subroutine save_integration_close(iout)
    integer :: iout
    real(dp), dimension(:,:), allocatable :: matr
    close(fi1)
    close(fi2)
    close(fi3)
    if (kin_model==4) then
        close(fi4)
    endif
    close(fi5)
    ! reallocate matrices for eta_rm and bub_vf functions
    ! interpolation doesn't work otherwise
    if (iout /= its) then
        allocate(matr(0:its,2))
        matr=etat
        deallocate(etat)
        allocate(etat(0:iout,2))
        etat=matr(0:iout,:)
        matr=port
        deallocate(port)
        allocate(port(0:iout,2))
        port=matr(0:iout,:)
        matr=init_bub_rad
        deallocate(init_bub_rad)
        allocate(init_bub_rad(0:iout,2))
        init_bub_rad=matr(0:iout,:)
    endif
end subroutine save_integration_close
!***********************************END****************************************


!********************************BEGINNING*************************************
!> loads old results
subroutine load_old_results
    integer :: i,j,ios
    real(dp), dimension(:,:), allocatable :: matrix
    j=0
    open(newunit(fi5),file=outputs_1d)
        do  !find number of points
            read(fi5,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        allocate(matrix(j,18))
        rewind(fi5)
        read(fi5,*)
        do i=2,j
            read(fi5,*) matrix(i,:)
        enddo
    close(fi5)
    allocate(times(j-1),Rt(j-1))
    times(1:j-1)=matrix(2:j,1)
    Rt(1:j-1)=matrix(2:j,2)
    deallocate(matrix)
end subroutine load_old_results
!***********************************END****************************************
end module in_out
