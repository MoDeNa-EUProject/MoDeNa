!> @file      foamConductivity/src/src/tests.f90
!! @ingroup   src_mod_foamConductivity
!! @author    Pavel Ferkl
!! @brief     Top level routines.
!! @details
!! Main subroutines and test functions.
module tests
    use constants
    use ioutils
    use foamprop, only: effrad
    use foamgeom, only: foam_morpholgy
    use conduction, only: effcond
    use condrad, only: equcond,cond,alpha,sigma
    implicit none
    private
    public loadParameters,eqcond,eqcond_por,eqcond_dcell,eqcond_strut,&
        eqcond_dfoam
    character(len=99) :: fileplacein_par='./inputs/'   !modena
    character(len=99) :: fileplacein_ref='../spectra/'  !modena
    character(len=99) :: fileplaceout='./'  !modena
    character(len=99) :: inputs='foamConductivity.json',spectra='spectra.out'
    character(len=99) :: after_foaming,after_foaming0='after_foaming.txt'
    character(len=99) :: bg_res='../../foamExpansion/results/bubbleGrowth/'
    character(len=99) :: qmom0D_res='../../foamExpansion/results/CFD0D/'
    character(len=99) :: qmom3D_res='../../foamExpansion/results/CFD3D/'
contains
!********************************BEGINNING*************************************
!> Calculates equivalent conductivity for one specific foam.
subroutine eqcond(regions)
    !> number of regions, on which the foam is divided.
    integer, intent(in) :: regions
    integer :: i,j,fi
    real(dp), dimension(:), allocatable :: regbound,regcond
    real(dp), dimension(:,:), allocatable :: regalpha,regsigma
    allocate(regbound(regions+1),regcond(regions),regalpha(regions,nbox),&
        regsigma(regions,nbox),alpha(nz,nbox),sigma(nz,nbox),cond(nz))
    do i=1,regions+1
        regbound(i)=(i-1)*dfoam/regions
    enddo
    do i=1,regions
        if (por < 0.8_dp) then
            write(*,*) 'Radiative prop. not calculated for porosity < 0.8'
            write(mfi,*) 'Radiative prop. not calculated for porosity < 0.8'
            krad=2e-3_dp
            call effcond
            open(newunit(fi),file='foamConductivity.out')
            write(fi,*) eqc_ross
            close(fi)
            return
        endif
        call foam_morpholgy
        if (testMode) then
            write(*,*) 'TESTING: radiative properties not calculated.'
            write(*,*) 'Disable test mode if you want more reasonable results.'
            krad=2e-3_dp
            call effcond
            open(newunit(fi),file='foamConductivity.out')
            write(fi,*) eqc_ross
            close(fi)
            return
        endif
        call effrad(spectra)
        call effcond
        regcond(i)=effc
        regalpha(i,:)=abscoeffbox
        regsigma(i,:)=scattcoeffbox
        deallocate(abscoeffbox,scattcoeffbox)
    enddo
    do i=1,nz
        do j=1,regions+1
            if ((i-0.5_dp)*dfoam/nz<regbound(j)) exit
        enddo
        j=j-1
        cond(i)=regcond(j)
        alpha(i,:)=regalpha(j,:)
        sigma(i,:)=regsigma(j,:)
    enddo
    call equcond
    deallocate(regbound,regcond,regalpha,regsigma,alpha,sigma,cond)
end subroutine eqcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculate dependance of equivalent conductivity on porosity.
subroutine eqcond_por
    integer :: fi,npoints,i
    real(dp) :: pormin,pormax,dpor
    pormin=0.90_dp
    pormax=0.99_dp
    npoints=19
    dpor=(pormax-pormin)/(npoints-1)
    open(newunit(fi),file='eqcond_por.csv')
    write(fi,*) 'porosity,foam_density,eq_cond,Ross_eq_cond,&
        kgas,ksol,krad'
    do i=1,npoints
        por=pormin+(i-1)*dpor
        rhof=(1-por)*rhos
        call eqcond(1)
        write(fi,'(1x,es23.15,7(",",es23.15))') &
            por,rhof,eqc,eqc_ross,kgas,ksol,krad
        write(*,*)
        write(mfi,*)
    enddo
    close(fi)
end subroutine eqcond_por
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculate dependance of equivalent conductivity on cell size.
subroutine eqcond_dcell
    integer :: fi,npoints,i
    real(dp) :: dcellmin,dcellmax,ddcell
    dcellmin=1e-6_dp
    dcellmax=1000e-6_dp
    npoints=7
    ddcell=log10(dcellmax/dcellmin)/(npoints-1)
    open(newunit(fi),file='eqcond_dcell.csv')
    write(fi,*) 'porosity,foam_density,eq_cond,Ross_eq_cond,&
        kgas,ksol,krad,dcell,fstrut'
    do i=1,npoints
        dcell=dcellmin*10**((i-1)*ddcell)
        call eqcond(1)
        write(fi,'(1x,es23.15,9(",",es23.15))') &
            por,rhof,eqc,eqc_ross,kgas,ksol,krad,dcell,fs
    enddo
    close(fi)
end subroutine eqcond_dcell
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculate dependance of equivalent conductivity on strut content.
subroutine eqcond_strut
    integer :: fi,npoints,i
    real(dp) :: strutmin,strutmax,dstrut
    strutmin=0.3_dp
    strutmax=0.9_dp
    npoints=19
    dstrut=(strutmax-strutmin)/(npoints-1)
    open(newunit(fi),file='eqcond_strut.csv')
    write(fi,*) 'porosity,foam_density,eq_cond,Ross_eq_cond,&
        kgas,ksol,krad,dcell,fstrut'
    do i=1,npoints
        fs=strutmin+(i-1)*dstrut
        call eqcond(1)
        write(fi,'(1x,es23.15,9(",",es23.15))') &
            por,rhof,eqc,eqc_ross,kgas,ksol,krad,dcell,fs
    enddo
    close(fi)
end subroutine eqcond_strut
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculate dependance of equivalent conductivity on cell size.
subroutine eqcond_dfoam
    integer :: fi,npoints,i
    real(dp) :: dfoammin,dfoammax,ddfoam
    dfoammin=1e-3_dp
    dfoammax=1000e-3_dp
    npoints=7
    ddfoam=log10(dfoammax/dfoammin)/(npoints-1)
    open(newunit(fi),file='eqcond_dfoam.csv')
    write(fi,*) 'porosity,foam_density,eq_cond,Ross_eq_cond,&
        kgas,ksol,krad,dcell,fstrut,dfoam'
    do i=1,npoints
        dfoam=dfoammin*10**((i-1)*ddfoam)
        call eqcond(1)
        write(fi,'(1x,es23.15,9(",",es23.15))') &
            por,rhof,eqc,eqc_ross,kgas,ksol,krad,dcell,fs,dfoam
    enddo
    close(fi)
end subroutine eqcond_dfoam
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Loads parameters, usually from text file.
subroutine loadParameters
    use physicalProperties
    use fson
    use fson_value_m, only: fson_value_get
    type(fson_value), pointer :: json_data
    integer :: fi,ios,i,j
    logical :: file_exists
    real(dp) :: matr(7),ttt
    real(dp), dimension(:), allocatable :: xgas
    character(len=80) :: strval !name of the file with morphology
    character(len=99) :: nspec
    character(len=99) :: kspec
    character(len=99) :: gasspec
    character(len=255), dimension(:), allocatable :: gasname
    call get_names(gasname)
    allocate(xgas(size(gasname)))
    inputs=TRIM(ADJUSTL(fileplacein_par))//TRIM(ADJUSTL(inputs))
    inquire(file=inputs,exist=file_exists) !first try current folder
    if (.not. file_exists) then
        inputs='../'//inputs !then try one folder up
        inquire(file=inputs,exist=file_exists)
        if (.not. file_exists) stop 'input file not found'
    endif
    spectra=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(spectra))
    json_data => fson_parse(inputs)
    call fson_get(json_data, "upperBoundary.temperature", temp1)
    call fson_get(json_data, "lowerBoundary.temperature", temp2)
    call fson_get(json_data, "upperBoundary.emittance", emi1)
    call fson_get(json_data, "lowerBoundary.emittance", emi2)
    call fson_get(json_data, "gasDensity", rhog)
    call fson_get(json_data, "solidDensity", rhos)
    call fson_get(json_data, "sourceOfProperty.porosity", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        rhof=matr(1)
        por=1-rhof/rhos
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        rhof=matr(1)
        por=1-rhof/rhos
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        rhof=matr(1)
        por=1-rhof/rhos
        close(fi)
    elseif (strval=="DirectInput") then
        call fson_get(json_data, "porosity", por)
    else
        write(*,*) 'unknown source for porosity'
        stop
    endif
    rhof=(1-por)*rhos
    call fson_get(json_data, "sourceOfProperty.cellSize", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        dcell=matr(2)
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        dcell=matr(2)
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        dcell=matr(2)
        close(fi)
    elseif (strval=="DirectInput") then
        call fson_get(json_data, "cellSize", dcell)
    else
        write(*,*) 'unknown source for cell size'
        stop
    endif
    call fson_get(json_data, "sourceOfProperty.gasComposition", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xgas = 0
        ! this works only for cyclopentane and CO2
        xgas(1) = matr(4)
        xgas(2) = matr(3)
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xgas = 0
        ! this works only for cyclopentane and CO2
        xgas(1) = matr(4)
        xgas(2) = matr(3)
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xgas = 0
        ! this works only for cyclopentane and CO2
        xgas(1) = matr(4)
        xgas(2) = matr(3)
        close(fi)
    elseif (strval=="DirectInput") then
        do i=1,size(xgas)
            call fson_get(&
                json_data,&
                "gasComposition."//TRIM(ADJUSTL(gasname(i))),&
                xgas(i))
        enddo
    else
        write(*,*) 'unknown source for gas composition'
        stop
    endif
    xgas = xgas / sum(xgas)
    call fson_get(json_data, "sourceOfProperty.strutContent", strval)
    if (strval=="StrutContent") then
        call strutContent(fs,rhof)
    elseif (strval=="DirectInput") then
        call fson_get(json_data, "strutContent", fs)
    else
        write(*,*) 'unknown source for strut content'
        stop
    endif
    call fson_get(json_data, "sourceOfProperty.wallThickness", strval)
    if (strval=="DirectInput") then
        call fson_get(json_data, "wallThickness", dwall)
    else
        write(*,*) 'unknown source for wall thickness'
        stop
    endif
    call fson_get(json_data, "morphologyInput", strval)
    if (strval=="wallThickness") then
        morph_input=1
    elseif (strval=="strutContent") then
        morph_input=2
    elseif (strval=="strutSize") then
        morph_input=3
    elseif (strval=="strutContent2") then
        morph_input=4
    else
        write(*,*) 'unknown choice of morphologyInput'
        stop
    endif
    call fson_get(json_data, "strutSize", dstrut)
    call fson_get(json_data, "foamThickness", dfoam)
    call fson_get(json_data, "spectra.polymer_n", nspec)
    call fson_get(json_data, "spectra.polymer_k", kspec)
    call fson_get(json_data, "spectra.gas_k", gasspec)
    call fson_get(json_data, "spatialDiscretization", nz)
    call fson_get(json_data, "useWallThicknessDistribution", wdist)
    if (wdist) then
        call fson_get(json_data, "wallThicknessStandardDeviation", wsdev)
    endif
    call fson_get(json_data, "numberOfGrayBoxes", nbox)
    call fson_get(json_data, "numericalEffectiveConductivity", numcond)
    if (numcond) then
        call fson_get(json_data, "structureName", structureName)
    endif
    call fson_get(json_data, "testMode", testMode)
    if (temp1<temp2) then
        tmean=temp1
        temp1=temp2
        temp2=tmean
    endif
    tmean=(temp1+temp2)/2
    call gasConductivity(cond1,tmean,xgas,gasname)
    call polymerConductivity(cond2,tmean)
    cond2 = 0.1455_dp+1e-6_dp*(tmean - 273.2_dp)**2
    n1=1
    k1=0
    write(*,*) 'System information:'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'higher temperature:',temp1,'K'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'lower temperature: ',temp2,'K'
    write(*,*) 'Phase properties:'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'gas conductivity:  ',cond1*1e3_dp,'mW/m/K'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'solid conductivity:',cond2*1e3_dp,'mW/m/K'
    write(mfi,*) 'System information:'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'higher temperature:',temp1,'K'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'lower temperature: ',temp2,'K'
    write(mfi,*) 'Phase properties:'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') &
        'gas conductivity:  ',cond1*1e3_dp,'mW/m/K'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') &
        'solid conductivity:',cond2*1e3_dp,'mW/m/K'

    j=0
    open(newunit(fi),file=TRIM(ADJUSTL(fileplacein_ref))//TRIM(ADJUSTL(nspec)))
        do  !find number of points
            read(fi,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        rewind(fi)
        allocate(lambdan(j),nwl(j))
        do i=1,j
            read(fi,*) lambdan(i),nwl(i)
        enddo
        lambdan=lambdan*1e-6_dp
    close(fi)

    j=0
    open(newunit(fi),file=TRIM(ADJUSTL(fileplacein_ref))//TRIM(ADJUSTL(kspec)))
        do  !find number of points
            read(fi,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        rewind(fi)
        allocate(lambdak(j),kwl(j))
        do i=1,j
            read(fi,*) lambdak(i),kwl(i)
        enddo
        lambdak=lambdak*1e-6_dp
    close(fi)

    j=0
    open(newunit(fi),&
        file=TRIM(ADJUSTL(fileplacein_ref))//TRIM(ADJUSTL(gasspec)))
        do  !find number of points
            read(fi,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        rewind(fi)
        allocate(lambdagas(j),acgas(j))
        do i=1,j
            read(fi,*) lambdagas(i),acgas(i)
        enddo
        lambdagas=lambdagas*1e-6_dp
    close(fi)
end subroutine loadParameters
!***********************************END****************************************
end module tests
