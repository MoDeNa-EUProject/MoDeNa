!> @file
!! subroutines for calculation of equivalent conductivity of foam, loading of
!! parameters and several parametric studies
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module tests
    use constants
    use ioutils
    use foamprop, only: effrad
    use foamgeom, only: foam_morpholgy
    use conduction, only: effcond
    use condrad, only: equcond,cond,alpha,sigma
    implicit none
    private
    public loadParameters,eqcond,eqcond_por,eqcond_dcell,eqcond_strut
    character(len=99) :: fileplacein_par='./'   !modena
    character(len=99) :: fileplacein_ref='../spectra/'  !modena
    character(len=99) :: fileplaceout='./'  !modena
    character(len=99) :: inputs='foamConductivity.json',spectra='spectra.out'
    character(len=99) :: nspec='spec_n.in'
    character(len=99) :: kspec='spec_k.in'
    character(len=99) :: gasspec='gasspec.in'
contains
!********************************BEGINNING*************************************
!> calculate equivalent conductivity for one specific foam
subroutine eqcond(regions)
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
            write(fi,*) effc
            close(fi)
            stop
        endif
        call foam_morpholgy
        if (testMode) then
            write(*,*) 'TESTING: radiative properties not calculated.'
            write(*,*) 'Ask Pavel if you want more reasonable results.'
            krad=2e-3_dp
            call effcond
            open(newunit(fi),file='foamConductivity.out')
            write(fi,*) effc
            close(fi)
            stop
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
!> calculate dependance of equivalent conductivity on porosity
subroutine eqcond_por
    integer :: fi,npoints,i
    real(dp) :: pormin,pormax,dpor
    pormin=0.90_dp
    pormax=0.995_dp
    npoints=20
    dpor=(pormax-pormin)/(npoints-1)
    open(newunit(fi),file='eqcond_por.csv')
    write(fi,*) 'porosity,foam_density,eq_cond,Ross_eq_cond,&
        kgas,ksol,krad'
    do i=1,npoints
        por=pormin+(i-1)*dpor
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
!> calculate dependance of equivalent conductivity on cell size
subroutine eqcond_dcell
    integer :: fi,npoints,i
    real(dp) :: dcellmin,dcellmax,ddcell
    dcellmin=1e-6_dp
    dcellmax=1000e-6_dp
    npoints=7
    ddcell=log10(dcellmax/dcellmin)/(npoints-1)
    open(newunit(fi),file='eqcond_dcell.out')
    write(fi,'(1000A23)') '#porosity','eq. conductivity'
    do i=1,npoints
        dcell=dcellmin*10**((i-1)*ddcell)
        call eqcond(1)
        write(fi,'(1000es23.15)') dcell,eqc
    enddo
    close(fi)
end subroutine eqcond_dcell
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate dependance of equivalent conductivity on strut content
subroutine eqcond_strut
    integer :: fi,npoints,i
    real(dp) :: strutmin,strutmax,dstrut
    strutmin=0.0_dp
    strutmax=0.9_dp
    npoints=19
    dstrut=(strutmax-strutmin)/(npoints-1)
    open(newunit(fi),file='eqcond_strut.out')
    write(fi,'(1000A23)') '#strut content','eq. conductivity','Ross. eq. cond.'
    do i=1,npoints
        fs=strutmin+(i-1)*dstrut
        call eqcond(1)
        write(fi,'(1000es23.15)') fs,eqc,eqc_ross
    enddo
    close(fi)
end subroutine eqcond_strut
!***********************************END****************************************


!********************************BEGINNING*************************************
!> loads parameters, usually from text file
subroutine loadParameters
    use physicalProperties
    use fson
    use fson_value_m, only: fson_value_get
    type(fson_value), pointer :: json_data
    integer :: fi,ios,i,j
    logical :: file_exists
    real(dp) :: xCO2,xAir,xCyP
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
    call fson_get(json_data, "gasComposition.CO2", xCO2)
    call fson_get(json_data, "gasComposition.Air", xAir)
    call fson_get(json_data, "gasComposition.Cyclopentane", xCyP)
    call fson_get(json_data, "gasDensity", rhog)
    call fson_get(json_data, "solidDensity", rhos)
    call fson_get(json_data, "porosity", por)
    call fson_get(json_data, "cellSize", dcell)
    call fson_get(json_data, "morphologyInput", morph_input)
    call fson_get(json_data, "wallThickness", dwall)
    call fson_get(json_data, "strutContent", fs)
    call fson_get(json_data, "strutSize", dstrut)
    call fson_get(json_data, "foamThickness", dfoam)
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
    call gasConductivity(cond1,tmean,xCO2,xAir,xCyP)
    call polymerConductivity(cond2,tmean)
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
