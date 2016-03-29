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
    character(len=99) :: inputs='foamConductivity.in',spectra='spectra.out'
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
        call foam_morpholgy
        if (testing) then
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
    open(newunit(fi),file='eqcond_por.out')
    write(fi,'(1000A23)') '#porosity','eq. conductivity','Ross. eq. cond.'
    do i=1,npoints
        por=pormin+(i-1)*dpor
        call eqcond(1)
        write(fi,'(1000es23.15)') por,eqc,eqc_ross
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
    integer :: fi,ios,i,j
    logical :: file_exists
    real(dp) :: xCO2,xAir,xCyP
    inputs=TRIM(ADJUSTL(fileplacein_par))//TRIM(ADJUSTL(inputs))
    spectra=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(spectra))
    inquire(file=inputs,exist=file_exists)
    if (file_exists) then
        open(newunit(fi),file=inputs)
    else
        open(newunit(fi),file='../'//inputs)
    endif
        read(fi,*) T1           !higher temperature
        read(fi,*) T2           !lower temperature
        read(fi,*) xCO2,xAir,xCyP
        call gasConductivity(cond1,(t1+t2)/2,xCO2,xAir,xCyP)
        ! read(fi,*) cond1        !gas conductivity
!        read(fi,*) cond2        !solid conductivity
        call polymerConductivity(cond2,(t1+t2)/2)
        read(fi,*) emi1         !emittance 1
        read(fi,*) emi2         !emittance 2
        read(fi,*) rho1         !gas density
        read(fi,*) rho2         !solid density
        read(fi,*) por          !porosity
        read(fi,*) dcell        !cell size
        read(fi,*) morph_input  !morphology input 1=wall thickness,
        ! 2=strut content, 3=strut diameter (3 is recommended others can have
        ! multiple solutions)
        read(fi,*) dwall        !wall thickness
        read(fi,*) fs           !strut content
        read(fi,*) dstrut       !strut diameter
        read(fi,*) dfoam        !foam thickness
        read(fi,*) nz           !spatial discretization
        read(fi,*) nrays        !number of testing rays
        read(fi,*) wdist        !use wall thickness distribution
        read(fi,*) wsdev        !wall thickness standard deviation
        read(fi,*) nbox         !number of gray boxes
        read(fi,*) numcond      !calcualte effective conductivity numerically
        read(fi,*) structureName!name of the file with morphology
    close(fi)
    tmean=(t1+t2)/2
    n1=1
    k1=0
    write(*,*) 'System information:'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'higher temperature:',T1,'K'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'lower temperature: ',T2,'K'
    write(*,*) 'Phase properties:'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'gas conductivity:  ',cond1*1e3_dp,'mW/m/K'
    write(*,'(2x,A,1x,es9.3,1x,A)') 'solid conductivity:',cond2*1e3_dp,'mW/m/K'
    write(mfi,*) 'System information:'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'higher temperature:',T1,'K'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'lower temperature: ',T2,'K'
    write(mfi,*) 'Phase properties:'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'gas conductivity:  ',cond1*1e3_dp,'mW/m/K'
    write(mfi,'(2x,A,1x,es9.3,1x,A)') 'solid conductivity:',cond2*1e3_dp,'mW/m/K'

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
    open(newunit(fi),file=TRIM(ADJUSTL(fileplacein_ref))//TRIM(ADJUSTL(gasspec)))
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
