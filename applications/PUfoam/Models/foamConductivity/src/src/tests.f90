!runs various parametric studies
!author: pavel.ferkl@vscht.cz
module tests
    use constants
    use ioutils
    use foamprop
    use foamgeom
    use conduction
    use condrad
    implicit none
    private
    public loadParameters,eqcond,eqcond_por,eqcond_dcell,eqcond_strut
    character(len=99) :: fileplacein_par='./'   !modena
    character(len=99) :: fileplacein_ref='../spectra/'  !modena
    character(len=99) :: fileplaceout='./'  !modena
    character(len=99) :: inputs='inputs.in',spectra='spectra.out'
    character(len=99) :: nspec='spec_n.in'
    character(len=99) :: kspec='spec_k.in'
    character(len=99) :: gasspec='gasspec.in'
contains
!********************************BEGINNING*************************************
!calculate equivalent conductivity
subroutine eqcond
    call foam_morpholgy
    call effrad(spectra)
    call effcond
    call equcond
end subroutine eqcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!calculate dependance of equivalent conductivity on porosity
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
        call eqcond
        write(fi,'(1000es23.15)') por,eqc,eqc_ross
    enddo
    close(fi)
end subroutine eqcond_por
!***********************************END****************************************


!********************************BEGINNING*************************************
!calculate dependance of equivalent conductivity on cell size
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
        call eqcond
        write(fi,'(1000es23.15)') dcell,eqc
    enddo
    close(fi)
end subroutine eqcond_dcell
!***********************************END****************************************


!********************************BEGINNING*************************************
!calculate dependance of equivalent conductivity on strut content
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
        call eqcond
        write(fi,'(1000es23.15)') fs,eqc,eqc_ross
    enddo
    close(fi)
end subroutine eqcond_strut
!***********************************END****************************************


!********************************BEGINNING*************************************
!loads parameters, usually from text file
subroutine loadParameters
    use physicalProperties
    integer :: fi,ios,i,j
    inputs=TRIM(ADJUSTL(fileplacein_par))//TRIM(ADJUSTL(inputs))
    spectra=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(spectra))
    open(newunit(fi),file=inputs)
        read(fi,*) T1           !higher temperature
        read(fi,*) T2           !lower temperature
        read(fi,*) cond1        !pore conductivity
!        read(fi,*) cond2        !solid conductivity
        call polymerConductivity(cond2,(t1+t2)/2)
        read(fi,*) emi1         !emittance 1
        read(fi,*) emi2         !emittance 2
        read(fi,*) rho1         !pore density
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
