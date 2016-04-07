!> @file
!! subroutines for growth of a single bubble and various parametric studies
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module tests
    use model
    implicit none
    private
    public onegrowth,eta_rm,bub_vf
contains
!********************************BEGINNING*************************************
!> simulates one growth of a bubble
subroutine onegrowth
    !HORNET windows
!    character(len=99) :: fileplacein='C:\Pavel\Dropbox\src\bubblegrowth_src\'
    !HORNET linux
!    character(len=99) :: fileplacein='/home/me/Dropbox/src/bubblegrowth_src/'
    !laptop windows
   ! character(len=99) :: fileplacein=&
   !      'C:\Users\pavel\Dropbox\src\bubblegrowth_src\'
!    character(len=99) :: fileplacein='./' !current folder
    character(len=99) :: fileplacein='../',& !modena
        fileplaceout='../results/',& !modena
        inputs='inputs.in',outputs_1d='outputs_1d.out',&
        outputs_GR='outputs_GR.out',outputs_GR_c='outputs_GR_c.out',&
        outputs_GR_p='outputs_GR_p.out',spar='GR_par.dat',concloc
    real(dp) :: tsave
    concloc=fileplaceout
    if (firstrun) then
        inputs=TRIM(ADJUSTL(fileplacein))//TRIM(ADJUSTL(inputs))
        outputs_1d=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_1d))
        outputs_GR=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR))
        outputs_GR_c=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR_c))
        outputs_GR_p=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR_p))
        spar=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(spar))
    else
        tsave=tend
    endif
    call read_inputs(inputs)
    if (.not. firstrun) then
        tend=tsave
    endif
!    call save_surrogate_parameters(spar)
    call bblpreproc
    call bblinteg(outputs_1d,outputs_GR,outputs_GR_c,outputs_GR_p,concloc)
    write(*,*) etat(1,1), etat(1,2), size(etat(:,1))
    write(*,*) eta_rm(80._dp)
    write(*,*) bub_vf(80._dp)
    deallocate(D,cbl,xgas,KH,fic,Mbl,dHv,mb,mb2,mb3,avconc,pressure,&
        diff_model,sol_model,cpblg,cpbll,RWORK,IWORK,dz,Y)
end subroutine onegrowth
!***********************************END****************************************


!********************************BEGINNING*************************************
!> viscosity of reaction mixture as function of time
real(dp) function eta_rm(t)
    use interpolation
    real(dp) :: t
    integer :: n
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=t
    call pwl_interp_1d ( size(etat(:,1)), etat(:,1), etat(:,2), ni, xi, yi )
    eta_rm=yi(1)
endfunction eta_rm
!***********************************END****************************************


!********************************BEGINNING*************************************
!> volume fraction of bubbles as function of time
real(dp) function bub_vf(t)
    use interpolation
    real(dp) :: t
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=t
    call pwl_interp_1d ( size(port(:,1)), port(:,1), port(:,2), ni, xi, yi )
    bub_vf=yi(1)
endfunction bub_vf
!***********************************END****************************************
end module tests
