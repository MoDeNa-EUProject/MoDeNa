!runs various parametric studies
!pavel.ferkl@vscht.cz
module tests
    use model
    implicit none
    private
    public onegrowth
contains
!********************************BEGINNING*************************************
!simulates one growth of a bubble
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
    concloc=fileplaceout
    inputs=TRIM(ADJUSTL(fileplacein))//TRIM(ADJUSTL(inputs))
    outputs_1d=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_1d))
    outputs_GR=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR))
    outputs_GR_c=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR_c))
    outputs_GR_p=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR_p))
    spar=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(spar))
    call read_inputs(inputs)
!    call save_surrogate_parameters(spar)
    call bblpreproc
    call bblinteg(outputs_1d,outputs_GR,outputs_GR_c,outputs_GR_p,concloc)
end subroutine onegrowth
!***********************************END****************************************
end module tests
