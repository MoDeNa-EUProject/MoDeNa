!> @file
!! subroutines for calculation of effective conductivity of the foam
!! (conduction-only)
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module conduction
    use constants
    use ioutils
    implicit none
    private
    public effcond
    integer :: &
        magmit=1000,&       !<maximum number of multigrid iterations
        tnode,&             !<total number of nodes
        tp,&                !<number of non-zero elements
        dimz,dimy,dimx     !<dimensions of the morphology
    real(dp) :: &
        mtol=1e-30_dp       !<multigrid toleration
    integer,  dimension(:), allocatable :: &
        tja,&       !<column index
        tia         !<column of first non-zero element in row
    real(dp), dimension(:), allocatable :: &
        ta,&        !<coefficients of temperature sparse matrix
        trhs,&      !<righthandside for temperature matrix
        tvec        !<temperature vector
    integer, dimension(:,:,:), allocatable :: &
        sfiel,&     !<structure field (0=solid, 1=pore)
        ind         !<matrix index
    real(dp), dimension(:,:,:), allocatable :: &
        dx,dy,dz,&  !<voxel sizes in x,y,z
        cfiel,&     !<conductivity field
        chfz,&      !<conductive heat flux in z
        tfiel       !<temperature field
contains
!********************************BEGINNING*************************************
!> determine effective conductivity of the foam
subroutine effcond
    real(dp) :: xw,xs,f
    character(len=30) :: fmt='(2x,A,1x,es9.3,1x,A)'
    write(*,*) 'Conduction - algebraic:'
    write(mfi,*) 'Conduction - algebraic:'
    xw=2*(1+cond1/(2*cond2))/3
    xs=(1+4*cond1/(cond1+cond2))/3
    f=(1-fs)*xw+fs*xs
    kgas=cond1*por/(por+(1-por)*f)
    ksol=cond2*f*(1-por)/(por+(1-por)*f)
    effc=kgas+ksol
    eqc_ross=krad+effc
    gcontr=kgas/eqc_ross
    scontr=ksol/eqc_ross
    rcontr=krad/eqc_ross
    write(*,fmt) 'effective conductivity:',effc*1e3_dp,'mW/m/K'
    write(*,fmt) 'equivalent conductivity:',eqc_ross*1e3_dp,'mW/m/K'
    write(*,fmt) 'contribution of gas:',gcontr*1e2_dp,'%'
    write(*,fmt) 'contribution of solid:',scontr*1e2_dp,'%'
    write(*,fmt) 'contribution of radiation:',rcontr*1e2_dp,'%'
    write(mfi,fmt) 'effective conductivity:',effc*1e3_dp,'mW/m/K'
    write(mfi,fmt) 'equivalent conductivity:',eqc_ross*1e3_dp,'mW/m/K'
    write(mfi,fmt) 'contribution of gas:',gcontr*1e2_dp,'%'
    write(mfi,fmt) 'contribution of solid:',scontr*1e2_dp,'%'
    write(mfi,fmt) 'contribution of radiation:',rcontr*1e2_dp,'%'
    if (numcond) then
        write(*,*) 'Conduction - numerical:'
        write(mfi,*) 'Conduction - numerical:'
        structureName='../'//structureName
        call loadStructure_vtk(structureName)
        call initialization_oc_ss
        call calculation_oc_ss
        call heatFlux_oc_ss
        write(*,fmt) 'effective conductivity:',effc_num*1e3_dp,'mW/m/K'
        write(mfi,fmt) 'effective conductivity:',effc_num*1e3_dp,'mW/m/K'
    endif
end subroutine effcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!> loads structure from text file (vtk)
!! Beware of dimension reading.
!! May give wrong results for very large or very small numbers
subroutine loadStructure_vtk(filename)
    implicit none
    character(len=80), intent(in) :: filename
    character(len=80) :: dummy_char,datatype
    real(dp) :: a,b,c
    integer, dimension(:), allocatable :: bmat
    integer :: i,j,k,l,fi
    open(newunit(fi),file=filename)
    read(fi,*)
    read(fi,*)
    read(fi,*)
    read(fi,*)
    read(fi,*) dummy_char,dimz,dimy,dimx
    tnode=dimz*dimy*dimx
    allocate(sfiel(dimz+1,0:dimy+1,0:dimx+1),dx(dimz,0:dimy+1,0:dimx+1),&
        dy(dimz,0:dimy+1,0:dimx+1),dz(dimz,0:dimy+1,0:dimx+1))
    dz=dfoam/dimz
    dy=dz
    dx=dz
    sfiel=1
    read(fi,*)
    read(fi,*)
    read(fi,*)
    read(fi,*) datatype
    if (datatype=='SCALARS') then
        read(fi,*)
        do i=1,dimz
            do j=1,dimy
                read(fi,*) (sfiel(i,j,k),k=1,dimx)
            enddo
            read(fi,*)
        enddo
    elseif (datatype=='COLOR_SCALARS') then
        allocate(bmat(tnode))
        ! 3 values on first data line
        read(fi,*) a,b,c
        bmat(1)=int(a+0.5)
        bmat(2)=int(b+0.5)
        bmat(3)=int(c+0.5)
        do i=1,tnode/2-2
            ! 2 values on next lines\
            read(fi,*) a,b
            bmat(2*i)=int(a+0.5)
            bmat(2*i+1)=int(b+0.5)
        enddo
        ! 1 value on last line
        read(fi,*) a
        bmat(tnode)=int(a+0.5)
        l=1
        do i=1,dimz
            do j=1,dimy
                do k=1,dimx
                    sfiel(i,j,k)=bmat(l)
                    l=l+1
                enddo
            enddo
        enddo
        deallocate(bmat)
    else
        stop 'unknown VTK format'
    endif
    close(fi)
    do i=1,dimz
        do j=1,dimy
            sfiel(i,j,0)=sfiel(i,j,dimx)
            sfiel(i,j,dimx+1)=sfiel(i,j,1)
            dx(i,j,0)=dx(i,j,dimx)
            dx(i,j,dimx+1)=dx(i,j,1)
            dy(i,j,0)=dy(i,j,dimx)
            dy(i,j,dimx+1)=dy(i,j,1)
            dz(i,j,0)=dz(i,j,dimx)
            dz(i,j,dimx+1)=dz(i,j,1)
        enddo
    enddo
    do i=1,dimz
        do j=1,dimx
            sfiel(i,0,j)=sfiel(i,dimy,j)
            sfiel(i,dimy+1,j)=sfiel(i,1,j)
            dx(i,0,j)=dx(i,dimy,j)
            dx(i,dimy+1,j)=dx(i,1,j)
            dy(i,0,j)=dy(i,dimy,j)
            dy(i,dimy+1,j)=dy(i,1,j)
            dz(i,0,j)=dz(i,dimy,j)
            dz(i,dimy+1,j)=dz(i,1,j)
        enddo
    enddo
end subroutine loadStructure_vtk
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves structure to text file (vtk)
subroutine saveStructure_vtk(filename)
    implicit none
    character(len=80), intent(in) :: filename
    integer :: i,j,k,fi
    open(newunit(fi),file=filename)
    write(fi,'(A26)') '# vtk DataFile Version 3.0'
    write(fi,'(A7)') 'vtkfile'
    write(fi,'(A5)') 'ASCII'
    write(fi,'(A25)') 'DATASET STRUCTURED_POINTS'
    write(fi,'(A10,2x,I8,2x,I4,2x,I4)') 'DIMENSIONS',dimz,dimy,dimx
    write(fi,'(A12)') 'ORIGIN 0 0 0'
    write(fi,'(A13)') 'SPACING 1 1 1'
    write(fi,'(A10,2x,I8)') 'POINT_DATA',tnode
    write(fi,'(A18)') 'SCALARS values int'
    write(fi,'(A20)') 'LOOKUP_TABLE default'
    do i=1,dimz
        do j=1,dimy
            write(fi,'(1000I2)') (sfiel(i,j,k),k=1,dimx)
        enddo
        write(fi,*)
    enddo
    close(fi)
end subroutine saveStructure_vtk
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves structure to text file (vtk)
subroutine save3DField_vtk(filename,field)
    implicit none
    character(len=80), intent(in) :: filename
    real(dp), dimension(:,:,:), intent(in) :: field
    integer :: i,j,k,fi
    open(newunit(fi),file=filename)
    write(fi,'(A26)') '# vtk DataFile Version 3.0'
    write(fi,'(A7)') 'vtkfile'
    write(fi,'(A5)') 'ASCII'
    write(fi,'(A25)') 'DATASET STRUCTURED_POINTS'
    write(fi,'(A10,2x,I8,2x,I4,2x,I4)') 'DIMENSIONS',dimz,dimy,dimx
    write(fi,'(A12)') 'ORIGIN 0 0 0'
    write(fi,'(A13)') 'SPACING 1 1 1'
    write(fi,'(A10,2x,I8)') 'POINT_DATA',tnode
    write(fi,'(A20)') 'SCALARS values float'
    write(fi,'(A20)') 'LOOKUP_TABLE default'
    do i=1,dimz
        do j=1,dimy
            write(fi,'(1000E12.5)') (field(i,j,k),k=1,dimx)
        enddo
        write(fi,*)
    enddo
    close(fi)
end subroutine save3DField_vtk
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Initializes vectors and fields for numerical calculation
subroutine initialization_oc_ss
    integer :: i,j,k,l
    !assign conductivity to voxels
    allocate(cfiel(dimz,0:dimy+1,0:dimx+1))
    do i=1,dimz
        do j=1,dimy
            do k=1,dimx
                if (sfiel(i,j,k)==0) then
                    cfiel(i,j,k)=cond2
                else
                    cfiel(i,j,k)=cond1
                endif
            enddo
        enddo
    enddo
    do i=1,dimz
        do j=1,dimy
            cfiel(i,j,0)=cfiel(i,j,dimx)
            cfiel(i,j,dimx+1)=cfiel(i,j,1)
        enddo
    enddo
    do i=1,dimz
        do j=1,dimx
            cfiel(i,0,j)=cfiel(i,dimy,j)
            cfiel(i,dimy+1,j)=cfiel(i,1,j)
        enddo
    enddo
    !initial temperature field not needed
    allocate(tfiel(dimz,dimy,dimx))
    tfiel=0
    !vector from field
    allocate(tvec(dimx*dimy*dimz))
    l=0
    do i=1,dimz
        do j=1,dimy
            do k=1,dimx
                l=l+1
                tvec(l)=tfiel(i,j,k)
            enddo
        enddo
    enddo
    !index
    allocate(ind(dimz,0:dimy+1,0:dimx+1))
    l=0
    do i=1,dimz
        do j=1,dimy
            do k=1,dimx
                l=l+1
                ind(i,j,k)=l
            enddo
        enddo
    enddo
    do i=1,dimz
        do j=1,dimy
            ind(i,j,0)=ind(i,j,dimx)
            ind(i,j,dimx+1)=ind(i,j,1)
        enddo
    enddo
    do i=1,dimz
        do j=1,dimx
            ind(i,0,j)=ind(i,dimy,j)
            ind(i,dimy+1,j)=ind(i,1,j)
        enddo
    enddo
end subroutine initialization_oc_ss
!***********************************END****************************************


!********************************BEGINNING*************************************
!> main loop for numerical simulation of effective conductivity
subroutine calculation_oc_ss
    integer :: iter !maximum number of multigrid iterations
    allocate(ta(7*tnode),tja(7*tnode),tia(tnode+1),trhs(tnode))
    call make_tmtrx_oc_ss   !fill temperature matrix
    call make_trhs_oc_ss    !make right hand side
    !set maximum number of iterations iter for dagmg (it changes in dagmg)
    iter=magmit
    call dagmg(tnode,ta(1:tp),tja(1:tp),tia,trhs,tvec,10,1,10,iter,mtol)
    if (iter<0) then
        write(*,*) 'temperature was not calculated'
        stop
    endif
    deallocate(ta,tja,tia,trhs,ind)
end subroutine calculation_oc_ss
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates matrix for calculation of temperature
subroutine make_tmtrx_oc_ss
    integer :: i,j,k,l,m
    l=1
    m=1
    do i=1,dimy
        do j=1,dimx
            ta(l)=1
            tja(l)=ind(1,i,j)
            tia(m)=l
            l=l+1
            m=m+1
        enddo
    enddo
    do i=2,dimz-1
        do j=1,dimy
            do k=1,dimx
                ta(l)=-cfiel(i-1,j,k)*cfiel(i,j,k)*(dz(i-1,j,k)+dz(i,j,k))/&
                    (cfiel(i-1,j,k)*dz(i,j,k)+cfiel(i,j,k)*dz(i-1,j,k))/&
                    (dz(i-1,j,k)/2+dz(i,j,k)/2)*dx(i,j,k)*dy(i,j,k)
                tja(l)=ind(i-1,j,k)
                tia(m)=l
                l=l+1
                m=m+1
                ta(l)=-cfiel(i,j-1,k)*cfiel(i,j,k)*(dy(i,j-1,k)+dy(i,j,k))/&
                    (cfiel(i,j-1,k)*dy(i,j,k)+cfiel(i,j,k)*dy(i,j-1,k))/&
                    (dy(i,j-1,k)/2+dy(i,j,k)/2)*dx(i,j,k)*dz(i,j,k)
                tja(l)=ind(i,j-1,k)
                l=l+1
                ta(l)=-cfiel(i,j,k-1)*cfiel(i,j,k)*(dx(i,j,k-1)+dx(i,j,k))/&
                    (cfiel(i,j,k-1)*dx(i,j,k)+cfiel(i,j,k)*dx(i,j,k-1))/&
                    (dx(i,j,k-1)/2+dx(i,j,k)/2)*dy(i,j,k)*dz(i,j,k)
                tja(l)=ind(i,j,k-1)
                l=l+1
                ta(l)=0 &
                    +cfiel(i-1,j,k)*cfiel(i,j,k)*(dz(i-1,j,k)+dz(i,j,k))/&
                        (cfiel(i-1,j,k)*dz(i,j,k)+cfiel(i,j,k)*dz(i-1,j,k))/&
                        (dz(i-1,j,k)/2+dz(i,j,k)/2)*dx(i,j,k)*dy(i,j,k) &
                    +cfiel(i,j-1,k)*cfiel(i,j,k)*(dy(i,j-1,k)+dy(i,j,k))/&
                        (cfiel(i,j-1,k)*dy(i,j,k)+cfiel(i,j,k)*dy(i,j-1,k))/&
                        (dy(i,j-1,k)/2+dy(i,j,k)/2)*dx(i,j,k)*dz(i,j,k) &
                    +cfiel(i,j,k-1)*cfiel(i,j,k)*(dx(i,j,k-1)+dx(i,j,k))/&
                        (cfiel(i,j,k-1)*dx(i,j,k)+cfiel(i,j,k)*dx(i,j,k-1))/&
                        (dx(i,j,k-1)/2+dx(i,j,k)/2)*dy(i,j,k)*dz(i,j,k) &
                    +cfiel(i,j,k)*cfiel(i,j,k+1)*(dx(i,j,k)+dx(i,j,k+1))/&
                        (cfiel(i,j,k)*dx(i,j,k+1)+cfiel(i,j,k+1)*dx(i,j,k))/&
                        (dx(i,j,k)/2+dx(i,j,k+1)/2)*dy(i,j,k)*dz(i,j,k) &
                    +cfiel(i,j,k)*cfiel(i,j+1,k)*(dy(i,j,k)+dy(i,j+1,k))/&
                        (cfiel(i,j,k)*dy(i,j+1,k)+cfiel(i,j+1,k)*dy(i,j,k))/&
                        (dy(i,j,k)/2+dy(i,j+1,k)/2)*dx(i,j,k)*dz(i,j,k) &
                    +cfiel(i,j,k)*cfiel(i+1,j,k)*(dz(i,j,k)+dz(i+1,j,k))/&
                        (cfiel(i,j,k)*dz(i+1,j,k)+cfiel(i+1,j,k)*dz(i,j,k))/&
                        (dz(i,j,k)/2+dz(i+1,j,k)/2)*dx(i,j,k)*dy(i,j,k)
                tja(l)=ind(i,j,k)
                l=l+1
                ta(l)=-cfiel(i,j,k)*cfiel(i,j,k+1)*(dx(i,j,k)+dx(i,j,k+1))/&
                    (cfiel(i,j,k)*dx(i,j,k+1)+cfiel(i,j,k+1)*dx(i,j,k))/&
                    (dx(i,j,k)/2+dx(i,j,k+1)/2)*dy(i,j,k)*dz(i,j,k)
                tja(l)=ind(i,j,k+1)
                l=l+1
                ta(l)=-cfiel(i,j,k)*cfiel(i,j+1,k)*(dy(i,j,k)+dy(i,j+1,k))/&
                    (cfiel(i,j,k)*dy(i,j+1,k)+cfiel(i,j+1,k)*dy(i,j,k))/&
                    (dy(i,j,k)/2+dy(i,j+1,k)/2)*dx(i,j,k)*dz(i,j,k)
                tja(l)=ind(i,j+1,k)
                l=l+1
                ta(l)=-cfiel(i,j,k)*cfiel(i+1,j,k)*(dz(i,j,k)+dz(i+1,j,k))/&
                    (cfiel(i,j,k)*dz(i+1,j,k)+cfiel(i+1,j,k)*dz(i,j,k))/&
                    (dz(i,j,k)/2+dz(i+1,j,k)/2)*dx(i,j,k)*dy(i,j,k)
                tja(l)=ind(i+1,j,k)
                l=l+1
            enddo
        enddo
    enddo
    do i=1,dimy
        do j=1,dimx
            ta(l)=1
            tja(l)=ind(dimz,i,j)
            tia(m)=l
            l=l+1
            m=m+1
        enddo
    enddo
    tia(m)=l
    tp=l-1
end subroutine make_tmtrx_oc_ss
!***********************************END****************************************


!********************************BEGINNING*************************************
!creates right hand side for calculation of temperature
subroutine make_trhs_oc_ss
    integer :: i,j,k,l
    l=1
    do i=1,dimy
        do j=1,dimx
            trhs(l)=temp2
            l=l+1
        enddo
    enddo
    do i=2,dimz-1
        do j=1,dimy
            do k=1,dimx
                trhs(l)=0
                l=l+1
            enddo
        enddo
    enddo
    do i=1,dimy
        do j=1,dimx
            trhs(l)=temp1
            l=l+1
        enddo
    enddo
end subroutine make_trhs_oc_ss
!***********************************END****************************************


!********************************BEGINNING*************************************
!calculates heat flux
subroutine heatFlux_oc_ss
    integer :: i,j,k,l
    real(dp) :: xxx,area
    character(len=80) :: filename
    allocate(chfz(dimz,dimy,dimx))
    l=1
    do i=1,dimz
        do j=1,dimy
            do k=1,dimx
                tfiel(i,j,k)=tvec(l)
                l=l+1
            enddo
        enddo
    enddo
    do i=2,dimz-1
        do j=1,dimy
            do k=1,dimx
                if (sfiel(i,j,k)==sfiel(i+1,j,k)) then
                    chfz(i,j,k)=-(tfiel(i,j,k)-tfiel(i+1,j,k))/&
                        (dz(i,j,k)/2+dz(i+1,j,k)/2)*cfiel(i,j,k)
                else
                    chfz(i,j,k)=(tfiel(i,j,k)-tfiel(i-1,j,k))/&
                        (dz(i,j,k)/2+dz(i-1,j,k)/2)*cfiel(i,j,k)
                endif
            enddo
        enddo
    enddo
    do j=1,dimy
        do k=1,dimx
            chfz(1,j,k)=-(tfiel(1,j,k)-tfiel(2,j,k))/&
                (dz(1,j,k)/2+dz(2,j,k)/2)*cfiel(1,j,k)
            chfz(dimz,j,k)=(tfiel(dimz,j,k)-tfiel(dimz-1,j,k))/&
                (dz(dimz,j,k)/2+dz(dimz-1,j,k)/2)*cfiel(dimz,j,k)
        enddo
    enddo
    xxx=0
    do i=1,dimz
        effc_num=0
        area=0
        do j=1,dimy
            do k=1,dimx
                effc_num=effc_num+chfz(i,j,k)*dx(i,j,k)*dy(i,j,k)
                area=area+dx(i,j,k)*dy(i,j,k)
            enddo
        enddo
        xxx=xxx+effc_num/area/(temp1-temp2)*dfoam
    enddo
    effc_num=xxx/dimz
    filename='../results/Temp.vtk'
    ! call save3DField_vtk(filename,tfiel)
    deallocate(chfz,tfiel,tvec,cfiel,sfiel,dx,dy,dz)
end subroutine heatFlux_oc_ss
!***********************************END****************************************
end module conduction
