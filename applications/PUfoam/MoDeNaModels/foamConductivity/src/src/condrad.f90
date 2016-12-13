!> @file
!! subroutines for calculation of equivalent conductivity of the foam
!! coupled conduction-radiation simulation
!! non-gray 1D P1-approximation
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module condrad
    use constants
    use ioutils
    implicit none
    private
    public equcond,cond,alpha,sigma
    !lapack variables
    character :: trans='n'
    character(len=30) :: fmt='(2x,A,1x,es9.3,1x,A)'
    integer, parameter :: kl=1,ku=1,ldab=2*kl+ku+1,nrhs=1
    integer :: info
    integer, dimension(:), allocatable :: tipiv,gipiv
    real(dp), dimension(:), allocatable :: trhs,grhs
    real(dp), dimension(:,:), allocatable :: tmatrix,gmatrix
    !end of lapack variables
    real(dp) :: dz
    real(dp), dimension(:), allocatable :: tvec,qcon,qrad,qtot,gqrad,cond
    real(dp), dimension(:,:), allocatable :: alpha,sigma
contains
!********************************BEGINNING*************************************
!> determines equivalent conductivity of the foam, main algorithm
subroutine equcond
    integer :: i,maxiter=20,fi
    real(dp) :: tol=1e-5_dp,res
    write(*,*) 'Conduction-Radiation:'
    write(mfi,*) 'Conduction-Radiation:'
    dz=dfoam/nz
    allocate(tmatrix(ldab,nz),gmatrix(ldab,nbox*nz),tipiv(nz),gipiv(nbox*nz),&
        trhs(nz),grhs(nbox*nz),tvec(nz),qcon(nz),qrad(nz),qtot(nz),gqrad(nz))
    forall (i=1:nz)
        tvec(i)=temp1+(i-1)*(temp2-temp1)/(nz-1) !initial temperature profile
    end forall
    call make_tmatrix
    call dgbtrf(nz,nz,kl,ku,tmatrix,ldab,tipiv,info)
    call make_gmatrix
    call dgbtrf(nbox*nz,nbox*nz,kl,ku,gmatrix,ldab,gipiv,info)
    do i=1,maxiter
        call make_trhs
        call dgbtrs(trans,nz,kl,ku,nrhs,tmatrix,ldab,tipiv,trhs,nz,info)
        res=abs(maxval(trhs-tvec))
        write(*,'(5x,A,1x,I2,A,1x,es9.3)') 'iteration:',i, ', residual:',res
        write(mfi,'(5x,A,1x,I2,A,1x,es9.3)') 'iteration:',i, ', residual:',res
        tvec=trhs
        if (res<tol) exit
    enddo
    call heatflux
    eqc=sum(qtot)/nz*dfoam/(temp1-temp2)
    rcontr=sum(qrad)/sum(qtot)
    gcontr=(1-rcontr)*kgas/(kgas+ksol)
    scontr=(1-rcontr)*ksol/(kgas+ksol)
    kgas=gcontr*eqc
    ksol=scontr*eqc
    krad=rcontr*eqc
    write(*,fmt) 'equivalent conductivity:',eqc*1e3_dp,'mW/m/K'
    write(*,fmt) 'contribution of gas:',gcontr*1e2_dp,'%'
    write(*,fmt) 'contribution of solid:',scontr*1e2_dp,'%'
    write(*,fmt) 'contribution of radiation:',rcontr*1e2_dp,'%'
    write(mfi,fmt) 'equivalent conductivity:',eqc*1e3_dp,'mW/m/K'
    write(mfi,fmt) 'contribution of gas:',gcontr*1e2_dp,'%'
    write(mfi,fmt) 'contribution of solid:',scontr*1e2_dp,'%'
    write(mfi,fmt) 'contribution of radiation:',rcontr*1e2_dp,'%'
    open(newunit(fi),file='foamConductivity.out')
    write(fi,*) eqc
    close(fi)
    deallocate(tmatrix,gmatrix,tipiv,gipiv,trhs,grhs,tvec,qcon,qrad,qtot,gqrad)
end subroutine equcond
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates matrix for calculation of incidence radiation
subroutine make_gmatrix
    integer :: i,j,k,l
    gmatrix=0
    gipiv=0
    l=0
    do k=1,nbox
        l=l+1
        do i=2,nz-1
            l=l+1
            j=l
            gmatrix(kl+ku+1+l-j,j)=2/((alpha(i,k)+sigma(i,k))*dz) + &
                3*alpha(i,k)*dz
            j=l-1
            gmatrix(kl+ku+1+l-j,j)=-1/((alpha(i,k)+sigma(i,k))*dz)
            j=l+1
            gmatrix(kl+ku+1+l-j,j)=-1/((alpha(i,k)+sigma(i,k))*dz)
        enddo
        l=l+1
        i=(k-1)*nz+1
        j=i
        gmatrix(kl+ku+1+i-j,j)=1/((alpha(1,k)+sigma(1,k))*dz) + &
            3*alpha(1,k)*dz + &
            9*emi1/(2-emi1)*alpha(1,k)/4/(alpha(1,k)+sigma(1,k))
        j=i+1
        gmatrix(kl+ku+1+i-j,j)=-1/((alpha(1,k)+sigma(1,k))*dz) - &
            3*emi1/(2-emi1)*alpha(1,k)/4/(alpha(1,k)+sigma(1,k))
        i=nz*k
        j=i-1
        gmatrix(kl+ku+1+i-j,j)=-1/((alpha(nz,k)+sigma(nz,k))*dz) - &
            3*emi2/(2-emi2)*alpha(nz,k)/4/(alpha(nz,k)+sigma(nz,k))
        j=i
        gmatrix(kl+ku+1+i-j,j)=1/((alpha(nz,k)+sigma(nz,k))*dz) + &
            3*alpha(nz,k)*dz + &
            9*emi2/(2-emi2)*alpha(nz,k)/4/(alpha(nz,k)+sigma(nz,k))
    enddo
end subroutine make_gmatrix
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates right hand side for calculation of incidence radiation
subroutine make_grhs
    integer :: i,j,k
    k=0
    do i=1,nbox
        do j=1,nz
            k=k+1
            grhs(k)=12*effn**2*sigmab*alpha(j,i)*tvec(j)**4*dz*fbepbox(i)
        enddo
        j=(i-1)*nz+1
        grhs(j)=grhs(j) + 6*emi1/(2-emi1)*alpha(1,i)*effn**2*sigmab*temp1**4/&
            (alpha(1,i)+sigma(1,i))*fbepbox(i)
        j=nz*i
        grhs(j)=grhs(j) + 6*emi2/(2-emi2)*alpha(nz,i)*effn**2*sigmab*temp2**4/&
            (alpha(nz,i)+sigma(nz,i))*fbepbox(i)
    enddo
end subroutine make_grhs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates matrix for calculation of temperature
subroutine make_tmatrix
    integer :: i,j
    tmatrix=0
    tipiv=0
    i=1
    j=i
    tmatrix(kl+ku+1+i-j,j)=(cond(i) + &
        cond(i)*cond(i+1)/(cond(i)+cond(i+1)))*2/dz
    j=i+1
    tmatrix(kl+ku+1+i-j,j)=-cond(i)*cond(i+1)/(cond(i)+cond(i+1))*2/dz
    do i=2,nz-1
        j=i-1
        tmatrix(kl+ku+1+i-j,j)=-cond(i-1)*cond(i)/(cond(i-1)+cond(i))*2/dz
        j=i
        tmatrix(kl+ku+1+i-j,j)=(cond(i-1)*cond(i)/(cond(i-1)+cond(i)) + &
            cond(i)*cond(i+1)/(cond(i)+cond(i+1)))*2/dz
        j=i+1
        tmatrix(kl+ku+1+i-j,j)=-cond(i)*cond(i+1)/(cond(i)+cond(i+1))*2/dz
    enddo
    i=nz
    j=i-1
    tmatrix(kl+ku+1+i-j,j)=-cond(i-1)*cond(i)/(cond(i-1)+cond(i))*2/dz
    j=i
    tmatrix(kl+ku+1+i-j,j)=(cond(i-1)*cond(i)/(cond(i-1)+cond(i)) + &
        cond(i))*2/dz
end subroutine make_tmatrix
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates right hand side for calculation of temperature
subroutine make_trhs
    call make_gqrad
    trhs=gqrad
    trhs(1)=trhs(1)+2*cond(1)*temp1/dz
    trhs(nz)=trhs(nz)+2*cond(nz)*temp2/dz
end subroutine make_trhs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates gradient of radiative heat flux (radiative source)
subroutine make_gqrad
    integer :: i,j,k
    call make_grhs
    call dgbtrs(trans,nbox*nz,kl,ku,nrhs,gmatrix,ldab,gipiv,grhs,nbox*nz,info)
    k=0
    gqrad=0
    do i=1,nbox
        do j=1,nz
            k=k+1
            gqrad(j)=gqrad(j) + alpha(j,i)*grhs(k)*dz - &
                4*effn**2*sigmab*alpha(j,i)*tvec(j)**4*dz*fbepbox(i)
        enddo
    enddo
end subroutine make_gqrad
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate heat flux
subroutine heatflux
    integer :: i,j,k
    qcon(1)=-2*cond(1)*cond(2)/(cond(1)+cond(2))*(tvec(2)-tvec(1))/dz
    do i=2,nz
        qcon(i)=-2*cond(i-1)*cond(i)/(cond(i-1)+cond(i))*(tvec(i)-tvec(i-1))/dz
    enddo
    qrad=0
    k=0
    do i=1,nbox
        j=(i-1)*nz+1
        qrad(1)=qrad(1)-(grhs(j+1)-grhs(j))/(3*(alpha(1,i)+sigma(1,i))*dz)
        k=k+1
        do j=2,nz
            k=k+1
            qrad(j)=qrad(j)-(grhs(k)-grhs(k-1))/(3*(alpha(j,i)+sigma(j,i))*dz)
        enddo
    enddo
    qtot=qcon+qrad
end subroutine heatflux
!***********************************END****************************************
end module condrad
