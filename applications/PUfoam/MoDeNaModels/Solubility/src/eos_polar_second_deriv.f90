
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine A_POLAR_drhoi_drhoj( n_comp, A_polar_rr )

  USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n_comp
  real, dimension(n_comp,n_comp), intent(out)  :: A_polar_rr

  !-----------------------------------------------------------------------------
  integer                                :: dipole
  integer                                :: quadrupole
  integer                                :: dipole_quad
  real, allocatable, dimension(:,:)      :: Add_rr, Aqq_rr, Adq_rr
  !-----------------------------------------------------------------------------

  A_polar_rr(:,:) = 0.0

  dipole      = 0
  quadrupole  = 0
  dipole_quad = 0
  if ( SUM( parame(1:ncomp,6) ) /= 0.0 ) dipole = 1
  if ( SUM( parame(1:ncomp,7) ) /= 0.0 ) quadrupole = 1
  if ( dipole == 1 .AND. quadrupole == 1 ) dipole_quad = 1

  !-----------------------------------------------------------------------------
  ! dipole-dipole term
  !-----------------------------------------------------------------------------
  if (dipole == 1) then

     allocate( Add_rr(ncomp, ncomp) )
     if (dd_term == 'GV') CALL A_rr_DD_GROSS_VRABEC( n_comp, Add_rr )
     if (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

     A_polar_rr = A_polar_rr + Add_rr
     deallocate( Add_rr )

  end if

  !-----------------------------------------------------------------------------
  ! quadrupole-quadrupole term
  !-----------------------------------------------------------------------------
  if (quadrupole == 1) then

     allocate( Aqq_rr(ncomp, ncomp) )
     if (qq_term == 'JG') CALL A_rr_QQ_GROSS( n_comp, Aqq_rr )
     if (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

     A_polar_rr = A_polar_rr + Aqq_rr
     deallocate( Aqq_rr )

  end if

  !-----------------------------------------------------------------------------
  ! dipole-quadrupole cross term
  !-----------------------------------------------------------------------------
  if (dipole_quad == 1) then

     allocate( Adq_rr(ncomp, ncomp) )
     if (dq_term == 'VG') CALL A_rr_DQ_VRABEC_GROSS( n_comp, Adq_rr )
     if (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

     A_polar_rr = A_polar_rr + Adq_rr
     deallocate( Adq_rr )

  end if

end subroutine A_POLAR_drhoi_drhoj


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine A_rr_DD_GROSS_VRABEC( n_comp, Add_rr )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t, dhs
  USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  integer, intent(in)                   :: n_comp
  real, dimension(n_comp,n_comp), INTENT(IN OUT)  :: Add_rr

  !-----------------------------------------------------------------------------
  integer                                :: i, j, k, l, m, o

  real                                   :: factor2, factor3, eij
  real                                   :: z3, z3_m
  real                                   :: fdd2, fdd3
  real, allocatable, dimension(:,:)      :: xijfa
  real, allocatable, dimension(:,:,:)    :: xijkfa

  real, dimension(nc)                    :: my2dd, z3_r
  real, dimension(nc)                    :: rhoi
  real, allocatable, dimension(:,:)      :: Idd2
  real, allocatable, dimension(:,:,:)    :: Idd3
  real, allocatable, dimension(:,:,:)    :: Idd2_r
  real, allocatable, dimension(:,:,:,:)  :: Idd3_r
  real, allocatable, dimension(:)        :: fdd2_r, fdd3_r, fdd_r
  real, allocatable, dimension(:,:)      :: Idd2_rkrl
  real, allocatable, dimension(:,:,:)    :: Idd3_rkrl
  real                                   :: fdd2_rkrl, fdd3_rkrl
  !-----------------------------------------------------------------------------

  z3 = eta
  do i = 1, ncomp
     if ( uij(i,i) == 0.0 ) write (*,*) 'A_rr_DD_GROSS_VRABEC: do not use dimensionless units'
     if ( uij(i,i) == 0.0 ) stop
     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
     z3_r(i) = PI/6.0 * mseg(i) * dhs(i)**3
  end do

  !-----------------------------------------------------------------------------
  ! some basic quantities
  !-----------------------------------------------------------------------------

  allocate( xijfa(ncomp,ncomp), xijkfa(ncomp,ncomp,ncomp) )
  allocate( Idd2(ncomp,ncomp), Idd2_r(ncomp,ncomp,ncomp) )
  allocate( Idd3(ncomp,ncomp,ncomp), Idd3_r(ncomp,ncomp,ncomp,ncomp) )
  allocate( fdd2_r(ncomp), fdd3_r(ncomp), fdd_r(ncomp) )
  allocate( Idd2_rkrl(ncomp,ncomp), Idd3_rkrl(ncomp,ncomp,ncomp) )

  rhoi( 1:ncomp ) = x (1:ncomp ) * rho

  factor2= -PI
  factor3= -4.0/3.0*PI*PI

  do i = 1, ncomp
     do j = 1, ncomp
           xijfa(i,j) = factor2* uij(i,i)*my2dd(i)*sig_ij(i,i)**3 /t  &
                                *uij(j,j)*my2dd(j)*sig_ij(j,j)**3 /t / sig_ij(i,j)**3 
           do l = 1, ncomp
              xijkfa(i,j,l)= factor3*uij(i,i)/t*my2dd(i)*sig_ij(i,i)**3   &
                                    *uij(j,j)/t*my2dd(j)*sig_ij(j,j)**3   &
                                    *uij(l,l)/t*my2dd(l)*sig_ij(l,l)**3   &
                                    /sig_ij(i,j)/sig_ij(i,l)/sig_ij(j,l)
           end do
     end do
  end do

  do i = 1, ncomp
     do j = 1, ncomp
        Idd2(i,j)  = 0.0
        if (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) then
           do m = 0, 4
              eij = (parame(i,3)*parame(j,3))**0.5
              Idd2(i,j)  =Idd2(i,j) + ( ddp2(i,j,m) + eij/t*ddp4(i,j,m))*z3**m
           end do
           do l = 1, ncomp
              Idd3(i,j,l)  = 0.0
              if (parame(l,6) /= 0.0) then
                 do m = 0, 4
                    Idd3(i,j,l) = Idd3(i,j,l) + ddp3(i,j,l,m)*z3**m
                 end do
              end if
           end do
        end if
     end do
  end do

  fdd2  = 0.0
  fdd3  = 0.0
  do i = 1, ncomp
     do j = 1, ncomp
        if (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) then
           fdd2 = fdd2 + rhoi(i)*rhoi(j)*xijfa(i,j) * Idd2(i,j)
           do l=1,ncomp
              if (parame(l,6) /= 0.0) then
                 fdd3 = fdd3 + rhoi(i)*rhoi(j)*rhoi(l)*xijkfa(i,j,l) * Idd3(i,j,l)
              end if
           end do
        end if
     end do
  end do


  !-----------------------------------------------------------------------------
  ! some first derivatives
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

  do i = 1, ncomp
     do j = 1, ncomp
        Idd2_r(k,i,j) = 0.0
        if (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) then
           do m = 1, 4
              z3_m = REAL(m) * z3**(m-1) * z3_r(k)
              eij = (parame(i,3)*parame(j,3))**0.5
              Idd2_r(k,i,j) = Idd2_r(k,i,j)+ ( ddp2(i,j,m) + eij/t * ddp4(i,j,m) ) * z3_m
           end do
           do o = 1, ncomp
              Idd3_r(k,i,j,o) = 0.0
              if (parame(o,6) /= 0.0) then
                 do m = 1, 4
                    Idd3_r(k,i,j,o)=Idd3_r(k,i,j,o) + ddp3(i,j,o,m)*REAL(m)*z3**(m-1)*z3_r(k)
                 end do
              end if
           end do
        end if
     end do
  end do

  fdd2_r(k) = 0.0
  fdd3_r(k) = 0.0
  do i = 1, ncomp
     fdd2_r(k) = fdd2_r(k) + 2.0*rhoi(i)*xijfa(i,k) * Idd2(i,k)
     do j = 1, ncomp
        if (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) then
           fdd2_r(k) = fdd2_r(k) + rhoi(i)*rhoi(j)*xijfa(i,j) * Idd2_r(k,i,j)
           fdd3_r(k) = fdd3_r(k) + 3.0 * rhoi(i)*rhoi(j)*xijkfa(i,j,k) * Idd3(i,j,k)
           do o = 1, ncomp
              if (parame(o,6) /= 0.0) then
                 fdd3_r(k) =fdd3_r(k) + rhoi(i)*rhoi(j)*rhoi(o)*xijkfa(i,j,o) * Idd3_r(k,i,j,o)
              end if
           end do
        end if
     end do
  end do

  if (fdd2 < -1.E-50 .AND. fdd3 /= 0.0 .AND. fdd2_r(k) /= 0.0 .AND. fdd3_r(k) /= 0.0) then

     fdd_r(k) = fdd2* (fdd2*fdd2_r(k) - 2.0*fdd3*fdd2_r(k) + fdd2*fdd3_r(k)) / (fdd2-fdd3)**2

  end if

  end do

  !-----------------------------------------------------------------------------
  ! second derivatives
  !-----------------------------------------------------------------------------

  do k = 1, ncomp
  do l = 1, ncomp

  do i = 1, ncomp
     do j = 1, ncomp
        Idd2_rkrl(i,j) = 0.0
        if (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) then
           do m = 2, 4
              z3_m = REAL(m-1) * REAL(m) * z3**(m-2) * z3_r(k) * z3_r(l)
              eij = (parame(i,3)*parame(j,3))**0.5
              Idd2_rkrl(i,j) = Idd2_rkrl(i,j)+ ( ddp2(i,j,m) + eij/t * ddp4(i,j,m) ) * z3_m
           end do
           do o = 1, ncomp
              Idd3_rkrl(i,j,o) = 0.0
              if (parame(o,6) /= 0.0) then
                 do m = 2, 4
                    Idd3_rkrl(i,j,o)=Idd3_rkrl(i,j,o) + ddp3(i,j,o,m)*REAL(m-1)*REAL(m)*z3**(m-2)*z3_r(k)*z3_r(l)
                 end do
              end if
           end do
        end if
     end do
  end do

  fdd2_rkrl = 2.0 * xijfa(k,l) * Idd2(k,l)
  fdd3_rkrl = 0.0
  do i = 1, ncomp
     fdd2_rkrl = fdd2_rkrl + 2.0 * rhoi(i) * ( xijfa(i,k) * Idd2_r(l,i,k) +  xijfa(i,l) * Idd2_r(k,i,l))
     fdd3_rkrl = fdd3_rkrl + 6.0 * rhoi(i) * xijkfa(i,k,l) * Idd3(i,k,l)
     do j = 1, ncomp
        if (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) then
           fdd2_rkrl = fdd2_rkrl + rhoi(i)*rhoi(j)*xijfa(i,j) * Idd2_rkrl(i,j)
           fdd3_rkrl = fdd3_rkrl + 3.0 * rhoi(i)*rhoi(j)* ( xijkfa(i,j,k) * Idd3_r(l,i,j,k)  &
                                                          + xijkfa(i,j,l) * Idd3_r(k,i,j,l) )
           do o = 1, ncomp
              if (parame(o,6) /= 0.0) then
                 fdd3_rkrl =fdd3_rkrl + rhoi(i)*rhoi(j)*rhoi(o)*xijkfa(i,j,o) * Idd3_rkrl(i,j,o)
              end if
           end do
        end if
     end do
  end do


  Add_rr(k,l) = ( 2.0*fdd2*fdd2_r(l)*( fdd2_r(k) + fdd3_r(k) ) + fdd2*fdd2*( fdd2_rkrl + fdd3_rkrl )  &
                         -2.0*( fdd2_r(k)*fdd2_r(l)*fdd3 + fdd2*fdd3_r(l)*fdd2_r(k) + fdd2*fdd3*fdd2_rkrl ) )  &
                / ( fdd2 - fdd3 )**2  &
                  + fdd_r(k) * 2.0 * ( fdd3_r(l) - fdd2_r(l) ) / ( fdd2 - fdd3 )

  end do
  end do

  deallocate( xijfa, xijkfa, Idd2, Idd2_r, Idd3, Idd3_r, fdd2_r, fdd3_r, fdd_r, Idd2_rkrl, Idd3_rkrl )

end subroutine A_rr_DD_GROSS_VRABEC




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine A_rr_QQ_GROSS( n_comp, Aqq_rr )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t, dhs
  USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  integer, intent(in)                   :: n_comp
  real, dimension(n_comp,n_comp), INTENT(IN OUT)  :: Aqq_rr

  !-----------------------------------------------------------------------------
  integer                                :: i, j, k, l, m, o

  real                                   :: factor2, factor3, eij
  real                                   :: z3, z3_m
  real                                   :: fqq2, fqq3
  real, allocatable, dimension(:,:)      :: xijfa
  real, allocatable, dimension(:,:,:)    :: xijkfa

  real, dimension(nc)                    :: qq2, z3_r
  real, dimension(nc)                    :: rhoi
  real, allocatable, dimension(:,:)      :: Iqq2
  real, allocatable, dimension(:,:,:)    :: Iqq3
  real, allocatable, dimension(:,:,:)    :: Iqq2_r
  real, allocatable, dimension(:,:,:,:)  :: Iqq3_r
  real, allocatable, dimension(:)        :: fqq2_r, fqq3_r, fqq_r
  real, allocatable, dimension(:,:)      :: Iqq2_rkrl
  real, allocatable, dimension(:,:,:)    :: Iqq3_rkrl
  real                                   :: fqq2_rkrl, fqq3_rkrl
  !-----------------------------------------------------------------------------

  z3 = eta
  do i = 1, ncomp
     if ( uij(i,i) == 0.0 ) write (*,*) 'A_rr_QQ_GROSS: do not use dimensionless units'
     if ( uij(i,i) == 0.0 ) stop
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
     z3_r(i) = PI/6.0 * mseg(i) * dhs(i)**3
  end do

  !-----------------------------------------------------------------------------
  ! some basic quantities
  !-----------------------------------------------------------------------------

  allocate( xijfa(ncomp,ncomp), xijkfa(ncomp,ncomp,ncomp) )
  allocate( Iqq2(ncomp,ncomp), Iqq2_r(ncomp,ncomp,ncomp) )
  allocate( Iqq3(ncomp,ncomp,ncomp), Iqq3_r(ncomp,ncomp,ncomp,ncomp) )
  allocate( fqq2_r(ncomp), fqq3_r(ncomp), fqq_r(ncomp) )
  allocate( Iqq2_rkrl(ncomp,ncomp), Iqq3_rkrl(ncomp,ncomp,ncomp) )

  rhoi( 1:ncomp ) = x (1:ncomp ) * rho

  factor2= -9.0/16.0*PI
  factor3=  9.0/16.0*PI**2

  do i = 1, ncomp
     do j = 1, ncomp
           xijfa(i,j) = factor2* uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t  &
                                *uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t / sig_ij(i,j)**7 
           do l = 1, ncomp
              xijkfa(i,j,l)= factor3*uij(i,i)/t*qq2(i)*sig_ij(i,i)**5   &
                                    *uij(j,j)/t*qq2(j)*sig_ij(j,j)**5   &
                                    *uij(l,l)/t*qq2(l)*sig_ij(l,l)**5   &
                                    / ( sig_ij(i,j)*sig_ij(i,l)*sig_ij(j,l) )**3
           end do
     end do
  end do

  do i = 1, ncomp
     do j = 1, ncomp
        Iqq2(i,j)  = 0.0
        if (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) then
           do m = 0, 4
              eij = (parame(i,3)*parame(j,3))**0.5
              Iqq2(i,j) = Iqq2(i,j) + ( qqp2(i,j,m) + eij/t * qqp4(i,j,m) ) * z3**m
           end do
           do l = 1, ncomp
              Iqq3(i,j,l)  = 0.0
              if (parame(l,7) /= 0.0) then
                 do m = 0, 4
                    Iqq3(i,j,l) = Iqq3(i,j,l) + qqp3(i,j,l,m)*z3**m
                 end do
              end if
           end do
        end if
     end do
  end do

  fqq2  = 0.0
  fqq3  = 0.0
  do i = 1, ncomp
     do j = 1, ncomp
        if (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) then
           fqq2 = fqq2 + rhoi(i)*rhoi(j)*xijfa(i,j) * Iqq2(i,j)
           do l=1,ncomp
              if (parame(l,7) /= 0.0) then
                 fqq3 = fqq3 + rhoi(i)*rhoi(j)*rhoi(l)*xijkfa(i,j,l) * Iqq3(i,j,l)
              end if
           end do
        end if
     end do
  end do


  !-----------------------------------------------------------------------------
  ! some first derivatives
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

  do i = 1, ncomp
     do j = 1, ncomp
        Iqq2_r(k,i,j) = 0.0
        if (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) then
           do m = 1, 4
              z3_m = REAL(m) * z3**(m-1) * z3_r(k)
              eij = (parame(i,3)*parame(j,3))**0.5
              Iqq2_r(k,i,j) = Iqq2_r(k,i,j)+ ( qqp2(i,j,m) + eij/t * qqp4(i,j,m) ) * z3_m
           end do
           do o = 1, ncomp
              Iqq3_r(k,i,j,o) = 0.0
              if (parame(o,7) /= 0.0) then
                 do m = 1, 4
                    Iqq3_r(k,i,j,o)=Iqq3_r(k,i,j,o) + qqp3(i,j,o,m)*REAL(m)*z3**(m-1)*z3_r(k)
                 end do
              end if
           end do
        end if
     end do
  end do

  fqq2_r(k) = 0.0
  fqq3_r(k) = 0.0
  do i = 1, ncomp
     fqq2_r(k) = fqq2_r(k) + 2.0*rhoi(i)*xijfa(i,k) * Iqq2(i,k)
     do j = 1, ncomp
        if (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) then
           fqq2_r(k) = fqq2_r(k) + rhoi(i)*rhoi(j)*xijfa(i,j) * Iqq2_r(k,i,j)
           fqq3_r(k) = fqq3_r(k) + 3.0 * rhoi(i)*rhoi(j)*xijkfa(i,j,k) * Iqq3(i,j,k)
           do o = 1, ncomp
              if (parame(o,7) /= 0.0) then
                 fqq3_r(k) =fqq3_r(k) + rhoi(i)*rhoi(j)*rhoi(o)*xijkfa(i,j,o) * Iqq3_r(k,i,j,o)
              end if
           end do
        end if
     end do
  end do

  if (fqq2 < -1.E-50 .AND. fqq3 /= 0.0 .AND. fqq2_r(k) /= 0.0 .AND. fqq3_r(k) /= 0.0) then

     fqq_r(k) = fqq2* (fqq2*fqq2_r(k) - 2.0*fqq3*fqq2_r(k) + fqq2*fqq3_r(k)) / (fqq2-fqq3)**2

  end if

  end do

  !-----------------------------------------------------------------------------
  ! second derivatives
  !-----------------------------------------------------------------------------

  do k = 1, ncomp
  do l = 1, ncomp

  do i = 1, ncomp
     do j = 1, ncomp
        Iqq2_rkrl(i,j) = 0.0
        if (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) then
           do m = 2, 4
              z3_m = REAL(m-1) * REAL(m) * z3**(m-2) * z3_r(k) * z3_r(l)
              eij = (parame(i,3)*parame(j,3))**0.5
              Iqq2_rkrl(i,j) = Iqq2_rkrl(i,j)+ ( qqp2(i,j,m) + eij/t * qqp4(i,j,m) ) * z3_m
           end do
           do o = 1, ncomp
              Iqq3_rkrl(i,j,o) = 0.0
              if (parame(o,7) /= 0.0) then
                 do m = 2, 4
                    Iqq3_rkrl(i,j,o)=Iqq3_rkrl(i,j,o) + qqp3(i,j,o,m)*REAL(m-1)*REAL(m)*z3**(m-2)*z3_r(k)*z3_r(l)
                 end do
              end if
           end do
        end if
     end do
  end do

  fqq2_rkrl = 2.0 * xijfa(k,l) * Iqq2(k,l)
  fqq3_rkrl = 0.0
  do i = 1, ncomp
     fqq2_rkrl = fqq2_rkrl + 2.0 * rhoi(i) * ( xijfa(i,k) * Iqq2_r(l,i,k) +  xijfa(i,l) * Iqq2_r(k,i,l))
     fqq3_rkrl = fqq3_rkrl + 6.0 * rhoi(i) * xijkfa(i,k,l) * Iqq3(i,k,l)
     do j = 1, ncomp
        if (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) then
           fqq2_rkrl = fqq2_rkrl + rhoi(i)*rhoi(j)*xijfa(i,j) * Iqq2_rkrl(i,j)
           fqq3_rkrl = fqq3_rkrl + 3.0 * rhoi(i)*rhoi(j)* ( xijkfa(i,j,k) * Iqq3_r(l,i,j,k)  &
                                                          + xijkfa(i,j,l) * Iqq3_r(k,i,j,l) )
           do o = 1, ncomp
              if (parame(o,7) /= 0.0) then
                 fqq3_rkrl =fqq3_rkrl + rhoi(i)*rhoi(j)*rhoi(o)*xijkfa(i,j,o) * Iqq3_rkrl(i,j,o)
              end if
           end do
        end if
     end do
  end do


  Aqq_rr(k,l) = ( 2.0*fqq2*fqq2_r(l)*( fqq2_r(k) + fqq3_r(k) ) + fqq2*fqq2*( fqq2_rkrl + fqq3_rkrl )  &
                         -2.0*( fqq2_r(k)*fqq2_r(l)*fqq3 + fqq2*fqq3_r(l)*fqq2_r(k) + fqq2*fqq3*fqq2_rkrl ) )  &
                / ( fqq2 - fqq3 )**2  &
                  + fqq_r(k) * 2.0 * ( fqq3_r(l) - fqq2_r(l) ) / ( fqq2 - fqq3 )

  end do
  end do

  deallocate( xijfa, xijkfa, Iqq2, Iqq2_r, Iqq3, Iqq3_r, fqq2_r, fqq3_r, fqq_r, Iqq2_rkrl, Iqq3_rkrl )

end subroutine A_rr_QQ_GROSS

!!$  do i=1,ncomp
!!$     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
!!$     myfac(i) = parame(i,3)/t*parame(i,2)**4 *my2dd(i)
!!$     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
!!$     q_fac(i) = parame(i,3)/t*parame(i,2)**4 *qq2(i)
!!$  end do


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine A_rr_DQ_VRABEC_GROSS( n_comp, Adq_rr )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t, dhs
  USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  integer, intent(in)                   :: n_comp
  real, dimension(n_comp,n_comp), INTENT(IN OUT)  :: Adq_rr

  !-----------------------------------------------------------------------------
  integer                                :: i, j, k, l, m, o

  real                                   :: factor2, factor3, eij
  real                                   :: z3, z3_m
  real                                   :: fdq2, fdq3
  real, allocatable, dimension(:,:)      :: xijfa
  real, allocatable, dimension(:,:,:)    :: xijkfa
  !real                                   :: xijkf_rk

  real, dimension(nc)                    :: z3_r
  real, dimension(nc)                    :: my2dd, myfac, qq2, q_fac
  real, dimension(nc)                    :: rhoi
  real, allocatable, dimension(:,:)      :: Idq2
  real, allocatable, dimension(:,:,:)    :: Idq3
  real, allocatable, dimension(:,:,:)    :: Idq2_r
  real, allocatable, dimension(:,:,:,:)  :: Idq3_r
  real, allocatable, dimension(:)        :: fdq2_r, fdq3_r, fdq_r
  real, allocatable, dimension(:,:)      :: Idq2_rkrl
  real, allocatable, dimension(:,:,:)    :: Idq3_rkrl
  real                                   :: fdq2_rkrl, fdq3_rkrl
  !-----------------------------------------------------------------------------

  z3 = eta
  do i = 1, ncomp
     if ( uij(i,i) == 0.0 ) write (*,*) 'A_rr_DQ_GROSS: do not use dimensionless units'
     if ( uij(i,i) == 0.0 ) stop
     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
     myfac(i) = parame(i,3)/t*parame(i,2)**4 *my2dd(i)
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
     q_fac(i) = parame(i,3)/t*parame(i,2)**4 *qq2(i)
     z3_r(i) = PI/6.0 * mseg(i) * dhs(i)**3
  end do

  !-----------------------------------------------------------------------------
  ! some basic quantities
  !-----------------------------------------------------------------------------

  allocate( xijfa(ncomp,ncomp), xijkfa(ncomp,ncomp,ncomp) )
  allocate( Idq2(ncomp,ncomp), Idq2_r(ncomp,ncomp,ncomp) )
  allocate( Idq3(ncomp,ncomp,ncomp), Idq3_r(ncomp,ncomp,ncomp,ncomp) )
  allocate( fdq2_r(ncomp), fdq3_r(ncomp), fdq_r(ncomp) )
  allocate( Idq2_rkrl(ncomp,ncomp), Idq3_rkrl(ncomp,ncomp,ncomp) )

  rhoi( 1:ncomp ) = x (1:ncomp ) * rho

  factor2 = -9.0/4.0*PI
  factor3 =  PI**2

  do i = 1, ncomp
     do j = 1, ncomp
           xijfa(i,j) = factor2* myfac(i) * q_fac(j) /sig_ij(i,j)**5    ! caution: asymmetric ...(i,j) /= ...(j,i)
           do l = 1, ncomp
              xijkfa(i,j,l)= factor3*( myfac(i)*q_fac(j)*myfac(l) + myfac(i)*q_fac(j)*q_fac(l)*1.193735 )  &
                    / (sig_ij(i,j)*sig_ij(i,l)*sig_ij(j,l))**2
           end do
     end do
  end do

  do i = 1, ncomp
     do j = 1, ncomp
        Idq2(i,j)  = 0.0
        do m = 0, 4
           eij = (parame(i,3)*parame(j,3))**0.5
           Idq2(i,j) = Idq2(i,j) + ( dqp2(i,j,m) + eij/t * dqp4(i,j,m) ) * z3**m
        end do
        do l = 1, ncomp
           Idq3(i,j,l)  = 0.0
           do m = 0, 4
              Idq3(i,j,l) = Idq3(i,j,l) + dqp3(i,j,l,m)*z3**m
           end do
        end do
     end do
  end do

  fdq2  = 0.0
  fdq3  = 0.0
  do i = 1, ncomp
     do j = 1, ncomp
        fdq2 = fdq2 + rhoi(i)*rhoi(j)*xijfa(i,j) * Idq2(i,j)
        do l=1,ncomp
           fdq3 = fdq3 + rhoi(i)*rhoi(j)*rhoi(l)*xijkfa(i,j,l) * Idq3(i,j,l)
        end do
     end do
  end do


  !-----------------------------------------------------------------------------
  ! some first derivatives
  !-----------------------------------------------------------------------------

  do k = 1, ncomp

  do i = 1, ncomp
     do j = 1, ncomp
        Idq2_r(k,i,j) = 0.0
        do m = 1, 4
           z3_m = REAL(m) * z3**(m-1) * z3_r(k)
           eij = (parame(i,3)*parame(j,3))**0.5
           Idq2_r(k,i,j) = Idq2_r(k,i,j)+ ( dqp2(i,j,m) + eij/t * dqp4(i,j,m) ) * z3_m
        end do
        do o = 1, ncomp
           Idq3_r(k,i,j,o) = 0.0
           do m = 1, 4
              Idq3_r(k,i,j,o)=Idq3_r(k,i,j,o) + dqp3(i,j,o,m)*REAL(m)*z3**(m-1)*z3_r(k)
           end do
        end do
     end do
  end do

  fdq2_r(k) = 0.0
  fdq3_r(k) = 0.0
  do i = 1, ncomp
     !xijf_r(k,i) = factor2 * ( myfac(k)*q_fac(i) + myfac(i)*q_fac(k) ) / sig_ij(i,k)**5
     !fdq2_r(k) = fdq2_r(k) +rhoi(i)*xijf_r(k,i) *Idq2(i,k)
     fdq2_r(k) = fdq2_r(k) +rhoi(i)*( xijfa(i,k)+xijfa(k,i) ) *Idq2(i,k)
     do j = 1, ncomp
        fdq2_r(k) = fdq2_r(k) + rhoi(i)*rhoi(j)*xijfa(i,j) * Idq2_r(k,i,j)
        !xijkf_rk=rhoi(i)*rhoi(j)  &
        !     *( myfac(i)*q_fac(j)*myfac(k) + myfac(i)*q_fac(k)*myfac(j) + myfac(k)*q_fac(i)*myfac(j)  &
        !      + (myfac(i)*q_fac(j)*q_fac(k) + myfac(i)*q_fac(k)*q_fac(j)  &
        !      + myfac(k)*q_fac(i)*q_fac(j))*1.193735 ) / (sig_ij(i,j)*sig_ij(i,k)*sig_ij(j,k))**2
        !fdq3_r(k) = fdq3_r(k) + factor3 * xijkf_rk * Idq3(i,j,k)
        fdq3_r(k) = fdq3_r(k) + rhoi(i)*rhoi(j)* (xijkfa(i,j,k)+xijkfa(i,k,j)+xijkfa(k,i,j)) *Idq3(i,j,k)
        do o = 1, ncomp
           fdq3_r(k) =fdq3_r(k) + rhoi(i)*rhoi(j)*rhoi(o)*xijkfa(i,j,o) * Idq3_r(k,i,j,o)
        end do
     end do
  end do

  if (fdq2 < -1.E-50 .AND. fdq3 /= 0.0 .AND. fdq2_r(k) /= 0.0 .AND. fdq3_r(k) /= 0.0) then

     fdq_r(k) = fdq2* (fdq2*fdq2_r(k) - 2.0*fdq3*fdq2_r(k) + fdq2*fdq3_r(k)) / (fdq2-fdq3)**2

  end if

  end do

  !-----------------------------------------------------------------------------
  ! second derivatives
  !-----------------------------------------------------------------------------

  do k = 1, ncomp
  do l = 1, ncomp

  do i = 1, ncomp
     do j = 1, ncomp
        Idq2_rkrl(i,j) = 0.0
        do m = 2, 4
           z3_m = REAL(m-1) * REAL(m) * z3**(m-2) * z3_r(k) * z3_r(l)
           eij = (parame(i,3)*parame(j,3))**0.5
           Idq2_rkrl(i,j) = Idq2_rkrl(i,j)+ ( dqp2(i,j,m) + eij/t * dqp4(i,j,m) ) * z3_m
        end do
        do o = 1, ncomp
           Idq3_rkrl(i,j,o) = 0.0
           do m = 2, 4
              Idq3_rkrl(i,j,o)=Idq3_rkrl(i,j,o) + dqp3(i,j,o,m)*REAL(m-1)*REAL(m)*z3**(m-2)*z3_r(k)*z3_r(l)
           end do
        end do
     end do
  end do


  !fdq2_rkrl = factor2 * ( myfac(k)*q_fac(l) + myfac(l)*q_fac(k) ) /sig_ij(k,l)**5 * Idq2(k,l)
  fdq2_rkrl = ( xijfa(l,k)+xijfa(k,l) ) * Idq2(k,l)
  fdq3_rkrl = 0.0
  do i = 1, ncomp
     ! fdq2_rkrl = fdq2_rkrl + rhoi(i) * ( xijf_r(k,i) * Idq2_r(l,i,k) +  xijf_r(l,i) * Idq2_r(k,i,l))
     fdq2_rkrl = fdq2_rkrl + rhoi(i) * ( ( xijfa(i,k)+xijfa(k,i) ) * Idq2_r(l,i,k)  &
                                      +  ( xijfa(i,l)+xijfa(l,i) ) * Idq2_r(k,i,l) )
     fdq3_rkrl = fdq3_rkrl + rhoi(i) * ( xijkfa(i,l,k)+xijkfa(i,k,l)+xijkfa(l,i,k)  &
                                        +xijkfa(k,i,l)+xijkfa(l,k,i)+xijkfa(k,l,i) ) * Idq3(i,k,l)
     do j = 1, ncomp
        fdq2_rkrl = fdq2_rkrl + rhoi(i)*rhoi(j)*xijfa(i,j) * Idq2_rkrl(i,j)
        fdq3_rkrl = fdq3_rkrl + rhoi(i)*rhoi(j)*( (xijkfa(i,j,k)+xijkfa(i,k,j)+xijkfa(k,i,j)) *Idq3_r(l,i,j,k)  &
                                                + (xijkfa(i,j,l)+xijkfa(i,l,j)+xijkfa(l,i,j)) *Idq3_r(k,i,j,l) )
        do o = 1, ncomp
           fdq3_rkrl =fdq3_rkrl + rhoi(i)*rhoi(j)*rhoi(o)*xijkfa(i,j,o) * Idq3_rkrl(i,j,o)
        end do
     end do
  end do

  Adq_rr(k,l) = ( 2.0*fdq2*fdq2_r(l)*( fdq2_r(k) + fdq3_r(k) ) + fdq2*fdq2*( fdq2_rkrl + fdq3_rkrl )  &
                         -2.0*( fdq2_r(k)*fdq2_r(l)*fdq3 + fdq2*fdq3_r(l)*fdq2_r(k) + fdq2*fdq3*fdq2_rkrl ) )  &
                / ( fdq2 - fdq3 )**2  &
                  + fdq_r(k) * 2.0 * ( fdq3_r(l) - fdq2_r(l) ) / ( fdq2 - fdq3 )

  end do
  end do

  deallocate( xijfa, xijkfa, Idq2, Idq2_r, Idq3, Idq3_r, fdq2_r, fdq3_r, fdq_r, Idq2_rkrl, Idq3_rkrl )

end subroutine A_rr_DQ_VRABEC_GROSS
