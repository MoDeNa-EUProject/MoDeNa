!> \file mod_DFT:DISP_WDA.F90
!!This file contains the subroutines which calculate the contribution of
!!dispersion to the Helmholtz energy functional. 







Module mod_DFT_DISP_WDA


Implicit None

Private

Public :: rhoi_disp_wd
Public :: a_disp_pcsaft
Public :: dF_disp_drho_wda





 Contains



 SUBROUTINE rhoi_disp_wd ( discret, fa_psi, fa_psi_max, psi_j, rhop, rhoi_disp,user )


 Use PetscManagement
 Use basic_variables, ONLY: ncomp, t, parame
 Use mod_DFT, Only: zp,dzp
 IMPLICIT NONE

#include <finclude/petscsys.h>

!
! ----------------------------------------------------------------------
 Type (userctx)      :: user
 PetscScalar         :: rhop(ncomp,user%gxs:user%gxe)
 INTEGER, INTENT(IN)                    :: discret
 INTEGER, INTENT(IN)                    :: fa_psi(ncomp)
 INTEGER, INTENT(IN)                    :: fa_psi_max
! REAL, INTENT(IN)                       :: dzp
 REAL, INTENT(IN)                       :: psi_j(ncomp)
! REAL, INTENT(IN)                       :: zp(user%gxs:user%gxe)
 REAL, INTENT(OUT)                      :: rhoi_disp(user%gxs:user%gxe,ncomp)
! ----------------------------------------------------------------------
 INTEGER                                :: ii, jj, icomp, nn
 REAL                                   :: zmin, zl, zr
 REAL                                   :: int1, zz1, xl, xh
 REAL, DIMENSION(700)                   :: y2, rx, ry1
! ----------------------------------------------------------------------

!write(*,*)'rhoi_wd'

 zmin = 1d-6
 DO icomp = 1, ncomp
   DO ii = (user%xs-fa_psi_max),(user%xe+fa_psi_max)!(-fa_psi_max), (discret+fa_psi_max)
     nn = 0
     zl = zp(ii) - psi_j(icomp)
     zr = zp(ii) + psi_j(icomp)
     DO jj = (ii-fa_psi(icomp)), (ii+fa_psi(icomp))
       IF ( zp(jj+1) > (zl+zmin) ) THEN
         ! first position: left side of the sphere: zl. Linear Interpolation of h's
         IF ( nn == 0 ) THEN
           nn = nn + 1
           rx(1)  = zl
           ry1(1) = 0.0  !  = 0.75*rhop(ii-fa,icomp)*( d.**2 -d.**2 )
!write(*,*)'FIRST',nn,rx(nn),ry1(nn)
         ! middle position: within the sphere: zl < jj < zr
         ELSE
           nn = nn + 1
           zz1 = zp(jj)-zp(ii)       ! distance z12 between 1 and 2
           rx(nn)  = zp(jj)
           ry1(nn) = rhop(icomp,jj) * (psi_j(icomp)*psi_j(icomp) - zz1*zz1)
!write(*,*)'MIDDLE',nn,rx(nn),ry1(nn)
           ! last position: right side of the sphere: zr. Linear Interpolation of h's
           IF ( zp(jj+1) > (zr-zmin) ) THEN
             nn = nn + 1
             rx(nn) = zr
             ry1(nn) = 0.0
!write(*,*)'LAST',nn,rx(nn),ry1(nn)
             EXIT
           END IF
         END IF
       END IF
     END DO
     xl = rx(1)
     xh = rx(nn)

     if( nn >= 700 ) then
       write(*,*) 'rhoi_disp_wd: bigger vectors rx, ry1, ry2, ... required!', nn
       pause
     end if

     CALL spline          ( rx(1:nn), ry1(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral ( rx(1:nn), ry1(1:nn), y2(1:nn), nn, xl, xh, int1 )
     rhoi_disp(ii,icomp) = int1 * 0.75/psi_j(icomp)**3

     if ( rhoi_disp(ii,icomp) < 0.0 ) then
        rhoi_disp(ii,icomp) = 0.0
        do jj = 2, nn
          rhoi_disp(ii,icomp) = rhoi_disp(ii,icomp) + (ry1(jj)+ry1(jj-1))/2.0 *(rx(jj)-rx(jj-1))
        end do
     end if
     if ( rhoi_disp(ii,icomp) < 0.0 ) rhoi_disp(ii,icomp) = rhop(icomp,ii)
   END DO
 END DO


END SUBROUTINE rhoi_disp_wd









!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE a_disp_pcsaft
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE a_disp_pcsaft( discret, fa_psi, fa_psi_max, rhoi_disp, rho_disp, adisp, mydisp, dadisp_dr,user )

 Use PetscManagement
!
 USE parameters, ONLY: PI, np, nc
 USE basic_variables, ONLY: ncomp, t, parame, xi, densta, ensemble_flag
 USE eos_variables, ONLY: dhs, sig_ij, uij
 USE eos_constants, ONLY: ap, bp
 IMPLICIT NONE

#include <finclude/petscsys.h>




!
! ----------------------------------------------------------------------
 Type (userctx)      :: user
 INTEGER, INTENT(IN)                    :: discret
 INTEGER, INTENT(IN)                    :: fa_psi(ncomp)
 INTEGER, INTENT(IN)                    :: fa_psi_max
 REAL, INTENT(IN)                       :: rhoi_disp(user%gxs:user%gxe,ncomp)
 REAL, INTENT(OUT)                      :: rho_disp(user%gxs:user%gxe)
 REAL, INTENT(OUT)                      :: adisp(user%gxs:user%gxe)
 REAL, INTENT(OUT)                      :: mydisp(user%gxs:user%gxe,ncomp)
 REAL, INTENT(OUT)                      :: dadisp_dr(user%gxs:user%gxe,ncomp)
! ----------------------------------------------------------------------
 INTEGER                                :: jj, m
 INTEGER                                :: icomp, jcomp, kcomp
 REAL, DIMENSION(ncomp)                 :: x_disp
 REAL                                   :: eta_disp, zms_disp
 REAL, DIMENSION(0:6)                   :: apar, bpar, apar_rk, bpar_rk
 REAL                                   :: m_mean, C1, I1, I2
 REAL                                   :: r2_ord1, r2_ord2
 REAL                                   :: eta_rk, m_mean_rk, zms2eta
 REAL                                   :: C1_m_mean, C1_eta, C1_rk
 REAL                                   :: I1_rk, I2_rk
 REAL                                   :: r2_ord1_rk, r2_ord2_rk
 REAL                                   :: term1, term2, term3
 INTEGER                                :: iphas
 REAL, DIMENSION(np,nc)                 :: mydisp2
! ----------------------------------------------------------------------

 DO icomp = 1, ncomp
   DO jj = (user%xs-fa_psi_max),(user%xe+fa_psi_max)!(-fa_psi_max), (discret+fa_psi_max)
     m_mean = 0.
     eta_disp = 0.
!      rho_disp(jj) = sum( rhop(jj,1:ncomp) )
     rho_disp(jj) = sum( rhoi_disp(jj,1:ncomp) )
     do kcomp=1, ncomp
!        x_disp(kcomp) = rhop(jj,kcomp) / rho_disp(jj)
       x_disp(kcomp) = rhoi_disp(jj,kcomp) / rho_disp(jj)
       m_mean = m_mean + x_disp(kcomp) * parame(kcomp,1)
!        eta_disp = eta_disp + rhop(jj,kcomp) * parame(kcomp,1) * dhs(kcomp)**3.
       eta_disp = eta_disp + rhoi_disp(jj,kcomp) * parame(kcomp,1) * dhs(kcomp)**3.
     end do
     eta_disp = eta_disp * PI / 6.
     zms_disp = 1. - eta_disp
     if( zms_disp <= 0. ) write(*,'(a56,i6,f12.5)') 'system too dense for disp contribution ( jj, eta_disp ):', jj, eta_disp 
   

     ! quantities of the dispersive free energy contribution
     I1   = 0.
     I2   = 0.
     do m = 0, 6
       apar(m) = ap(m,1) + (1.-1./m_mean)*ap(m,2) + (1.-1./m_mean)*(1.-2./m_mean)*ap(m,3)
       bpar(m) = bp(m,1) + (1.-1./m_mean)*bp(m,2) + (1.-1./m_mean)*(1.-2./m_mean)*bp(m,3)
       I1 = I1 + apar(m) * eta_disp**m
       I2 = I2 + bpar(m) * eta_disp**m
     end do

     r2_ord1 = 0.
     r2_ord2 = 0.
     do kcomp = 1, ncomp
       do jcomp = 1, ncomp
         r2_ord1 = r2_ord1 + rhoi_disp(jj,kcomp)*rhoi_disp(jj,jcomp)*parame(kcomp,1) & 
                          *parame(jcomp,1)*sig_ij(kcomp,jcomp)**3 * uij(kcomp,jcomp)/t
         r2_ord2 = r2_ord2 + rhoi_disp(jj,kcomp)*rhoi_disp(jj,jcomp)*parame(kcomp,1) & 
                          *parame(jcomp,1)*sig_ij(kcomp,jcomp)**3 * (uij(kcomp,jcomp)/t)**2
       end do
     end do

     C1 = (1.-m_mean)*(20.*eta_disp-27.*eta_disp**2 +12.*eta_disp**3 -2.*eta_disp**4 )/(zms_disp*(2.-eta_disp))**2
     C1 = 1. + m_mean*(8.*eta_disp-2.*eta_disp**2)/zms_disp**4 + C1
     C1 = 1. / C1


     ! dispersive free energy contribution
     adisp(jj) = -2.*PI*I1*r2_ord1/rho_disp(jj) - PI*C1*m_mean*I2*r2_ord2/rho_disp(jj)


     ! quantity-derivatives of the dispersive free energy contribution
     m_mean_rk = ( parame(icomp,1) - m_mean ) / rho_disp(jj)

     do m = 0, 6
       apar_rk(m) = m_mean_rk/m_mean**2 * ( ap(m,2) + (3. - 4./m_mean) * ap(m,3) )
       bpar_rk(m) = m_mean_rk/m_mean**2 * ( bp(m,2) + (3. - 4./m_mean) * bp(m,3) )
     end do
     eta_rk = parame(icomp,1) * dhs(icomp)**3. * PI / 6.
     I1_rk = apar_rk(0) + apar(1)*eta_rk + apar_rk(1)*eta_disp
     I2_rk = bpar_rk(0) + bpar(1)*eta_rk + bpar_rk(1)*eta_disp
     do m = 2, 6
       I1_rk = I1_rk + apar(m)*REAL(m)*eta_disp**REAL(m-1)*eta_rk + apar_rk(m)*eta_disp**REAL(m)
       I2_rk = I2_rk + bpar(m)*REAL(m)*eta_disp**REAL(m-1)*eta_rk + bpar_rk(m)*eta_disp**REAL(m)
     end do

     r2_ord1_rk  = 0.
     r2_ord2_rk  = 0.
     do kcomp = 1,ncomp
       r2_ord1_rk = r2_ord1_rk + rhoi_disp(jj,kcomp) * parame(kcomp,1) * sig_ij(icomp,kcomp)**3 * uij(icomp,kcomp)/t
       r2_ord2_rk = r2_ord2_rk + rhoi_disp(jj,kcomp) * parame(kcomp,1) * sig_ij(icomp,kcomp)**3 *(uij(icomp,kcomp)/t)**2 
     end do
     r2_ord1_rk = 2. * parame(icomp,1) * r2_ord1_rk
     r2_ord2_rk = 2. * parame(icomp,1) * r2_ord2_rk

     zms2eta = zms_disp * (2.-eta_disp)
     C1_m_mean = ( 8.*eta_disp - 2.*eta_disp*eta_disp ) / zms_disp**4 &
               - ( 20.*eta_disp - 27.*eta_disp*eta_disp + 12.*eta_disp**3 - 2.*eta_disp**4 ) / zms2eta**2
     C1_eta = m_mean * ( 8. + 20.*eta_disp - 4.*eta_disp*eta_disp ) / zms_disp**5 &
            + (1. - m_mean) * ( 40. - 48.*eta_disp + 12.*eta_disp*eta_disp + 2.*eta_disp**3 ) / zms2eta**3
     C1_rk = - C1 * C1 * ( eta_rk * C1_eta + m_mean_rk * C1_m_mean )


     ! chemical potential and derivative of adisp (analytically)
     term1 = - 2. * PI * ( I1 * r2_ord1_rk + r2_ord1 * I1_rk )
     term2 = - PI * C1 * m_mean * ( I2 * r2_ord2_rk + r2_ord2 * I2_rk )
     term3 = - PI * I2 * r2_ord2 * ( m_mean * C1_rk + C1 * m_mean_rk )
     mydisp(jj,icomp) = term1 + term2 + term3
     dadisp_dr(jj,icomp) = ( mydisp(jj,icomp) - adisp(jj) ) / rho_disp(jj) 


     ! chemical potential and derivative of adisp (numerically)
!      do kcomp=1, ncomp
!        xi(1,kcomp) = x_disp(kcomp)
!      end do
!      iphas = 1
!      ensemble_flag = 'tv'
!      densta(1) = eta_disp
!      call ONLY_ONE_TERM_EOS_NUMERICAL ( 'disp_term', 'PC-SAFT  ' )
!      if( rho_disp(jj) /= 0. ) CALL FUGACITY(mydisp2)
!      call RESTORE_PREVIOUS_EOS_NUMERICAL
!      mydisp(jj,1:ncomp) = mydisp2(iphas,1:ncomp)
!      dadisp_dr(jj,icomp) = ( mydisp(jj,icomp) - adisp(jj) ) / rho_disp(jj)


     if( rho_disp(jj) == 0. ) then
       adisp(jj) = 0.
       mydisp(jj,icomp) = 0.
       dadisp_dr(jj,icomp) = 0.
     end if

   END DO
 END DO


END SUBROUTINE a_disp_pcsaft







!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE dF_disp_drho_wda
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE dF_disp_drho_wda( ii, WDA_var, fa_psi, psi_j, rhop, rhop_sum, &
                         rhoi_disp, rho_disp, adisp, mydisp, dadisp_dr, dF_drho_att,user )

 Use PetscManagement
 Use basic_variables, ONLY: ncomp
 Use mod_DFT, Only: zp
 IMPLICIT NONE

#include <finclude/petscsys.h>

!
! ----------------------------------------------------------------------
 Type (userctx)      :: user
 PetscScalar         :: rhop(ncomp,user%gxs:user%gxe)
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: ii, WDA_var
 INTEGER, INTENT(IN)                    :: fa_psi(ncomp)
 REAL, INTENT(IN)                       :: psi_j(ncomp)
! REAL, INTENT(IN)                       :: zp(user%gxs:user%gxe)
 REAL, INTENT(IN)                       :: adisp(user%gxs:user%gxe)
 REAL, INTENT(IN)                       :: mydisp(user%gxs:user%gxe,ncomp)
 REAL, INTENT(IN)                       :: rhop_sum(user%gxs:user%gxe)
 REAL, INTENT(IN)                       :: rhoi_disp(user%gxs:user%gxe,ncomp)
 REAL, INTENT(IN)                       :: rho_disp(user%gxs:user%gxe)
 REAL, INTENT(IN)                       :: dadisp_dr(user%gxs:user%gxe,ncomp)
 REAL, INTENT(OUT)                      :: dF_drho_att(ncomp)
! ----------------------------------------------------------------------
 INTEGER                                :: jj, icomp, nn
 REAL                                   :: int3
 REAL                                   :: zmin, zl, zr, zz1, xl, xh
 REAL, DIMENSION(700)                   :: y2, hx, hy3
! ----------------------------------------------------------------------

 zmin = 1d-6
 DO icomp = 1, ncomp
   nn = 0
   zl = zp(ii) - psi_j(icomp)
   zr = zp(ii) + psi_j(icomp)
   DO jj = (ii-fa_psi(icomp)), (ii+fa_psi(icomp))
     IF ( zp(jj+1) > (zl+zmin) ) THEN
       ! first position: left side of the sphere: zl. Linear Interpolation of h's
       IF ( nn == 0 ) THEN
         nn = nn + 1
         hx(1) = zl
         hy3(1) = 0.
!write(*,*)'FIRST',nn,hx(nn),hy3(nn)
       ! middle position: within the sphere: zl < jj < zr
       ELSE 
         nn = nn + 1
         zz1 = zp(jj) - zp(ii)                 ! distance z12 between 1 and 2
         hx(nn)  = zp(jj)
         if(WDA_var == 1) hy3(nn) = mydisp(jj,icomp) * ( psi_j(icomp)*psi_j(icomp) - zz1**2 )
         if(WDA_var == 2) hy3(nn) = rhop_sum(jj) * dadisp_dr(jj,icomp) * ( psi_j(icomp)*psi_j(icomp) - zz1**2 )
!write(*,*)'MIDDLE',nn,hx(nn),hy3(nn)
         ! last position: right side of the sphere: zr. Linear Interpolation of h's
         IF ( zp(jj+1) > (zr-zmin) ) THEN
           nn = nn + 1
           hx(nn)  = zr
           hy3(nn) = 0.
!write(*,*)'LAST',nn,hx(nn),hy3(nn)

           EXIT
         END IF
       END IF
     END IF
   END DO
   xl = hx(1)
   xh = hx(nn)

   if( nn >= 700 ) then
     write(*,*) 'Surface Tension Code: dF_disp_wd_pcsaft: bigger vectors hx, hy3, ... required!', nn
     stop 5
   end if
   
   CALL spline         ( hx(1:nn), hy3(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
   CALL splint_integral( hx(1:nn), hy3(1:nn), y2(1:nn), nn, xl, xh, int3 )
   
   if(WDA_var == 1) dF_drho_att(icomp) = int3 * 0.75 / psi_j(icomp)**3.
   if(WDA_var == 2) dF_drho_att(icomp) = int3 * 0.75 / psi_j(icomp)**3. + adisp(ii)
 END DO


END SUBROUTINE dF_disp_drho_wda



End Module mod_DFT_DISP_WDA


! ! ! 
! ! ! Subroutine DISP_Weighted_Densities(rhop, rhop_wd, user) 
! ! ! 
! ! ! Use PetscManagement
! ! ! Use BASIC_VARIABLES, Only: ncomp
! ! ! Use mod_DFT, Only: fa_disp,ab_disp,zp,dzp
! ! ! Implicit None
! ! ! 
! ! ! #include <finclude/petscsys.h>
! ! ! 
! ! ! !passed
! ! ! Type (userctx)    :: user
! ! ! PetscScalar       :: rhop(ncomp,user%gxs:user%gxe)
! ! ! REAL, INTENT(OUT) :: rhop_wd(user%gxs:user%gxe,ncomp)
! ! ! 
! ! ! 
! ! ! !local
! ! ! INTEGER :: i,j,k
! ! ! INTEGER :: n
! ! ! INTEGER, parameter :: NMAX = 800
! ! ! REAL    :: x_int(NMAX),y_int(NMAX),y2(NMAX) !Fehlermeldung einbauen, falls dim > 400!! 
! ! ! REAL    :: zmin,dz,zz
! ! ! REAL    :: xlo,xhi,int1
! ! ! INTEGER :: fa_disp_max
! ! ! 
! ! ! 
! ! ! 
! ! ! zmin = 1d-6
! ! ! fa_disp_max = maxval(fa_disp(1:ncomp))
! ! ! 
! ! ! Do k = 1,ncomp
! ! !    
! ! !    Do i = user%xs - fa_disp_max , user%xe + fa_disp_max
! ! !    n     = 1
! ! !    x_int = 0.0
! ! !    y_int = 0.0
! ! ! 
! ! !      Do j = i-fa_disp(k),i+fa_disp(k)
! ! ! 
! ! !          If( ( zp(i)-zp(j+1) ) < ab_disp(k)  .and. ( zp(i) - zp(j) ) >= ab_disp(k)   ) Then !the position of j+1 is already within i-d/2 while j is still outside this range in this case, the integration steplength (dz) is just the distance, which j+1 overlaps with i-d/2 and what is integrated is the interpolated value of the integrand
! ! ! 
! ! !            If(n/= 1) stop 'n /=1 in DISP_Weighted_Densities' !here always n=1!
! ! !            zz = zp(j) - zp(i)                    !distance between grid points j and i
! ! !            dz = zp(j+1) - (zp(i) - ab_disp(k))           !the part of the intervall between zp(j) and zp(j+1) which is already within i-d/2
! ! !            !if(dz < epsilon(dz)) dz = epsilon(dz) !bei unguenstiger Kombination von sig und ngrid kann dz unter Machinengenauigkeit epsilon liegen, dann ist x(2) = x(1) + dz = x(1) -> das fuehrt zu Abbruch in Spline Interpolation           
! ! !            x_int(n) = 0.0              !array containing x-values for spline integration 
! ! !            y_int(n) = 0.0
! ! !       
! ! !       
! ! !          Else If (zp(j) > (zp(i)-ab_disp(k)) .and. zp(j) <= (zp(i)+ab_disp(k))) Then !grid point j  within i+-d2 
! ! !            
! ! !            n = n + 1
! ! !            x_int(n) = x_int(n-1) +  dz !first time in this If condition, dz is stil the old value from above!
! ! !            zz = zp(j) - zp(i)
! ! !            dz = dzp
! ! !            y_int(n) = rhop(k,j) * (ab_disp(k)*ab_disp(k) - zz*zz  )
! ! !           
! ! !           
! ! !            If (zp(j) < (zp(i)+ab_disp(k)) .and. zp(j+1) >= (zp(i)+ab_disp(k)) ) Then !zp(j) is still within zp(i)+d2 but zp(j+1) is already outside zp(i)+d2
! ! !              dz = zp(i) + ab_disp(k) - zp(j)
! ! !              !If(dz <= epsilon(dz)) exit   !wie oben, kann auch hier bei ungluecklicher Wahl von sig und ngrid dz < eps werden und somit x(n) = x(n-1) -> Abbruch in Spline interpolation. Dann einfach ngrid aendern!
! ! !              zz = zp(j) - zp(i)
! ! !              n = n + 1
! ! !              x_int(n) = x_int(n-1) + dz
! ! !              y_int(n) = 0.0
! ! !         
! ! !         
! ! !            End If           
! ! !           End If    
! ! !      End Do
! ! !      
! ! !      
! ! !      xlo = x_int(1)
! ! !      xhi = x_int(n)     
! ! ! 
! ! !       If(n > NMAX) stop 'Increase NMAX in DISP_Weighted_Densities (auch in AD Routine!!)'
! ! ! 
! ! ! write(*,*)'n',n
! ! ! pause
! ! ! write(*,*)'rx',x_int(1:n)
! ! ! pause
! ! ! write(*,*)'ry',y_int(1:n)
! ! ! pause
! ! ! 
! ! ! 
! ! ! 
! ! !      CALL spline( x_int(1:n), y_int(1:n), n, 1.E30, 1.E30, y2(1:n) )
! ! !      CALL splint_integral ( x_int(1:n), y_int(1:n), y2(1:n), n, xlo, xhi, int1 )
! ! !      rhop_wd(i,k) = 0.75 * int1 / ab_disp(k)**3
! ! !      
! ! !       if ( rhop_wd(i,k) < 0.0 ) then
! ! !         rhop_wd(i,k) = 0.0
! ! !         do j = 2, n
! ! !           rhop_wd(i,k) = rhop_wd(i,k) + (y_int(j)+y_int(j-1))/2.0 *(x_int(j)-x_int(j-1))
! ! !         end do
! ! !      end if
! ! !      if ( rhop_wd(i,k) < 0.0 ) rhop_wd(i,k) = rhop(k,i)
! ! !      
! ! !    End Do
! ! ! End Do
! ! ! 
! ! ! 
! ! ! End Subroutine DISP_Weighted_Densities
! ! ! 
! ! ! 
! ! ! 
! ! ! Subroutine DISP_mu(rhop_wd,f_disp,my_disp,df_disp_drk,user)
! ! ! 
! ! ! Use PetscManagement
! ! ! Use PARAMETERS, Only: PI
! ! ! Use EOS_CONSTANTS, Only: ap,bp
! ! ! Use BASIC_VARIABLES, Only: ncomp,t,parame
! ! ! Use EOS_VARIABLES, Only: dhs,sig_ij,uij
! ! ! Use mod_DFT, Only: fa_disp
! ! ! Implicit None
! ! ! 
! ! ! #include <finclude/petscsys.h>
! ! ! 
! ! ! !passed
! ! ! Type (userctx)    :: user
! ! ! REAL, INTENT(IN)  :: rhop_wd(user%gxs:user%gxe,ncomp)
! ! ! REAL, INTENT(OUT) :: my_disp(user%gxs:user%gxe,ncomp) !chemPot_disp / kT
! ! ! REAL, INTENT(OUT) :: f_disp(user%gxs:user%gxe)  !F_disp / NkT
! ! ! REAL, INTENT(OUT) :: df_disp_drk(user%gxs:user%gxe,ncomp) ! d(F/NkT) / drho_k = (mu/kT - F_disp / NkT)/rho                    
! ! ! 
! ! ! ! ! !local
! ! ! ! ! INTEGER :: k,ii,kk,m
! ! ! ! ! REAL :: m_mean, m_rk(ncomp)
! ! ! ! ! REAL :: apar(0:6),bpar(0:6)
! ! ! ! ! REAL :: ap_rk(ncomp,0:6),bp_rk(ncomp,0:6)
! ! ! ! ! REAL :: xi(ncomp),z3
! ! ! ! ! REAL :: I1,I2,I1_rk,I2_rk
! ! ! ! ! REAL :: order1,order2,ord1_rk,ord2_rk
! ! ! ! ! REAL :: c1_con,c2_con, c1_rk,rho2
! ! ! ! ! REAL :: rhop_wd_sum, eta_disp, eta_rk, zms
! ! ! ! ! INTEGER :: fa_disp_max
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! fa_disp_max = maxval(fa_disp(1:ncomp))
! ! ! ! ! 
! ! ! ! ! Do k = 1,ncomp
! ! ! ! !   Do ii = user%xs-fa_disp_max,user%xe+fa_disp_max
! ! ! ! !      
! ! ! ! !     rhop_wd_sum = SUM(rhop_wd(ii,1:ncomp))
! ! ! ! !     
! ! ! ! !     m_mean   = 0.0
! ! ! ! !     eta_disp = 0.0
! ! ! ! !     Do kk = 1,ncomp
! ! ! ! !        xi(kk) = rhop_wd(ii,kk) / rhop_wd_sum
! ! ! ! !        m_mean = m_mean + xi(kk)*parame(kk,1)
! ! ! ! !        eta_disp = eta_disp + rhop_wd(ii,kk)*parame(kk,1)*dhs(kk)**3
! ! ! ! !     End Do
! ! ! ! !     
! ! ! ! !     eta_disp = eta_disp * PI / 6.0
! ! ! ! !     eta_rk   = parame(k,1) * dhs(k)**3 * PI / 6.0
! ! ! ! !        
! ! ! ! !     m_rk(k) = ( parame(k,1) - m_mean ) / rhop_wd_sum
! ! ! ! !     
! ! ! ! !     
! ! ! ! !     DO m = 0, 6
! ! ! ! !          apar(m) = ap(m,1) + (1.0-1.0/m_mean)*ap(m,2)  &
! ! ! ! !             + (1.0-1.0/m_mean)*(1.0-2.0/m_mean)*ap(m,3)
! ! ! ! !          bpar(m) = bp(m,1) + (1.0-1.0/m_mean)*bp(m,2)  &
! ! ! ! !             + (1.0-1.0/m_mean)*(1.0-2.0/m_mean)*bp(m,3)
! ! ! ! ! 
! ! ! ! !         ! --- derivatives of apar, bpar to rho_k ---------------------------
! ! ! ! !         ap_rk(k,m) = m_rk(k)/m_mean**2 * ( ap(m,2) + (3.0 -4.0/m_mean) *ap(m,3) )
! ! ! ! !         bp_rk(k,m) = m_rk(k)/m_mean**2 * ( bp(m,2) + (3.0 -4.0/m_mean) *bp(m,3) )
! ! ! ! !     END DO
! ! ! ! ! 
! ! ! ! !     
! ! ! ! !     I1    = 0.0
! ! ! ! !     I2    = 0.0
! ! ! ! !     I1_rk = 0.0
! ! ! ! !     I2_rk = 0.0
! ! ! ! !     DO m = 0, 6
! ! ! ! !       I1  = I1 + apar(m)*eta_disp**REAL(m)
! ! ! ! !       I2  = I2 + bpar(m)*eta_disp**REAL(m)
! ! ! ! !       I1_rk = I1_rk + apar(m)*REAL(m)*eta_disp**REAL(m-1)*eta_rk + ap_rk(k,m)*eta_disp**REAL(m)
! ! ! ! !       I2_rk = I2_rk + bpar(m)*REAL(m)*eta_disp**REAL(m-1)*eta_rk + bp_rk(k,m)*eta_disp**REAL(m)
! ! ! ! !     END DO
! ! ! ! !  
! ! ! ! !  
! ! ! ! !     ord1_rk  = 0.0
! ! ! ! !     ord2_rk  = 0.0
! ! ! ! !     order1   = 0.0
! ! ! ! !     order2   = 0.0
! ! ! ! !     DO kk = 1,ncomp
! ! ! ! !       !sig_ij(kk,k) = 0.5 * ( dhs(kk) + dhs(k) )
! ! ! ! !       !uij(kk,k)    = (1.0 - kij(kk,k)) * SQRT( eps(kk) * eps(k) )
! ! ! ! !       order1 = order1 + xi(kk)*xi(k)* parame(kk,1)*parame(k,1)*sig_ij(kk,k)**3 * uij(kk,k)/t
! ! ! ! !       order2 = order2 + xi(kk)*xi(k)* parame(kk,1)*parame(k,1)*sig_ij(kk,k)**3 * (uij(kk,k)/t)**2
! ! ! ! !       ord1_rk = ord1_rk + 2.0*parame(k,1)*rhop_wd_sum*xi(kk)*parame(kk,1)*sig_ij(kk,k)**3  *uij(kk,k)/t
! ! ! ! !       ord2_rk = ord2_rk + 2.0*parame(k,1)*rhop_wd_sum*xi(kk)*parame(kk,1)*sig_ij(kk,k)**3 *(uij(kk,k)/t)**2 
! ! ! ! !     END DO
! ! ! ! !     
! ! ! ! !     z3  = eta_disp
! ! ! ! !     zms = 1.0 - z3
! ! ! ! !     
! ! ! ! !     c1_con= 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3*z3)/zms**4   &
! ! ! ! !                     + (1.0 - m_mean)*(20.0*z3-27.0*z3*z3 +12.0*z3**3 -2.0*z3**4 )  &
! ! ! ! !                     /(zms*(2.0-z3))**2  )
! ! ! ! !     c2_con= - c1_con*c1_con *(  m_mean*(-4.0*z3*z3+20.0*z3+8.0)/zms**5   &
! ! ! ! !                                 + (1.0 - m_mean) *(2.0*z3**3 +12.0*z3*z3-48.0*z3+40.0)  &
! ! ! ! !                                   /(zms*(2.0-z3))**3  )
! ! ! ! !     c1_rk= c2_con*eta_rk - c1_con*c1_con*m_rk(k)   *  ( (8.0*z3-2.0*z3*z3)/zms**4   &
! ! ! ! !            - (-2.0*z3**4 +12.0*z3**3 -27.0*z3*z3+20.0*z3) / (zms*(2.0-z3))**2  )
! ! ! ! !     
! ! ! ! !     
! ! ! ! !     rho2 = rhop_wd_sum * rhop_wd_sum 
! ! ! ! !     
! ! ! ! !     my_disp(ii,k) = -2.0*PI* ( order1*rho2*I1_rk + ord1_rk*I1 )  &
! ! ! ! !               -    PI* c1_con*m_mean * ( order2*rho2*I2_rk + ord2_rk*I2 )  &
! ! ! ! !               -    PI* ( c1_con*m_rk(k) + c1_rk*m_mean ) * order2*rho2*I2
! ! ! ! ! 
! ! ! ! !     f_disp(ii)  = -2.0 * PI * rhop_wd_sum * I1 * order1 &             !hier * rho_wd_sum, bei Elmar/rho_wd_sum, da in order1 bei Elmar mit rhoi, hier mit xi!!
! ! ! ! !                   -PI * rhop_wd_sum * c1_con * m_mean * I2 * order2
! ! ! ! !               
! ! ! ! !     df_disp_drk(ii,k) = ( my_disp(ii,k) - f_disp(ii) ) / rhop_wd_sum          
! ! ! ! !               
! ! ! ! !    End Do
! ! ! ! !               
! ! ! ! ! End Do
! ! ! 
! ! ! 
! ! ! End Subroutine DISP_mu
! ! ! 
! ! ! 
! ! ! 
! ! ! Subroutine DISP_dFdrho_wda(ii,rhop,rhop_wd,my_disp,f_disp,df_disp_drk,dF_drho_disp,user)
! ! ! 
! ! ! Use PetscManagement
! ! ! 
! ! ! Use BASIC_VARIABLES, Only: ncomp
! ! ! Use mod_DFT, Only: zp,dzp,fa_disp,ab_disp
! ! ! 
! ! ! Implicit None
! ! ! #include <finclude/petscsys.h>
! ! ! 
! ! ! !passed
! ! ! INTEGER, INTENT(IN) :: ii
! ! ! Type (userctx)      :: user
! ! ! PetscScalar         :: rhop(ncomp,user%gxs:user%gxe)
! ! ! REAL, INTENT(IN)    :: rhop_wd(user%gxs:user%gxe,ncomp)
! ! ! REAL, INTENT(IN)    :: my_disp(user%gxs:user%gxe,ncomp) !chemPot_disp / kT
! ! ! REAL, INTENT(IN)    :: f_disp(user%gxs:user%gxe)  !F_disp / NkT
! ! ! REAL, INTENT(IN)    :: df_disp_drk(user%gxs:user%gxe,ncomp) ! d(F/NkT) / drho_k = (mu/kT - F_disp / NkT)/rho                    
! ! ! REAL, INTENT(OUT)   :: dF_drho_disp(ncomp)
! ! ! 
! ! ! ! ! !local
! ! ! ! ! INTEGER :: n,icomp,jj
! ! ! ! ! REAL    :: zmin,zz,dz
! ! ! ! ! INTEGER, parameter :: NMAX = 800
! ! ! ! ! REAL    :: x_int(NMAX), y_int(NMAX), y2(NMAX)
! ! ! ! ! REAL    :: xhi,xlo,int2
! ! ! ! ! REAL    :: rhop_sum
! ! ! ! ! 
! ! ! ! !  zmin = 1d-6
! ! ! ! !  DO icomp = 1, ncomp
! ! ! ! !    n = 1
! ! ! ! !    x_int = 0.0
! ! ! ! !    y_int = 0.0
! ! ! ! !     
! ! ! ! !    DO jj = (ii-fa_disp(icomp)), (ii+fa_disp(icomp))
! ! ! ! ! 
! ! ! ! !    !IF ( zp(jj+1) > (zl+zmin) ) THEN
! ! ! ! !    If( ( zp(ii)-zp(jj+1) ) < ab_disp(icomp)  .and. ( zp(ii) - zp(jj) ) >= ab_disp(icomp)   ) Then !the position of j+1 is already within i-d/2 while j is still outside this range in this case, the integration steplength (dz) is just the distance, which j+1 overlaps with i-d/2 and what is integrated is the interpolated value of the integrand
! ! ! ! ! 
! ! ! ! !            If(n/= 1) stop 'n /=1 in DISP_Weighted_Densities' !here always n=1!
! ! ! ! !            zz = zp(jj) - zp(ii)                    !distance between grid points j and i
! ! ! ! !            dz = zp(jj+1) - (zp(ii) - ab_disp(icomp))           !the part of the intervall between zp(j) and zp(j+1) which is already within i-d/2
! ! ! ! !            !if(dz < epsilon(dz)) dz = epsilon(dz) !bei unguenstiger Kombination von sig und ngrid kann dz unter Machinengenauigkeit epsilon liegen, dann ist x(2) = x(1) + dz = x(1) -> das fuehrt zu Abbruch in Spline Interpolation           
! ! ! ! !            x_int(n) = 0.0              !array containing x-values for spline integration 
! ! ! ! !            y_int(n) = 0.0
! ! ! ! !   
! ! ! ! !     Else If (zp(jj) > (zp(ii)-ab_disp(icomp)) .and. zp(jj) <= (zp(ii)+ab_disp(icomp))) Then !grid point j  within i+-d2 
! ! ! ! !            
! ! ! ! !            n = n + 1
! ! ! ! !            x_int(n) = x_int(n-1) +  dz !first time in this If condition, dz is stil the old value from above!
! ! ! ! !            zz = zp(jj) - zp(ii)
! ! ! ! !            dz = dzp
! ! ! ! !            y_int(n) = my_disp(jj,icomp) * (ab_disp(icomp)*ab_disp(icomp) - zz*zz  )
! ! ! ! !       
! ! ! ! !            !rhop_sum = sum(rhop(1:ncomp,jj))
! ! ! ! !            !y_int(n) = rhop_sum * df_disp_drk(jj,icomp) * (ab_disp(icomp)*ab_disp(icomp) - zz*zz) 
! ! ! ! !       
! ! ! ! !            If (zp(jj) < (zp(ii)+ab_disp(icomp)) .and. zp(jj+1) >= (zp(ii)+ab_disp(icomp)) ) Then !zp(j) is still within zp(i)+d2 but zp(j+1) is already outside zp(i)+d2
! ! ! ! !              dz = zp(ii) + ab_disp(icomp) - zp(jj)
! ! ! ! !              !If(dz <= epsilon(dz)) exit   !wie oben, kann auch hier bei ungluecklicher Wahl von sig und ngrid dz < eps werden und somit x(n) = x(n-1) -> Abbruch in Spline interpolation. Dann einfach ngrid aendern!
! ! ! ! !              zz = zp(jj) - zp(ii)
! ! ! ! !              n = n + 1
! ! ! ! !              x_int(n) = x_int(n-1) + dz
! ! ! ! !              y_int(n) = 0.0
! ! ! ! !                        
! ! ! ! !            End If       
! ! ! ! !     End If             
! ! ! ! !    END DO
! ! ! ! !    xlo = x_int(1)
! ! ! ! !    xhi = x_int(n)
! ! ! ! !    
! ! ! ! !    If(n > NMAX) stop 'Increase NMAX in DISP_dFdrho_wda  (auch in AD Routine!!)'
! ! ! ! ! 
! ! ! ! !    CALL spline         ( x_int(1:n), y_int(1:n), n, 1.E30, 1.E30, y2(1:n) )
! ! ! ! !    CALL splint_integral( x_int(1:n), y_int(1:n), y2(1:n), n, xlo, xhi, int2 )
! ! ! ! !    
! ! ! ! !    dF_drho_disp(icomp) = int2 * 0.75 / ab_disp(icomp)**3
! ! ! ! !    !dF_drho_disp(icomp) = int2 * 0.75 / ab_disp(icomp)**3 + f_disp(ii)
! ! ! ! !  END DO
! ! ! ! ! 
! ! ! 
! ! ! End Subroutine DISP_dFdrho_wda
! ! ! 

