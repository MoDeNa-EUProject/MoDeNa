!> \file mod_DFT_CHAIN.F90
!!This file contains the subroutines which calculate the contribution of
!!chain formation to the Helmholtz energy functional. 


Module mod_DFT_CHAIN

Implicit None
private


public :: Chain_aux
public :: Chain_dFdrho


 contains
 

 
Subroutine Chain_aux(rhop,rhobar,lambda,user)

Use BASIC_VARIABLES, Only: ncomp
Use EOS_VARIABLES, Only: dhs,rho
Use mod_DFT, Only: zp,dzp,fa


!PETSc module
Use PetscManagement

#include <finclude/petscsys.h>

!passed
type (userctx)   :: user
PetscScalar      :: rhop(ncomp,user%gxs:user%gxe)
REAL,INTENT(OUT) :: rhobar(user%gxs:user%gxe,ncomp)
REAL,INTENT(OUT) :: lambda(user%gxs:user%gxe,ncomp)

!local
INTEGER :: i,j,k
REAL    :: dhsk
INTEGER :: fak,n
REAL    :: zz,dz,xlo,xhi,integral_lamb,integral_rb
INTEGER, parameter :: NMAX = 800
REAL,dimension(NMAX)    :: x_int, lamb_int, rb_int 
REAL,dimension(NMAX)    :: y2_lamb, y2_rb 
REAL    :: rhopjk,rhopjp1k


!fak = maxval(fa(1:ncomp))

Do k = 1,ncomp
   dhsk = dhs(k)
   fak  = fa(k) + 1

   Do i = user%xs-fak,user%xe+fak  !lambda und rhobar werden bis i+-sig gebraucht -> Schleife bis +- fa
      n = 1        !this is the index of the arrays that will be passed to the spline integration routines
      x_int    = 0.0
      lamb_int = 0.0
      rb_int   = 0.0
   
      Do j = i-fak,i+fak      !um lambda bei i zu berechnen, muss bis +- sig um i integriert werden -> Schleife bis +- fa 
       rhopjk   = rhop(k,j)
       rhopjp1k = rhop(k,j+1)
   

       If( ( zp(i)-zp(j+1) ) < dhsk  .and. ( zp(i) - zp(j) ) >= dhsk   ) Then !the position of j+1 is already within i-d while j is still outside this range in this case, the integration steplength (dz) is just the distance, which j+1 overlaps with i-d and what is integrated is the interpolated value of the integrand

           If(n/= 1) Then
              write(*,*) 'Surface Tension Code: n /=1 in Chain_aux' !here always n=1!
              stop 5
           End If
              
           zz = zp(j) - zp(i)                    !distance between grid points j and i
           dz = zp(j+1) - (zp(i) - dhsk)           !the part of the intervall between zp(j) and zp(j+1) which is already within i-d
           !if(dz < epsilon(dz)) dz = epsilon(dz) !bei unguenstiger Kombination von sig und ngrid kann dz unter Machinengenauigkeit epsilon liegen, dann ist x(2) = x(1) + dz = x(1) -> das fuehrt zu Abbruch in Spline Interpolation           
           x_int(n)    = 0.0              !array containing x-values for spline integration 
           
           lamb_int(n) = rhopjk + (rhopjp1k - rhopjk)/dzp * (dzp-dz)  !lineare interpolation analog zum FMT Teil
           rb_int(n)   = 0.0 !erklärung analog wie bei n3_int in FMT Teil
          
          
       Else If (zp(j) > (zp(i)-dhsk) .and. zp(j) <= (zp(i)+dhsk)) Then !grid point j  within i+-d 
           
           n = n + 1
           x_int(n) = x_int(n-1) +  dz !first time in this If condition, dz is stil the old value from above!
           zz = zp(j) - zp(i)
           dz = dzp
           lamb_int(n) = rhopjk 
           rb_int(n)   = rhopjk * ( dhsk**2 - zz**2 )
           
           If (zp(j) < (zp(i)+dhsk) .and. zp(j+1) >= (zp(i)+dhsk) ) Then !zp(j) is still within zp(i)+d but zp(j+1) is already outside zp(i)+d

           dz = zp(i) + dhsk - zp(j)
           !If(dz <= epsilon(dz)) exit   !wie oben, kann auch hier bei ungluecklicher Wahl von sig und ngrid dz < eps werden und somit x(n) = x(n-1) -> Abbruch in Spline interpolation. Dann einfach ngrid aendern!
           zz = zp(j) - zp(i)
           n = n + 1
           x_int(n)    = x_int(n-1) + dz
           rb_int(n)   = 0.0
           lamb_int(n) = rhopjk + (rhopjp1k - rhopjk)/dzp * dz
!            If(x_int(n) == x_int(n-1)) Then
!              n = n - 1
!              exit
!            End If  

           End If           
       End If
      End Do
      
      
      
    !spline integration
    xlo = x_int(1)
    xhi = x_int(n)

    If(n > NMAX) Then
        write(*,*)'Surface Tension Code: Increase NMAX in Chain_aux (also in AD Routine!!)'
        stop 5
    End If 
    
    call spline ( x_int, lamb_int, n, 1.E30, 1.E30, y2_lamb )
    call spline ( x_int, rb_int, n, 1.E30, 1.E30, y2_rb )
    
    call splint_integral ( x_int, lamb_int, y2_lamb, n, xlo, xhi, integral_lamb )
    call splint_integral ( x_int, rb_int, y2_rb, n, xlo, xhi, integral_rb )
    lambda(i,k) = 0.5  * integral_lamb / dhsk
    rhobar(i,k) = 0.75 * integral_rb / dhsk**3
    
    If(lambda(i,k) < epsilon(dz) ) lambda(i,k) = epsilon(dz) 
    !If ( lambda(i,k) < 0.5*xx(k)*rho ) write (*,*) 'warning: lambda too low',i,lambda(i,k)
    !If ( k == 1 .AND. lambda(i,k) < 0.5*rhob(2,k) ) lambda(i,k) = 0.5*rhob(2,k)

    
   End Do
   
End Do


End Subroutine Chain_aux








Subroutine Chain_dFdrho(i,rhop,lambda,rhobar,dF_drho_CHAIN,f_ch,user)

Use BASIC_VARIABLES, Only: ncomp,parame
Use EOS_VARIABLES, Only: dhs
Use mod_DFT, Only: zp,dzp,fa

!PETSc module
Use PetscManagement

#include <finclude/petscsys.h>

!passed
INTEGER, INTENT(IN) :: i !the grid point to calculate dFdrho at
type (userctx)   :: user
PetscScalar      :: rhop(ncomp,user%gxs:user%gxe)
REAL,INTENT(IN) :: rhobar(user%gxs:user%gxe,ncomp)
REAL,INTENT(IN) :: lambda(user%gxs:user%gxe,ncomp)
REAL,INTENT(OUT) :: dF_drho_CHAIN(user%xs:user%xe,ncomp)
REAL, INTENT(OUT) :: f_ch


!local
INTEGER :: j,k,n
REAL    :: dhsk
INTEGER :: fak
REAL    :: rhopjk,rhopjp1k,logrho,xlo,xhi
REAL    :: ycorr(ncomp),dlny(ncomp,ncomp)
INTEGER, parameter :: NMAX = 800
REAL,dimension(NMAX)    :: x_int, int_1,int_2 !Fehlerwarnung falls 800 ueberschritten einbauen
REAL,dimension(NMAX)    :: y2_1, y2_2 
REAL    :: dz,zz,integral_1,integral_2
REAL    :: lamy


 call Cavity_mix(rhobar(i,1:ncomp),ycorr,dlny)
 
 f_ch = 0.0 

 Do k = 1,ncomp
   fak  = fa(k)
   dhsk = dhs(k)
   n = 1
   x_int = 0.0
   int_1 = 0.0
   int_2 = 0.0

   
   Do j = i-fak,i+fak      !es muss bis +- sig um i integriert werden 
       rhopjk   = rhop(k,j)
       rhopjp1k = rhop(k,j+1)
       
       If( ( zp(i)-zp(j+1) ) < dhsk  .and. ( zp(i) - zp(j) ) >= dhsk   ) Then !the position of j+1 is already within i-d while j is still outside this range in this case, the integration steplength (dz) is just the distance, which j+1 overlaps with i-d and what is integrated is the interpolated value of the integrand

           If(n/= 1) Then
              write(*,*) 'Surface Tension Code: n /=1 in Chain_dFdrho' !here always n=1!
              stop 5 
           End If
           
           zz = zp(j) - zp(i)                    !distance between grid points j and i
           dz = zp(j+1) - (zp(i) - dhsk)           !the part of the intervall between zp(j) and zp(j+1) which is already within i-d
           !if(dz < epsilon(dz)) dz = epsilon(dz) !bei unguenstiger Kombination von sig und ngrid kann dz unter Machinengenauigkeit epsilon liegen, dann ist x(2) = x(1) + dz = x(1) -> das fuehrt zu Abbruch in Spline Interpolation           
           x_int(n) = 0.0              !array containing x-values for spline integration 
           int_1(n) = 0.0 !erklärung analog wie bei n3_int in FMT Teil
           int_2(n) = rhopjk*lambda(j,k) + (rhopjp1k*lambda(j+1,k) - rhopjk*lambda(j,k) )/dzp * (dzp-dz)  !lineare interpolation analog zum FMT Teil
          
       Else If (zp(j) > (zp(i)-dhsk) .and. zp(j) <= (zp(i)+dhsk)) Then !grid point j  within i+-d 
           
           n = n + 1
           x_int(n) = x_int(n-1) +  dz !first time in this If condition, dz is stil the old value from above!
           zz = zp(j) - zp(i)
           dz = dzp
           int_1(n) =  SUM( ( parame(1:ncomp,1)-1.0 ) * rhop(1:ncomp,j)*dlny(1:ncomp,k) ) &
                  * 0.75/dhsk**3 * (dhsk**2-zz**2)
           int_2(n) = rhopjk / lambda(j,k)
           
           If (zp(j) < (zp(i)+dhsk) .and. zp(j+1) >= (zp(i)+dhsk) ) Then !zp(j) is still within zp(i)+d but zp(j+1) is already outside zp(i)+d

           dz = zp(i) + dhsk - zp(j)
           !If(dz <= epsilon(dz)) exit   !wie oben, kann auch hier bei ungluecklicher Wahl von sig und ngrid dz < eps werden und somit x(n) = x(n-1) -> Abbruch in Spline interpolation. Dann einfach ngrid aendern!
           zz = zp(j) - zp(i)
           n = n + 1
           x_int(n)    = x_int(n-1) + dz
           int_1(n) = 0.0
           int_2(n) = rhopjk*lambda(j,k) + (rhopjp1k*lambda(j+1,k) - rhopjk*lambda(j,k) )/dzp * dz 
           
!            If(x_int(n) == x_int(n-1)) Then
!              n = n - 1
!              exit
!            End If  

           End If           
       End If
      End Do
      
      !Spline Integration
      xlo = x_int(1)
      xhi = x_int(n)   

      If(n > NMAX) Then
         write(*,*)'Surface Tension Code: Increase NMAX in Chain_dFdrho (also in AD Routine!!)'
         stop 5
      End If 
 
      CALL spline         ( x_int, int_1, n, 1.E30, 1.E30, y2_1 )
      CALL splint_integral( x_int, int_1, y2_1, n, xlo, xhi, integral_1 )
      CALL spline         ( x_int, int_2, n, 1.E30, 1.E30, y2_2 )
      CALL splint_integral( x_int, int_2, y2_2, n, xlo, xhi, integral_2 )

      If(rhop(k,i) < epsilon(dz))  rhop = epsilon(dz)
      
      dF_drho_CHAIN(i,k) = (parame(k,1) - 1.0)*log(rhop(k,i)) &
                           - (parame(k,1) - 1.0) * ( log(ycorr(k)*lambda(i,k))-1.0 + 0.5*integral_2/dhsk )  - integral_1


      f_ch = f_ch + ( parame(k,1) - 1.0 ) * rhop(k,i) * ( log(rhop(k,i))   - 1.0 )  &
             - ( parame(k,1) - 1.0 ) * rhop(k,i) * ( log(ycorr(k)*lambda(i,k)) - 1.0 )




       
End Do












End Subroutine Chain_dFdrho
















!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE Cavity_mix ( rhoi, ycorr, dlnydr )
!
 USE PARAMETERS, ONLY: PI
 USE BASIC_VARIABLES, ONLY: ncomp,parame
 USE EOS_VARIABLES, Only: dhs
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN)                       :: rhoi(ncomp)
 REAL, INTENT(OUT)                      :: ycorr(ncomp)
 REAL, INTENT(OUT)                      :: dlnydr(ncomp,ncomp)     ! this is: d( ln( yij ) ) / d( rho(k) ) used with i=j
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k
 REAL                                   :: z0, z1, z2, z3, zms, z1_rk, z2_rk, z3_rk
 REAL, DIMENSION(ncomp,ncomp)                 :: dij_ab, gij, gij_rk
! ----------------------------------------------------------------------

 z0 = PI / 6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) )
 z1 = PI / 6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp) )
 z2 = PI / 6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**2 )
 z3 = PI / 6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )

 zms    = 1.0 - z3

 DO i = 1,ncomp
    DO j=1,ncomp
       dij_ab(i,j)=dhs(i)*dhs(j)/(dhs(i)+dhs(j))
    ENDDO
 END DO

 DO k = 1, ncomp
    DO i = 1, ncomp
       z1_rk = PI/6.0 * parame(k,1) * dhs(k)
       z2_rk = PI/6.0 * parame(k,1) * dhs(k)*dhs(k)
       z3_rk = PI/6.0 * parame(k,1) * dhs(k)**3 
       !DO j = 1, ncomp
       j = i
       gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms  &
            + 2.0*(dij_ab(i,j)*z2)**2 /zms**3 
       !dgijdz(i,j)= 1.0/zms/zms + 3.0*dij_ab(i,j)*z2*(1.0+z3)/z3/zms**3   &
       !           + (dij_ab(i,j)*z2/zms/zms)**2 *(4.0+2.0*z3)/z3
       gij_rk(i,j) = z3_rk/zms/zms  &
            + 3.0*dij_ab(i,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
            + dij_ab(i,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
       !END DO

       ycorr(i)  = gij(i,i)
       dlnydr(i,k) = gij_rk(i,i) / gij(i,i)

    END DO
 END DO

END SUBROUTINE Cavity_mix










End Module mod_DFT_CHAIN