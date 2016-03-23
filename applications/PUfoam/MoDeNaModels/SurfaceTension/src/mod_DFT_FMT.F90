!> \file mod_DFT_FMT.F90
!!This file contains the subroutines which calculate the contribution of
!!volume exclusion to the Helmholtz energy functional. 




Module mod_DFT_FMT

Implicit None

Private

Public :: FMT_Weighted_Densities
Public :: FMT_dFdrho

 Contains





Subroutine FMT_Weighted_Densities(rhop,n0,n1,n2,n3,nv1,nv2,phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2,user)

Use PARAMETERS, ONLY: PI
Use BASIC_VARIABLES, Only: ncomp,parame
Use EOS_VARIABLES, Only: dhs
Use mod_DFT, Only: zp,dzp,fa

!PETSc module
Use PetscManagement

#include <finclude/petscsys.h>

!passed
Type (userctx)                                         :: user
PetscScalar                                            :: rhop(ncomp,user%gxs:user%gxe)
REAL,dimension(user%gxs:user%gxe),Intent(OUT)          :: n0,n1,n2,n3,nv1,nv2  !ngp muss groesser als fa+fa/2 sein!!
REAL,dimension(user%gxs:user%gxe),Intent(OUT)          :: phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2

!local
Integer :: k,i,j
Integer :: fa2
REAL    :: dz,d2,zz
INTEGER :: n
INTEGER, parameter :: NMAX = 800
REAL,dimension(NMAX)    :: x_int,n2_int,n3_int,nv2_int !Fehlerwarnung falls 200 ueberschritten einbauen
REAL,dimension(NMAX)    :: y2_n2, y2_n3,y2_nv2 
REAL                   :: int_n2,int_n3,int_nv2,xhi,xlo
REAL                   :: zms,zms2,zms3,logzms 
REAL                   :: nn0,nn1,nn2,nn3,nnv1,nnv2
REAL                   :: rhopjk, rhopjp1k


 n0 = 0.0
 n1 = 0.0
 n2 = 0.0
 n3 = 0.0
 nv1 = 0.0
 nv2 = 0.0

 fa2 = maxval(fa(1:ncomp) + 5) / 2
 
Do k=1,ncomp
 !fa2 = (fa(k) + 5) / 2   !grid points in sig/2
 d2  = dhs(k) / 2.0      !half of dhs [A] 

 
 Do i=user%xs-fa2,user%xe+fa2 !to evaluate dF/drho at any point, the derivatives dphi/dn have to be evaluated at +-d/2 around this point
    n = 1        !this is the index of the arrays that will be passed to the spline integration routines
    x_int   = 0.0
    n2_int  = 0.0
    n3_int  = 0.0 
    nv2_int = 0.0

    Do j=i-fa2,i+fa2 !to evaluate dphi/dn at a given point, the weighted densities are needed at +-d/2 around this point

       rhopjk   = rhop(k,j)   * parame(k,1)
       rhopjp1k = rhop(k,j+1) * parame(k,1)
    
       If( ( zp(i)-zp(j+1) ) < d2  .and. ( zp(i) - zp(j) ) >= d2   ) Then !the position of j+1 is already within i-d/2 while j is still outside this range in this case, the integration steplength (dz) is just the distance, which j+1 overlaps with i-d/2 and what is integrated is the interpolated value of the integrand

           If(n/= 1) Then 
               write(*,*)'Surface Tension Code: n /=1 in FMT_Weighted_Densities' !here always n=1!
               stop 5
           End If    

           zz = zp(j) - zp(i)                    !distance between grid points j and i
           dz = zp(j+1) - (zp(i) - d2)           !the part of the intervall between zp(j) and zp(j+1) which is already within i-d/2
           !if(dz < epsilon(dz)) dz = epsilon(dz) !bei unguenstiger Kombination von sig und ngrid kann dz unter Machinengenauigkeit epsilon liegen, dann ist x(2) = x(1) + dz = x(1) -> das fuehrt zu Abbruch in Spline Interpolation           
           x_int(n)    = 0.0              !array containing x-values for spline integration 
           n2_int(n)   = rhopjk + (rhopjp1k-rhopjk) / (dzp) * (dzp-dz)                 !integrand für n2 (=rhop) linear interpoliert für den Punk zp(i)-d/2             
           n3_int(n)   = 0.0                                                              !integrand für n3: rhop*(d2**2 - z'**2), da hier gerade z' = d2 -> integrand wird hier = 0!! 
           nv2_int(n) = rhopjk*zz +(rhopjp1k*(zp(j+1)-zp(i))-rhopjk*zz)/(dzp)*(dzp-dz) !analog
          
          
       Else If (zp(j) > (zp(i)-d2) .and. zp(j) <= (zp(i)+d2)) Then !grid point j  within i+-d2 
           
           n = n + 1
           x_int(n) = x_int(n-1) +  dz !first time in this If condition, dz is stil the old value from above!
           zz = zp(j) - zp(i)
           dz = dzp
           n2_int(n)  = rhopjk
           n3_int(n)  = rhopjk * (d2**2 - zz**2)
           nv2_int(n) = rhopjk * zz 
           
           If (zp(j) < (zp(i)+d2) .and. zp(j+1) >= (zp(i)+d2) ) Then !zp(j) is still within zp(i)+d2 but zp(j+1) is already outside zp(i)+d2

           dz = zp(i) + d2 - zp(j)
           !If(dz <= epsilon(dz)) exit   !wie oben, kann auch hier bei ungluecklicher Wahl von sig und ngrid dz < eps werden und somit x(n) = x(n-1) -> Abbruch in Spline interpolation. Dann einfach ngrid aendern!
           zz = zp(j) - zp(i)
           n = n + 1
           x_int(n)   = x_int(n-1) + dz
!            If(x_int(n) == x_int(n-1)) Then
!              n = n - 1
!              exit
!            End If  
           n2_int(n)  = rhopjk + (rhopjp1k-rhopjk) / (dzp) * dz
           n3_int(n)  = 0.0 !Begründung wie oben
           nv2_int(n) = rhopjk*zz + (rhopjp1k*(zp(j+1)-zp(i)) - rhopjk*zz) / (dzp) * dz 
    
           End If
           
       End If


    End Do

    !spline integration
    xlo = x_int(1)
    xhi = x_int(n)

    If(n > NMAX) Then
        write(*,*) 'Increase NMAX in FMT_Weighted_Densities (also in AD Routine!!)'
        stop 5 
    End If    

    call spline ( x_int, n2_int, n, 1.E30, 1.E30, y2_n2 )
    call spline ( x_int, n3_int, n, 1.E30, 1.E30, y2_n3 )
    call spline ( x_int, nv2_int, n, 1.E30, 1.E30, y2_nv2 )
    
    call splint_integral ( x_int, n2_int, y2_n2, n, xlo, xhi, int_n2 )
    call splint_integral ( x_int, n3_int, y2_n3, n, xlo, xhi, int_n3 )
    call splint_integral ( x_int, nv2_int, y2_nv2, n, xlo, xhi, int_nv2 )

    !weighted densities
    n2(i)  = n2(i) + PI * dhs(k) * int_n2
    n1(i)  = n1(i) + 0.5 * int_n2
    n0(i)  = n0(i) + int_n2/dhs(k)
    n3(i)  = n3(i) + PI * int_n3
    nv2(i) = nv2(i) -2.0 * PI * int_nv2    
    nv1(i) = nv1(i) -int_nv2 / dhs(k)
   
!     If(i < 5) Then 
!     write(*,*)'i dens',i,nv2(i),n1(i),n0(i)
!     write(*,*)'i dens',i,n3(i),nv1(i),nv2(i)
!     end if
!    
!    If(i > 95) then 
!       write(*,*)'i dens',i,nv2(i),n1(i),n0(i)
!       write(*,*)'i dens',i,n3(i),nv1(i),nv2(i)
!    End If 
!    
    End Do
   
  ! pause
   
End Do    


!derivatives of FMT helmholtz energy density w.r.t. weighted densities 
Do i=user%xs-maxval((fa(1:ncomp)+1)/2),user%xe+maxval((fa(1:ncomp)+1)/2)     
    !weils kürzer ist
    nn0 = n0(i)
    nn1 = n1(i)
    nn2 = n2(i)
    nn3 = n3(i)
    nnv1 = nv1(i)
    nnv2 = nv2(i) 
    
    zms = 1.0 - nn3
    zms2 = zms*zms
    zms3 = zms2*zms
    logzms = log(zms)
    if(isnan(logzms)) Then
        write(*,*)'Surface Tension Code: zms < 0, log(zms) undefined FMT_Weighted_Densities'
        stop 5     
    End If
    phi_dn0(i) = -logzms
    phi_dn1(i) = nn2/zms 
    phi_dn2(i) = nn1/zms + 3.0*(nn2*nn2-nnv2*nnv2) * (nn3+zms2*logzms) / (36.0*PI*nn3*nn3*zms2)  
    
    phi_dn3(i) = nn0/zms + (nn1*nn2-nnv1*nnv2)/zms2 - (nn2**3-3.0*nn2*nnv2*nnv2) * (nn3*(nn3**2-5.0*nn3+2.0)+2.0*zms3*logzms) &
                       / (36.0*PI*nn3**3*zms3) 

    phi_dnv1(i) = -nnv2/zms

    phi_dnv2(i) = -nnv1/zms - 6.0*nn2*nnv2*(nn3+zms2*logzms)/(36.0*PI*nn3**2*zms2)  

End Do

End Subroutine FMT_Weighted_Densities







Subroutine FMT_dFdrho(i,fa,user,phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2,dF_drho_FMT)

Use BASIC_VARIABLES, Only: ncomp,parame
Use EOS_VARIABLES, Only: dhs
Use mod_DFT, Only: zp,dzp
Use PARAMETERS, ONLY: PI

!PETSc module
Use PetscManagement

!passed
INTEGER, INTENT(IN) :: i         !the grid point at which to calculate dFdrho
INTEGER, INTENT(IN) :: fa(ncomp)
Type (userctx)                                     :: user
REAL,dimension(user%gxs:user%gxe),INTENT(IN)       :: phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2
REAL,dimension(user%xs:user%xe,ncomp), INTENT(OUT)       :: dF_drho_FMT


!local
INTEGER :: j,k,n
INTEGER :: fa2
REAL    :: dz,d2,zz, zz_jp1
INTEGER, parameter :: NMAX = 800
REAL,dimension(NMAX) :: x_int,y_int,y0_int,y1_int,y2_int,y3_int,yv1_int,yv2_int,y2
REAL      :: xhi,xlo,integral,int0,int1,int2,int3,intv1,intv2
REAL      :: at_j, at_jp1


!Das Integral (Gleichung A1 in Gross DFT 2009) wird hier in einem Schlag berechnet!
!Falls die einzelnen Terme einzeln integriert werden sollen, einfach die auskommentierte Version verwenden

Do k=1,ncomp !das einzige, das hier von k abhaengt, sind fa und dhs!! die Ableitungen phi_dn... sind nicht Komponentenspez, da die 
             !gewichteten Dichten ja uch nicht mehr Komponentenspez sind (n_i = sum n_i(k))

  n = 1
  fa2 = ( fa(k) + 5 ) / 2 !number of grid points within dhs/2
  d2  = dhs(k)/2.0

  Do j = i-fa2,i+fa2  
   
       If( ( zp(i)-zp(j+1) ) < d2  .and. ( zp(i) - zp(j) ) >= d2   ) Then !the position of j+1 is already within i-d/2 while j is still outside this range in this case, the integration steplength (dz) is just the distance, which j+1 overlaps with i-d/2 and what is integrated is the interpolated value of the integrand
          x_int(n) = 0.0
          If(n/=1) Then
             write(*,*) 'Surface Tension Code: error in FMT_dFdrho, n should be 1 here!'
             stop 5
          End If
          dz = zp(j+1) - (zp(i) - d2)
          zz = zp(j) - zp(i)
          zz_jp1 = zp(j+1) - zp(i)     
          at_j = phi_dn0(j)/dhs(k) + 0.5*phi_dn1(j) + PI*dhs(k)*phi_dn2(j) &
                        + PI*phi_dn3(j)*(d2**2 - zz**2) + phi_dnv1(j)*zz/dhs(k) + 2.0*PI*phi_dnv2(j)*zz

          at_jp1 = phi_dn0(j+1)/dhs(k) + 0.5*phi_dn1(j+1) + PI*dhs(k)*phi_dn2(j+1) &
                        + PI*phi_dn3(j+1)*(d2**2 - zz_jp1**2) + phi_dnv1(j+1)*zz_jp1/dhs(k) + 2.0*PI*phi_dnv2(j+1)*zz_jp1

          y_int(n) = at_j + (at_jp1-at_j)/dzp * (dzp-dz)   !lineare interpolation genau, wie in FMT_Weighted_Densities         
!           y0_int(n) = phi_dn0(j) + (phi_dn0(j+1) -phi_dn0(j))/dzp * (dzp-dz)
!           y1_int(n) = phi_dn1(j) + (phi_dn1(j+1) -phi_dn1(j))/dzp * (dzp-dz)          
!           y2_int(n) = phi_dn2(j) + (phi_dn2(j+1) -phi_dn2(j))/dzp * (dzp-dz)          
!           y3_int(n) = phi_dn3(j)*(d2**2-zz**2) + (phi_dn3(j+1)*(d2**2-zz_jp1**2) - phi_dn3(j)*(d2**2-zz**2) )/dzp * (dzp-dz)
!           yv1_int(n) = phi_dnv1(j)*zz + (phi_dnv1(j+1)*zz_jp1 - phi_dnv1(j)*zz)/dzp * (dzp-dz)          
!           yv2_int(n) = phi_dnv2(j)*zz + (phi_dnv2(j+1)*zz_jp1 - phi_dnv2(j)*zz)/dzp * (dzp-dz)
!           
          
       Else If (zp(j) > (zp(i)-d2) .and. zp(j) <= (zp(i)+d2)) Then !grid points j and j+1 are completely within i+-d2 
          
          n = n + 1
          zz = zp(j) - zp(i)
          x_int(n) = x_int(n-1) + dz 
          
          
          y_int(n) = phi_dn0(j)/dhs(k) + 0.5*phi_dn1(j) + PI*dhs(k)*phi_dn2(j) &
                        + PI*phi_dn3(j)*(d2**2 - zz**2) + phi_dnv1(j)*zz/dhs(k) + 2.0*PI*phi_dnv2(j)*zz            
          
          
!           y0_int(n) = phi_dn0(j)
!           y1_int(n) = phi_dn1(j)
!           y2_int(n) = phi_dn2(j)
!           y3_int(n) = phi_dn3(j)*(d2**2-zz**2)
!           yv1_int(n) = phi_dnv1(j)*zz
!           yv2_int(n) = phi_dnv2(j)*zz
!           
          
          dz = dzp          
          
          If (zp(j) < (zp(i)+d2) .and. zp(j+1) > (zp(i)+d2) ) Then !zp(j) is still within zp(i)+d2 but zp(j+1) is already out side zp(i)+d2

          n = n + 1
          zz = zp(j) - zp(i)
          zz_jp1 = zp(j+1) - zp(i)
          dz = zp(i) + d2 - zp(j)
          x_int(n) = x_int(n-1) + dz
          at_j = phi_dn0(j)/dhs(k) + 0.5*phi_dn1(j) + PI*dhs(k)*phi_dn2(j) &
                        + PI*phi_dn3(j)*(d2**2 - zz**2) + phi_dnv1(j)*zz/dhs(k) + 2.0*PI*phi_dnv2(j)*zz
          at_jp1 = phi_dn0(j+1)/dhs(k) + 0.5*phi_dn1(j+1) + PI*dhs(k)*phi_dn2(j+1) &
                        + PI*phi_dn3(j+1)*(d2**2 - zz_jp1**2) + phi_dnv1(j+1)*zz_jp1/dhs(k) + 2.0*PI*phi_dnv2(j+1)*zz_jp1
          y_int(n) = at_j + (at_jp1-at_j)/dzp * dz

          
!           y0_int(n) = phi_dn0(j) + (phi_dn0(j+1) - phi_dn0(j))/dzp * dz
!           y1_int(n) = phi_dn1(j) + (phi_dn1(j+1) - phi_dn1(j))/dzp * dz
!           y2_int(n) = phi_dn2(j) + (phi_dn2(j+1) - phi_dn2(j))/dzp * dz
!           y3_int(n) = phi_dn3(j)*(d2**2 - zz**2) + (phi_dn3(j+1)*(d2**2 - zz_jp1**2) - phi_dn3(j)*(d2**2 - zz**2))/dzp * dz
!           yv1_int(n) = phi_dnv1(j)*zz + (phi_dnv1(j+1)*zz_jp1 - phi_dnv1(j)*zz)/dzp * dz  
!           yv2_int(n) = phi_dnv2(j)*zz + (phi_dnv2(j+1)*zz_jp1 - phi_dnv2(j)*zz)/dzp * dz  
!           
          
                  
          End If
          
       End If


  End Do

  !spline integration
  xlo = x_int(1)
  xhi = x_int(n)

  If(n > NMAX) Then
      write(*,*)'Surface Tension code: Increase NMAX in FMT_dFdrho (also in AD Routine!!)'
      stop 5 
  End If
  
  call spline(x_int,y_int,n,1.E30,1.E30,y2)
  call splint_integral(x_int,y_int,y2,n,xlo,xhi,integral)
  
!   call spline ( x_int, y0_int, n, 1.E30, 1.E30, y2 )
!   call splint_integral ( x_int, y0_int, y2, n, xlo, xhi, int0 )
!   call spline ( x_int, y1_int, n, 1.E30, 1.E30, y2 )
!   call splint_integral ( x_int, y1_int, y2, n, xlo, xhi, int1 )
!   call spline ( x_int, y2_int, n, 1.E30, 1.E30, y2 )
!   call splint_integral ( x_int, y2_int, y2, n, xlo, xhi, int2 )
!   call spline ( x_int, y3_int, n, 1.E30, 1.E30, y2 )
!   call splint_integral ( x_int, y3_int, y2, n, xlo, xhi, int3 )
!   call spline ( x_int, yv1_int, n, 1.E30, 1.E30, y2 )
!   call splint_integral ( x_int, yv1_int, y2, n, xlo, xhi, intv1 )
!   call spline ( x_int, yv2_int, n, 1.E30, 1.E30, y2 )
!   call splint_integral ( x_int, yv2_int, y2, n, xlo, xhi, intv2 )
!   
  
  !dF_drho_FMT(i,k) = int0/dhs(k) + 0.5*int1 + PI*dhs(k)*int2 + PI*int3 + intv1/dhs(k) + 2.0*PI*intv2
  dF_drho_FMT(i,k) = integral*parame(k,1)
  
End Do


End Subroutine FMT_dFdrho


End Module mod_DFT_FMT 
