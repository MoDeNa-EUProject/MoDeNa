!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! This module ........... 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
Module rdf_variables

  implicit none
  save

  real, dimension(0:20)                         :: fac(0:20)

End Module rdf_variables




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE rdf_int
!
! this subroutine calculates the radial pair correlation function of
! hard spheres (mseg=1) or tangent-sphere hard chains (with "mseg"
! spherical segments). The pair correlation functions are for the
! Percus-Yevick closure and the analytic solution was proposed by
! Y. Tang, B.C.-Y. Lu, J. Chem. Phys. 105.1996.8262.
! The density "dense" is the packing fraction, the radius "radius" is
! dimensionless in terms of the hard-sphere diameter. The output is
! the segment-segment radial distribution function g(r).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE rdf_int ( dense, mseg, radius, rdf )

  USE rdf_variables
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: dense
  REAL, INTENT(IN)                       :: mseg
  REAL, INTENT(IN)                       :: radius
  REAL, INTENT(OUT)                      :: rdf

  !-----------------------------------------------------------------------------
  INTEGER                                :: n
  REAL                                   :: rg, omega, radind, cfun, facul
  !-----------------------------------------------------------------------------
  ! OPEN (10,file='rdf.xlo')

  ! determine the faculty of integer numbers from 0 to 20. If this
  ! subroutine is called repeatedly, it makes sense to do this only once.

  DO n = 0, 20
     fac(n) = facul(n)
  END DO

  omega = ( mseg - 1.0 ) / mseg
  rg = 0.0

  DO  n = 0,4

     radind = radius - 1.0 - REAL(n)
     rg = rg + (-12.0*dense)**n *( (1.0+dense/2.0-omega*(1.0-dense))  &
          *cfun(dense,omega,2,n,(n+1),radind)  &
          +(1.0+dense*2.0-(1.0+dense/2.0)*omega +(1.0-dense)/2.0*omega**2)  &
          *cfun(dense,omega,1,n,(n+1),radind)  &
          +((1.0+2.0*dense)*omega-(1.0-dense)*omega**2 )  &
          *cfun(dense,omega,0,n,(n+1),radind)  &
          -((1.0+dense/2.0)*omega-(1.0-dense)/2.0*omega**2 )  &
          *cfun(dense,omega,1,n,(n+1),(radind-1.0))  &
          -((1.0+2.0*dense)*omega-(1.0-dense)*omega**2 )  &
          *cfun(dense,omega,0,n,(n+1),(radind-1.0))   )

  END DO

  rdf = rg / radius
  ! write (*,*) 'r, g(r)',radius,rg/radius

END SUBROUTINE rdf_int


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION lnfac(xx)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: xx
  !-----------------------------------------------------------------------------
  INTEGER                                :: j
  REAL, SAVE                             :: cof(6)
  REAL, SAVE                             :: stp
  REAL                                   :: ser,tmp,x,y
  !-----------------------------------------------------------------------------
  DATA cof,stp/76.18009172947146,-86.50532032941677,  &
       24.01409824083091,-1.231739572450155,.1208650973866179E-2,  &
       -.5395239384953E-5,2.5066282746310005/

  x = xx
  y = x

  tmp = x + 5.5
  tmp = ( x + 0.5 ) * LOG( tmp ) - tmp
  ser = 1.000000000190015
  DO  j = 1, 6
     y = y + 1.0
     ser = ser + cof(j) / y
  END DO
  lnfac = tmp + LOG( stp * ser / x )

END FUNCTION lnfac



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION facul(n)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: n

  !-----------------------------------------------------------------------------
  INTEGER                                :: j
  INTEGER, SAVE                          :: ntop
  REAL, SAVE                             :: table(165)
  REAL                                   :: lnfac,mult
  !-----------------------------------------------------------------------------
  DATA ntop,table(1)    / 0, 1.0 /

  IF ( n < 0 ) THEN
     write (*,*) ' Negative factorial in FUNCTION fac'
     stop
  ELSE IF ( n <= ntop ) THEN
     facul = table(n+1)
  ELSE IF ( n <= 165 ) THEN
     DO  j = ntop+1, n
        mult = REAL(j)
        table(j+1)=mult*table(j)
     END DO
     ntop = n
     facul = table(n+1)
  ELSE
     facul = LOG( lnfac( REAL(n)+1.0 ) )
  END IF

END FUNCTION facul



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

COMPLEX FUNCTION anroot(dens,omeg,m)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: dens
  REAL, INTENT(IN OUT)                   :: omeg
  INTEGER, INTENT(IN)                    :: m

  !-----------------------------------------------------------------------------
  REAL                                   :: coeff(3),qq,rr,aa,bb,ohm, pi
  REAL, SAVE                             :: roh2
  COMPLEX, SAVE                          :: anr0,anr1,anr2
  !-----------------------------------------------------------------------------
  DATA roh2,anr0,anr1,anr2    / 0.0, (0.0, 0.0), (0.0, 0.0), (0.0, 0.0) /


  pi   = 3.14159265359

  IF (dens /= roh2) THEN
     roh2 = dens
     coeff(1)=(6.0*dens-(1.0-dens)*omeg) /(1.0-dens)
     coeff(2)=(18.0*dens**2 - 6.0*dens*(1.0-dens)*omeg) /((1.0-dens)**2)
     coeff(3)=(-12.0*dens*(1.0+2.0*dens) +12.0*dens*(1.0-dens)*omeg) /((1.0-dens)**2)
     qq = ( coeff(1)**2 - 3.0*coeff(2) ) / 9.0
     rr = ( 2.0*coeff(1)**3 - 9.0*coeff(1)*coeff(2) + 27.0*coeff(3) ) / 54.0

     IF ( rr**2 < qq**3 ) THEN

        ohm = ACOS( rr/qq**(3.0/2.0) )
        anr0=-2.0* SQRT(qq)* COS( ohm/3.0        )-coeff(1)/3.0
        anr1=-2.0* SQRT(qq)* COS((ohm+2.0*pi)/3.0)-coeff(1)/3.0
        anr2=-2.0* SQRT(qq)* COS((ohm-2.0*pi)/3.0)-coeff(1)/3.0

     ELSE

        aa = - ( ABS(rr) + SQRT( rr**2 -qq**3 ) )**(1.0/3.0)
        aa = rr/ ABS(rr) * aa
        IF ( aa /= 0.0 ) THEN
           bb = qq/aa
        ELSE
           bb = 0.0
        END IF

        anr0 = CMPLX(   aa+bb      -coeff(1)/3.0 ,  0.0 )
        anr1 = CMPLX((-(aa+bb)/2.0 -coeff(1)/3.0),  (SQRT(3.0)/2.0*(aa-bb)))
        anr2 = CMPLX((-(aa+bb)/2.0 -coeff(1)/3.0), -(SQRT(3.0)/2.0*(aa-bb)))

     END IF

  END IF

  IF (m == 0) THEN
     anroot = anr0
  ELSE IF (m == 1) THEN
     anroot = anr1
  ELSE IF (m == 2) THEN
     anroot = anr2
  ELSE
     WRITE (*,*) ' Only 3 roots do exist. Asked for: ',m+1
     STOP
  END IF



  !      write (*,*) '1',anroot,(anroot**3
  !     &          + coeff(1)*anroot**2 + coeff(2)*anroot + coeff(3))
  !      DO 2 i = 1,3
  !      polyn =anroot**3 +coeff(1)*anroot**2 +coeff(2)*anroot+coeff(3)
  !      deriv =3.0*anroot**2 +2.0*coeff(1)*anroot+coeff(2)
  !      realrt = anroot
  !      imagrt = AIMAG(anroot)
  !      realpt = polyn/deriv
  !      imagpt = AIMAG(polyn/deriv)
  !      realrt = realrt - realpt
  !      imagrt = imagrt - imagpt
  !c      anroot = anroot - polyn/deriv
  !      anroot = CMPLX(realrt,imagrt)
  !      dummy = anroot
  !      write (*,*) '2',polyn/deriv,anroot
  ! 2    CONTINUE
  !      IF (m.eq.2) stop

END FUNCTION anroot



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

COMPLEX FUNCTION afun (dense,om,an1,an2,ak,aalpha)

  USE rdf_variables
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i,j,an1,an2,ak
  REAL                                   :: dense,om
  COMPLEX                                :: aalpha, af
  !-----------------------------------------------------------------------------

  af = 0.0

  DO  i = 0, an2
     DO  j = 0, an2
        IF ((i+j) <= an2 .AND. (i+j+j) >= (ak-an1)) THEN

           af = af + fac(an2)*fac(i+an1+j+j)  &
                /(fac(i)*fac(j)*fac(an2-i-j)*fac(i+j+j+an1-ak))  &
                * (1.0+dense/2.0-(1.0-dense)/2.0*om)**REAL(i)  &
                * (1.0+2.0*dense-(1.0-dense)*om)**REAL(an2-i-j)  &
                * ((1.0-dense)**2 *om/(12.0*dense))**REAL(j)  &
                * aalpha**REAL(an1+i+j+j-ak)

        END IF
     END DO
  END DO

  afun = af

END FUNCTION afun



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

COMPLEX FUNCTION bfun (dense,om,bn1,bn2,bn3,bi,bal)

  USE rdf_variables
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: k,l,bn1,bn2,bn3,bi,bal,beta,gamma
  REAL                                   :: dense, om, imsum,factor
  COMPLEX                                :: talpha, anroot, afun, bf,comfra
  !-----------------------------------------------------------------------------

  IF (bal == 0) THEN
     beta = 1
     gamma= 2
  ELSE IF (bal == 1) THEN
     beta = 0
     gamma= 2
  ELSE IF (bal == 2) THEN
     beta = 0
     gamma= 1
  ELSE
     WRITE(*,*) ' no proper value for bal in BFUN !!'
     STOP
  END IF

  !      write(*,*) ' BFUN ',dense,om,bn1,bn2,bn3,bi,bal
  bf = (0.0, 0.0)
  imsum = 0.0

  DO  k = 0, (bn3-bi)
     DO  l = 0, (bn3-bi-k)
        talpha = anroot(dense,om,bal)
        !          bf = bf
        !     &  +(-1.0)**REAL(bn3-bi-k)*fac(bn3-1+l)*fac(bn3+bn3-1-bi-k-l)
        !     &      /( fac(k)*fac(l)*fac(bn3-bi-k-l)*fac(bn3-1)*fac(bn3-1) )
        !     &      *AFUN(dense,om,bn1,bn2,k,talpha)
        !     &      /( (talpha
        !     &           - anroot(dense,om,beta))**REAL(bn3+l)
        !     &       *(talpha
        !     &           -  anroot(dense,om,gamma))**REAL(bn3+bn3-bi-k-l) )

        factor=(-1.0)**REAL(bn3-bi-k)*fac(bn3-1+l)*fac(bn3+bn3-1-bi-k-l)  &
             /( fac(k)*fac(l)*fac(bn3-bi-k-l)*fac(bn3-1)*fac(bn3-1) )
        comfra=   afun(dense,om,bn1,bn2,k,talpha) /( (talpha  &
             - anroot(dense,om,beta))**REAL(bn3+l) *(talpha  &
             -  anroot(dense,om,gamma))**REAL(bn3+bn3-bi-k-l) )
        imsum = imsum + factor*AIMAG(comfra)
        bf = bf + factor * comfra

     END DO
  END DO

  bfun = 1.0 / (1.0-dense)**REAL(bn3+bn3) * bf

END FUNCTION bfun



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION cfun (dense,om,cn1,cn2,cn3,rad)

  USE rdf_variables
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: cn1,cn2,cn3, i, alph
  REAL                                   :: dense,om,rad,fact
  COMPLEX                                :: anroot,bfun,cf,root,swing,compl
  !-----------------------------------------------------------------------------

  cf = CMPLX( 0.0, 0.0 )
  IF (ABS(rad) < 1.e-10) rad = 0.0

  DO  alph = 0,2
     DO  i = 1,cn3
        IF (rad >= 0.) THEN
           root = anroot(dense,om,alph)
           swing = root*rad
           compl =bfun(dense,om,cn1,cn2,cn3,i,alph)*CEXP(swing)
           fact   = rad**REAL(i-1) / fac(i-1)
           cf = cf + compl * fact
        END IF
     END DO
  END DO
  cfun = cf

  IF ( (cn1+cn2) < (cn3+cn3+cn3) ) THEN
     cfun = cf
  ELSE IF ( (cn1+cn2) == (cn3+cn3+cn3) ) THEN
     IF (rad == 0.0) THEN
        cfun = cf + (1.0+dense/2.0)**REAL(cn2) /((1.0-dense)**REAL(cn3+cn3))
        WRITE (*,*) 'in Schleife',cfun
     ELSE
        cfun = cf
     END IF
  ELSE
     WRITE(*,*) ' no proper case in CFUN !!'
     STOP
  END IF

END FUNCTION cfun
