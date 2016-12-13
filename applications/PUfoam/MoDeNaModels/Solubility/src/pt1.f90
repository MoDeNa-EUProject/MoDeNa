!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE f_pt1 ( fres )

  USE EOS_VARIABLES, ONLY: PI, t, parame, mseg, rho
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      :: fres
  !-----------------------------------------------------------------------------
  REAL                                   :: I1_dft, I2_dft, order1
  !-----------------------------------------------------------------------------

  order1 = mseg(1)*mseg(1)*parame(1,2)**3 *parame(1,3) / t

  CALL f_dft (I1_dft,I2_dft)

  fres  = 2.0 * PI * rho * I1_dft * order1

END SUBROUTINE f_pt1



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE phipt1 ( my_pt1, p_pt1 )

  USE EOS_VARIABLES, ONLY: KBOL, t, rho
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      :: my_pt1
  REAL, INTENT(IN OUT)                   :: p_pt1

  !-----------------------------------------------------------------------------
  REAL                                   :: fres, pgesdz, pgesd2, pgesd3
  REAL                                   :: zres
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! density iteration
  !-----------------------------------------------------------------------------
  CALL pressure_pt1 ( p_pt1, pgesdz, pgesd2, pgesd3 )

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------
  zres = (p_pt1 * 1.E-30) / (kbol*t*rho)
  !      zres = zges - 1.0

  !-----------------------------------------------------------------------------
  ! f res
  !-----------------------------------------------------------------------------
  CALL f_pt1 ( fres )

  my_pt1 = fres + zres

END SUBROUTINE phipt1


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE pressure_pt1 ( pges, pgesdz, pgesd2, pgesd3 )

  USE EOS_VARIABLES, ONLY: KBOL, t, rho, eta, tfr
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      :: pges
  REAL, INTENT(OUT)                      :: pgesdz
  REAL, INTENT(OUT)                      :: pgesd2
  REAL, INTENT(OUT)                      :: pgesd3

  !-----------------------------------------------------------------------------
  REAL                                   :: dzetdv, dicht, dist, fact, z3t
  REAL                                   :: fres1, fres2, fres3, fres4, fres5, fres
  REAL                                   :: df_dr, df_drdr
  REAL                                   :: tfr_1, tfr_2, tfr_3, tfr_4, tfr_5
  !-----------------------------------------------------------------------------

  dzetdv = eta*rho
  z3t = eta/rho

  IF (eta > 1D-1) THEN
     fact = 1.0
  ELSE IF (eta <= 0.1.AND.eta > 0.01) THEN
     fact = 10.0
  ELSE
     fact = 100.0
  END IF
  dist = eta*3.E-3 *fact

  dicht  = eta
  eta  = dicht - 2.0*dist
  rho = eta/z3t
  CALL f_pt1 ( fres )
  fres1  = fres
  tfr_1  = tfr
  eta  = dicht - dist
  rho = eta/z3t
  CALL f_pt1 ( fres )
  fres2  = fres
  tfr_2  = tfr
  eta  = dicht + dist
  rho = eta/z3t
  CALL f_pt1 ( fres )
  fres3  = fres
  tfr_3  = tfr
  eta  = dicht + 2.0*dist
  rho = eta/z3t
  CALL f_pt1 ( fres )
  fres4  = fres
  tfr_4  = tfr
  eta  = dicht
  rho = eta/z3t
  CALL f_pt1 ( fres )
  fres5  = fres
  tfr_5  = tfr

  !-----------------------------------------------------------------------------
  !      ptfr   = (-tfr_4+8.0*tfr_3-8.0*tfr_2+tfr_1)/(12.0*dist)
  !     &           *dzetdv*(KBOL*T)/1.E-30
  !      ztfr =ptfr /( rho * (KBOL*t) / 1.E-30)
  !      ptfrdz = (-tfr_4+16.0*tfr_3-3.E1*tfr_5+16.0*tfr_2-tfr_1)
  !     &             /(12.0*(dist**2 ))* dzetdv*(KBOL*T)/1.E-30
  !     &         + (-tfr_4+8.0*tfr_3-8.0*tfr_2+tfr_1)
  !     &            /(12.0*dist) * 2.0 *eta*6.0/PI/D
  !     &                               *(KBOL*T)/1.E-30
  !      ztfrdz=ptfrdz/( rho*(kbol*T)/1.E-30 ) -  ztfr/eta
  !      write (*,*) eta,ztfr,ztfrdz

  !      dtfr_dr   = (-tfr_4+8.0*tfr_3-8.0*tfr_2+tfr_1)/(12.0*dist)
  !      write (*,*) eta,dtfr_dr
  !      stop
  !-----------------------------------------------------------------------------

  df_dr   = (-fres4+8.0*fres3-8.0*fres2+fres1)/(12.0*dist)
  df_drdr = (-fres4+16.0*fres3-3.E1*fres5+16.0*fres2-fres1)  &
       /(12.0*(dist**2 ))

  pges   = (-fres4+8.0*fres3-8.0*fres2+fres1)  &
       /(12.0*dist) *dzetdv*(kbol*t)/1.E-30

  pgesdz = (-fres4+16.0*fres3-3.E1*fres5+16.0*fres2-fres1)  &
       /(12.0*(dist**2 ))* dzetdv*(kbol*t)/1.E-30  &
       + (-fres4+8.0*fres3-8.0*fres2+fres1) /(12.0*dist) * 2.0 *rho  &
       *(kbol*t)/1.E-30

  pgesd2 = (fres4-2.0*fres3+2.0*fres2-fres1) /(2.0*dist**3 )  &
       * dzetdv*(kbol*t)/1.E-30  &
       + (-fres4+16.0*fres3-3.E1*fres5+16.0*fres2-fres1) /(12.0*(dist**2 ))  &
       * 4.0 *rho *(kbol*t)/1.E-30 + (-fres4+8.0*fres3-8.0*fres2+fres1)  &
       /(12.0*dist) * 2.0 /z3t *(kbol*t)/1.E-30
  pgesd3 = (fres4-4.0*fres3+6.0*fres5-4.0*fres2+fres1) /(dist**4)  &
       * dzetdv*(kbol*t)/1.E-30 + (fres4-2.0*fres3+2.0*fres2-fres1)  &
       /(2.0*dist**3 ) * 6.0 *rho *(kbol*t)/1.E-30  &
       + (-fres4+16.0*fres3-3.E1*fres5+16.0*fres2-fres1)  &
       /(12.0*dist**2 )* 6.0 /z3t *(kbol*t)/1.E-30


END SUBROUTINE pressure_pt1
