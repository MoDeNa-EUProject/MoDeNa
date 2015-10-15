Module mod_ChemPot



Implicit None



private

public :: Chemical_Potential
public :: PCSAFT_const
REAL, public, allocatable :: ChemPot_tot(:)
REAL, public, allocatable :: ChemPot_res(:)

 contains




Subroutine Chemical_Potential

Use BASIC_VARIABLES
Use EOS_VARIABLES
Use PARAMETERS
USe EOS_CONSTANTS

Integer :: k,i,j
REAL :: z0t,z1t,z2t,z3t
REAL :: z0,z1,z2,z3,zms
REAL :: z0_rk,z1_rk,z2_rk,z3_rk
REAL, allocatable :: mhs(:), mhc(:)
REAL, allocatable :: gij(:,:), gij_rk(:,:), dij_ab(:,:)

!DISP HERE!!!
INTEGER :: m
REAL :: m_mean,I1,I2,I1_rk,I2_rk
REAL :: ord1_rk,ord2_rk
REAL :: c1_con,c2_con,c1_rk
REAL :: order1,order2
REAL :: apar(0:6),bpar(0:6)
REAL, allocatable :: m_rk(:),mdsp(:)
REAL, allocatable :: ap_rk(:,:),bp_rk(:,:)

Allocate(mdsp(ncomp),m_rk(ncomp),ap_rk(ncomp,0:6),bp_rk(ncomp,0:6))
Allocate(sig_ij(ncomp,ncomp),uij(ncomp,ncomp))

kij = 0.0



Allocate(mhs(ncomp), mhc(ncomp), ChemPot_res(ncomp), ChemPot_tot(ncomp)  )
Allocate(gij(ncomp,ncomp), gij_rk(ncomp,ncomp), dij_ab(ncomp,ncomp) )


!belege dichteunabhängige Parameter (z0z,z1t,z2t,z3t)
z0t = PI / 6.0 * SUM( xx(1:ncomp) * mseg(1:ncomp) )
z1t = PI / 6.0 * SUM( xx(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp) )
z2t = PI / 6.0 * SUM( xx(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**2 )
z3t = PI / 6.0 * SUM( xx(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )



! --- Eq.(A.8) ---------------------------------------------------------
rho = eta / z3t   !total number density [particles/A^3]
z0  = z0t * rho
z1  = z1t * rho
z2  = z2t * rho
z3  = z3t * rho

zms    = 1.0 - eta

call PCSAFT_const(ap,bp) !get PCSAFT constants

DO  k = 1, ncomp
 
  z0_rk = PI/6.0 * mseg(k)
  z1_rk = PI/6.0 * mseg(k) * dhs(k)
  z2_rk = PI/6.0 * mseg(k) * dhs(k)*dhs(k)
  z3_rk = PI/6.0 * mseg(k) * dhs(k)**3
  
  ! --------------------------------------------------------------------
  ! d(f)/d(x) : hard sphere contribution
  ! --------------------------------------------------------------------
  mhs(k) =  6.0/PI* (  3.0*(z1_rk*z2+z1*z2_rk)/zms + 3.0*z1*z2*z3_rk/zms/zms  &
                + 3.0*z2*z2*z2_rk/z3/zms/zms + z2**3 *z3_rk*(3.0*z3-1.0)/z3/z3/zms**3   &
                + ((3.0*z2*z2*z2_rk*z3-2.0*z2**3 *z3_rk)/z3**3 -z0_rk) *LOG(zms)  &
                + (z0-z2**3 /z3/z3)*z3_rk/zms  )

                
  ! --------------------------------------------------------------------
  ! d(f)/d(x) : chain term
  ! --------------------------------------------------------------------
  DO i = 1, ncomp
    DO j = 1, ncomp
      dij_ab(i,j) = dhs(i)*dhs(j) / ( dhs(i) + dhs(j) )
    END DO
  END DO
  
  DO i = 1, ncomp
    DO j = 1, ncomp
      gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 /zms**3 
      gij_rk(i,j) = z3_rk/zms/zms  &
                      + 3.0*dij_ab(i,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
                      + dij_ab(i,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
    END DO
  END DO
  
  mhc(k) = 0.0
  DO i = 1, ncomp
    mhc(k) = mhc(k) + xx(i)*rho * (1.0-mseg(i)) / gij(i,i) * gij_rk(i,i)
  END DO
  mhc(k) = mhc(k) + ( 1.0-mseg(k)) * LOG( gij(k,k) )

  ! --------------------------------------------------------------------
  ! PC-SAFT:  d(f)/d(x) : dispersion contribution
  ! --------------------------------------------------------------------
    
    m_mean = z0t / (PI/6.0)
    m_rk(k) = ( mseg(k) - m_mean ) / rho
    
    
    DO m = 0, 6
         apar(m) = ap(m,1) + (1.0-1.0/m_mean)*ap(m,2)  &
            + (1.0-1.0/m_mean)*(1.0-2.0/m_mean)*ap(m,3)
         bpar(m) = bp(m,1) + (1.0-1.0/m_mean)*bp(m,2)  &
            + (1.0-1.0/m_mean)*(1.0-2.0/m_mean)*bp(m,3)

        ! --- derivatives of apar, bpar to rho_k ---------------------------
        ap_rk(k,m) = m_rk(k)/m_mean**2 * ( ap(m,2) + (3.0 -4.0/m_mean) *ap(m,3) )
        bp_rk(k,m) = m_rk(k)/m_mean**2 * ( bp(m,2) + (3.0 -4.0/m_mean) *bp(m,3) )
    END DO

    I1    = 0.0
    I2    = 0.0
    I1_rk = 0.0
    I2_rk = 0.0
    DO m = 0, 6
      I1  = I1 + apar(m)*eta**REAL(m)
      I2  = I2 + bpar(m)*eta**REAL(m)
      I1_rk = I1_rk + apar(m)*REAL(m)*eta**REAL(m-1)*z3_rk + ap_rk(k,m)*eta**REAL(m)
      I2_rk = I2_rk + bpar(m)*REAL(m)*eta**REAL(m-1)*z3_rk + bp_rk(k,m)*eta**REAL(m)
    END DO
    
    ord1_rk  = 0.0
    ord2_rk  = 0.0
    order1   = 0.0
    order2   = 0.0
    DO i = 1,ncomp
      sig_ij(i,k) = 0.5 * ( dhs(i) + dhs(k) )
      uij(i,k)    = (1.0 - kij(i,k)) * SQRT( eps(i) * eps(k) )
      order1 = order1 + xx(i)*xx(k)* mseg(i)*mseg(k)*sig_ij(i,k)**3 * uij(i,k)/t
      order2 = order2 + xx(i)*xx(k)* mseg(i)*mseg(k)*sig_ij(i,k)**3 * (uij(i,k)/t)**2
      ord1_rk = ord1_rk + 2.0*mseg(k)*rho*xx(i)*mseg(i)*sig_ij(i,k)**3  *uij(i,k)/t
      ord2_rk = ord2_rk + 2.0*mseg(k)*rho*xx(i)*mseg(i)*sig_ij(i,k)**3 *(uij(i,k)/t)**2 
    END DO
    
    
    c1_con= 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3*z3)/zms**4   &
                    + (1.0 - m_mean)*(20.0*z3-27.0*z3*z3 +12.0*z3**3 -2.0*z3**4 )  &
                    /(zms*(2.0-z3))**2  )
    c2_con= - c1_con*c1_con *(  m_mean*(-4.0*z3*z3+20.0*z3+8.0)/zms**5   &
                                + (1.0 - m_mean) *(2.0*z3**3 +12.0*z3*z3-48.0*z3+40.0)  &
                                  /(zms*(2.0-z3))**3  )
    c1_rk= c2_con*z3_rk - c1_con*c1_con*m_rk(k)   *  ( (8.0*z3-2.0*z3*z3)/zms**4   &
           - (-2.0*z3**4 +12.0*z3**3 -27.0*z3*z3+20.0*z3) / (zms*(2.0-z3))**2  )
    
    mdsp(k) = -2.0*PI* ( order1*rho*rho*I1_rk + ord1_rk*I1 )  &
              -    PI* c1_con*m_mean * ( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
              -    PI* ( c1_con*m_rk(k) + c1_rk*m_mean ) * order2*rho*rho*I2
     
                            
End Do

!residual 
 ChemPot_res(1:ncomp) = mhs(1:ncomp) + mhc(1:ncomp) + mdsp(1:ncomp)                 ! /kT [-]
!total 
 ChemPot_tot(1:ncomp) = ChemPot_res(1:ncomp) + log( xx(1:ncomp)*rho ) ! /kT [-]

write(*,*)'rhob',xx(1:ncomp)*rho
 
End Subroutine Chemical_Potential









Subroutine PCSAFT_const(ap,bp)

Implicit None

!passed
REAL, INTENT(OUT) :: ap(0:6,3)
REAL, INTENT(OUT) :: bp(0:6,3)

! --- dispersion term constants ----------------------------------------
ap(0,1) =  0.91056314451539
ap(0,2) = -0.30840169182720
ap(0,3) = -0.09061483509767
ap(1,1) =  0.63612814494991
ap(1,2) =  0.18605311591713
ap(1,3) =  0.45278428063920
ap(2,1) =  2.68613478913903
ap(2,2) = -2.50300472586548
ap(2,3) =  0.59627007280101
ap(3,1) = -26.5473624914884
ap(3,2) =  21.4197936296668
ap(3,3) = -1.72418291311787
ap(4,1) =  97.7592087835073
ap(4,2) = -65.2558853303492
ap(4,3) = -4.13021125311661
ap(5,1) = -159.591540865600
ap(5,2) =  83.3186804808856
ap(5,3) =  13.7766318697211
ap(6,1) =  91.2977740839123
ap(6,2) = -33.7469229297323
ap(6,3) = -8.67284703679646

bp(0,1) =  0.72409469413165
bp(0,2) = -0.57554980753450
bp(0,3) =  0.09768831158356
bp(1,1) =  1.11913959304690  *2.0
bp(1,2) =  0.34975477607218  *2.0
bp(1,3) = -0.12787874908050  *2.0
bp(2,1) = -1.33419498282114  *3.0
bp(2,2) =  1.29752244631769  *3.0
bp(2,3) = -3.05195205099107  *3.0
bp(3,1) = -5.25089420371162  *4.0
bp(3,2) = -4.30386791194303  *4.0
bp(3,3) =  5.16051899359931  *4.0
bp(4,1) =  5.37112827253230  *5.0
bp(4,2) =  38.5344528930499  *5.0
bp(4,3) = -7.76088601041257  *5.0
bp(5,1) =  34.4252230677698  *6.0
bp(5,2) = -26.9710769414608  *6.0
bp(5,3) =  15.6044623461691  *6.0
bp(6,1) = -50.8003365888685  *7.0
bp(6,2) = -23.6010990650801  *7.0
bp(6,3) = -4.23812936930675  *7.0


End Subroutine PCSAFT_const



End Module mod_ChemPot
