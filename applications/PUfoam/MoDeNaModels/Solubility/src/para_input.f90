!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE para_input
!
! This subroutine provides pure component parameters and kij parameters.
! The following syntax applies:
!
! compna(i)                  component name
!
! parame(i,k)                pure comp. parameter:
!              parame(i,1):  segment number  [/]
!              parame(i,2):  segment diameter "sigma" [Angstrom]
!              parame(i,3):  segment energy param. epsilon/k [K]
!              parame(i,4):  model parameter; not used for PC-SAFT (=0)
!                            it is 10K most of the time for SAFT [K]
!              parame(i,5):  Param. for T-dependent segment diameter [/]
!              parame(i,6):  dipolar moment [debye]
!              parame(i,7):  quadrupolar moment [debye]
!              parame(i,8):  number of segments that are part of a branching 4-mer [/]
!              parame(i,9):
!              parame(i,10): ionic charge number (positiv or negativ) [/]
!              parame(i,11): polarizability [A**3]
!              parame(i,12): number of association sites [/]
!              parame(i,13): (=kap_hb, see below) [/]
!              parame(i,14 to 25): (=eps_hb, see below) [K]
! nhb_typ(i)                 number of different types of association sites (comp. i)
! nhb_no(i,k)                number of association sites of type k
! eps_hb                     depth of association potential [K]
! kap_hb                     effective width of assoc. potential (angle-averg.)
! mm                         molec. mass
! scaling                    param. for roughly scaling the set of objective functions
!
! As opposed to low-molec mass compounds, the molecular mass of a
! polymer is not obtained from this routine. Rather, it is a
! user-specification given in the file INPUT.INP
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE para_input

  USE BASIC_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  !-----------------------------------------------------------------------------

  IF (eos == 1) THEN
     CALL pcsaft_par
  ELSE
     write (*,*) 'Solubility Code: no pure component parameters defined for EOS=',eos
     stop 5
  END IF

  DO i = 1, ncomp
     IF ( mm(i) >= 1.0       .AND. mm(i) < 45.0  ) THEN
        scaling(i) = 10000.0
     ELSE IF( mm(i) >= 45.0  .AND. mm(i) < 90.0  ) THEN
        scaling(i) = 1000.0
     ELSE IF( mm(i) >= 90.0  .AND. mm(i) < 150.0 ) THEN
        scaling(i) = 100.0
     ELSE IF( mm(i) >= 150.0 .AND. mm(i) < 250.0 ) THEN
        scaling(i) = 10.0
     ELSE
        scaling(i) = 1.0
     END IF
     IF (parame(i,10) /= 0.0) scaling(i) = scaling(i) / 1.E4   ! Electrolytes
  END DO

END SUBROUTINE para_input


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pcsaft_par
!
! pure component parameters and kij parameters
! (as described in SUBROUTINE para_input)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE pcsaft_par

  USE BASIC_VARIABLES
  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, no
  INTEGER, DIMENSION(nc)                 :: nhb_typ
  INTEGER, DIMENSION(nc,nsite)           :: nhb_no
  REAL, DIMENSION(nc,nc,nsite,nsite)     :: eps_hb
  REAL, DIMENSION(nc,nc)                 :: kap_hb
  !-----------------------------------------------------------------------------


  DO  i = 1, ncomp
     parame(i,4) = 0.0     ! T correct. required for SAFT, not PC-SAFT
     parame(i,5) = 0.12    ! Param. for T-dependent segment diameter
     parame(i,6) = 0.0     ! dipolar moment
     parame(i,7) = 0.0     ! quadrupolar moment
     parame(i,8) = 0.0     ! number of segments that are part of a branching 4-mer
     parame(i,9) = 0.0
     parame(i,10)= 0.0     ! ionic charge number
     parame(i,11)= 0.0     ! polarizability
     lli(i)      = 0.0
     phi_criti(i)= 0.0
     chap(i)     = 0.0

     nhb_typ(i)  = 0
     kap_hb(i,i) = 0.0
     ! irgendwann sollten nhb_typ und kap_hb durch parame(i,12) und (i,13)
     ! ersetzt werden.

     SELECT CASE ( compna(i) )

     !MoDeNa components

    CASE ('14-butandiol')
      mm(i)       = 90.12
      parame(i,1) = 4.35923557
      parame(i,2) = 3.02947364
      parame(i,3) = 197.11998863
      
   CASE ('surfactant') 
     mm(i) =  2655.24078
     parame(i,1) = 78.5859962
     parame(i,2) = 4.17006833
     parame(i,3) = 230.284526
     parame(i,6) = 17.9645
      

    CASE('mdi')
      mm(i) =  2.50252E+02
      parame(i,1) = mm(i)*0.030769
      parame(i,2) = 2.886003
      parame(i,3) = 283.052778

    CASE ('air')
      mm(i)       = 28.899	!n2 and o2 according to mole fractions
      parame(i,1) = 1.18938	!n2 and o2 according to mole fractions (weighted artihm. avg)
      parame(i,2) = 3.28694	!n2 and o2 according to mole fractions (weighted artihm. avg)
      parame(i,3) = 95.672	!n2 and o2 according to mole fractions (weighted artihm. avg)

    CASE ('pu') !parameters obtained from MD simulations
!      mm(i) =  2042.22 !pu n = 5
!      parame(i,1) = mm(i)*0.008845
!      parame(i,2) = 5.680270
!      parame(i,3) = 497.997594
!      mm(i) =  340.37 !pu n = 0
!      parame(i,1) = mm(i)*0.043312
!      parame(i,2) = 3.008359
!      parame(i,3) = 273.445205
!     mm(i) =  680.74 !pu n = 1
!     parame(i,1) = mm(i)*0.024106
!     parame(i,2) = 3.744327
!     parame(i,3) = 321.486386
      mm(i) =  1021.11 !pu n = 2
      parame(i,1) = mm(i)*0.015076
      parame(i,2) = 4.537837
      parame(i,3) = 400.036950

   CASE ('tpg')
     mm(i) =  192.25
     parame(i,1) = mm(i)*0.01239
     parame(i,2) =  4.549
     parame(i,3) = 148.678
     parame(i,6) = 0.41

    nhb_typ(i)  = 2             ! no. of different association sites
    nhb_no(i,1) = 2            ! no. of sites of type 1
    nhb_no(i,2) = 2            ! no. of sites of type 2

    eps_hb(i,i,1,2)= 5597.844
    eps_hb(i,i,2,1)= 5597.844
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 0.03

   CASE ('polyol')      !parameters from homosegment GC-PC-SAFT

     mm(i) =  350.0
     parame(i,1) = 9.5746020799
     parame(i,2) =  3.5933857111
     parame(i,3) = 240.5792908908
     parame(i,6) = 2.7439527108

     nhb_typ(i)  = 2             ! no. of different association sites
     nhb_no(i,1) = 2            ! no. of sites of type 1
     nhb_no(i,2) = 2            ! no. of sites of type 2

     eps_hb(i,i,1,2)= 2517.00932
     eps_hb(i,i,2,1)= 2517.00932
     eps_hb(i,i,1,1)= 0.0
     eps_hb(i,i,2,2)= 0.0
     kap_hb(i,i)    = 0.0068253298

     CASE ('ps')
        parame(i,1) = mm(i)*1.9E-2
        parame(i,2) = 4.10705961
        parame(i,3) = 267.0
     CASE ('pg2') !Polyglycerol 2
        mm(i) = 2000.0
        parame(i,1) = mm(i)*2.37E-2 ! from figure 5 PCSAFT paper
        parame(i,2) = 3.8           ! from figure 5 PCSAFT paper
        parame(i,3) = 270.0         ! starting value for iteration
        ! this is the extra parameter
        parame(i,8) = mm(i)*2.37E-2

        nhb_typ(i)  = 2             ! no. of different association sites
        nhb_no(i,1) = 27            ! no. of sites of type 1
        nhb_no(i,2) = 27            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2544.6     ! taken from butanol (same M/OH)
        eps_hb(i,i,2,1)= 2544.6
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)= .00489087833   ! taken from butanol (same M/OH)
     CASE ('peva')
        parame(i,1) = mm(i)*2.63E-2
        ! --   0 Gew.% VA-------------
        ! parame(i,2) = 4.021767
        ! parame(i,3) = 249.5
        ! -- 7.5 Gew.% VA-------------
        ! parame(i,2) = 4.011
        ! parame(i,3) = 248.1864
        ! parame(i,3) = 247.6286
        ! ---12.7 Gew.% VA------------
        ! parame(i,2) = 4.0028
        ! parame(i,3) = 247.2075
        ! parame(i,3) = 246.24454
        ! ---27.3 Gew.% VA------------
        ! parame(i,2) = 3.9762
        ! parame(i,3) = 244.114
        ! parame(i,3) = 241.9345
        ! ---31.8 Gew.% VA------------
        parame(i,2) = 3.9666
        parame(i,3) = 243.0436
        ! parame(i,3) = 240.46
        ! ---42.7 Gew.% VA------------
        ! parame(i,2) = 3.9400
        ! parame(i,3) = 240.184
        ! parame(i,3) = 236.62
        ! ---------------
     CASE ('pp')
        parame(i,1) = mm(i)*2.2E-2
        parame(i,2) = 4.2
        parame(i,3) = 220.0

        parame(i,1) = mm(i)*0.0230487701
        parame(i,2) = 4.1
        parame(i,3) = 217.0
     CASE ('pe')     ! hdpe
        parame(i,1) = mm(i)*2.622E-2
        parame(i,2) = 4.021767
        parame(i,3) = 252.0
        ! HDPE: extrapolated from pure comp. param. of n-alkane series!
        ! parame(i,1) = mm(i)*2.4346E-2
        ! parame(i,2) = 4.07182
        ! parame(i,3) = 269.67
        !!  parame(i,3) = 252.5
     CASE ('ldpe')
        parame(i,1) = mm(i)*2.63E-2
        parame(i,2) = 4.021767
        parame(i,3) = 249.5
     CASE ('pba')
        parame(i,1) = mm(i)*2.5872E-2
        parame(i,2) = 3.95
        parame(i,3) = 229.0
     CASE ('dextran')
        parame(i,1) = mm(i)*2.E-2
        parame(i,2) = 4.0
        parame(i,3) = 300.0
     CASE ('glycol-ethers')
        ! mm(i) = 218.0
        ! parame(i,1) = 7.4044
        ! parame(i,2) = 3.61576
        ! parame(i,3) = 244.0034598
        mm(i) = 222.0
        parame(i,1) = 7.994
        parame(i,2) = 3.445377778
        parame(i,3) = 234.916506
     CASE ('LJ')
        mm(i)       = 1.0
        parame(i,1)  = 4.0
        parame(i,2) = 1.0
        parame(i,3) = 1.0
     CASE ('LJ1205')
        mm(i)       = 1.0
        parame(i,1)  = 1.0
        parame(i,2) = 1.0
        parame(i,3) = 140.0
     CASE ('adamantane')
        mm(i)       =   136.235000000000
        parame(i,1) =   4.81897145432221
        parame(i,2) =   3.47128575274660
        parame(i,3) =   266.936967922521
     CASE ('methane')
        mm(i)       = 16.043
        parame(i,1)  = 1.0
        parame(i,2) = 3.70388767
        parame(i,3) = 150.033987
        ! LLi(i)      = 1.185*parame(i,2)
        ! phi_criti(i)= 11.141
        ! chap(i)     = 0.787
        lli(i)      = 1.398*parame(i,2)
        phi_criti(i)= 16.01197
        chap(i)     = 0.6
        IF (pol == 2) parame(i,11)= 2.593
        ! --- adjusted to Tc, Pc und omega ---
        ! mm(i)       = 16.0430000000000
        ! parame(i,1) = 1.03353666429362
        ! parame(i,2) = 3.64824920605089
        ! parame(i,3) = 147.903953522994
        lli(i)      = 2.254442763775*parame(i,2)
        phi_criti(i)= 42.060975627454
        chap(i)     = 0.704895924
        lli(i)      = 1.935801125833*parame(i,2)
        phi_criti(i)= 26.363325937261
        chap(i)     = 0.700112854298
        lli(i)      = 2.610103087662*parame(i,2)
        phi_criti(i)= 38.192854403173
        chap(i)     = 0.812100472735
        ! 2.122960316503   34.937141524804    0.734513223627
        ! 2.082897379591   33.036391564859    0.877578492999
     CASE ('ethane')
        mm(i)       = 30.070
        parame(i,1) = 1.60684
        parame(i,2) = 3.52059
        parame(i,3) = 191.4238
        lli(i)      = 1.40*parame(i,2)
        phi_criti(i)= 15.38
        chap(i)     = 0.520
        IF (pol == 2) parame(i,11)= 4.3
        ! --- adjusted to Tc, Pc und omega ---
        ! mm(i)       = 30.069
        ! parame(i,1) = 1.74034548122
        ! parame(i,2) = 3.4697441893134
        ! parame(i,3) = 181.90770083591
        IF (pol >= 1) mm(i)       = 30.0700000000000
        IF (pol >= 1) parame(i,1) = mm(i)*  5.341907666260094E-002
        IF (pol >= 1) parame(i,2) = 3.52104466654628
        IF (pol >= 1) parame(i,3) = 191.449300423694
        IF (pol >= 1) parame(i,7) = 0.650000000000000
        IF (pol >= 1) lli(i)      = 0.0
        IF (pol >= 1) phi_criti(i)= 0.0
        IF (pol >= 1) chap(i)     = 0.0
     CASE ('propane')
        mm(i)       = 44.096
        parame(i,1) = mm(i)*  .0453970622
        parame(i,2) = 3.61835302
        parame(i,3) = 208.110116
        lli(i)      = 1.8*parame(i,2)
        phi_criti(i)= 21.0
        chap(i)     = 1.0
        lli(i)      = 1.63*parame(i,2)
        phi_criti(i)= 20.37
        chap(i)     = 0.397
        IF (pol == 2) parame(i,11)= 6.29
     CASE ('butane')
        mm(i)       = 58.123
        parame(i,1) = mm(i)*  .0401146927
        parame(i,2) = 3.70860139
        parame(i,3) = 222.877405
        lli(i)      = 1.75*parame(i,2)
        phi_criti(i)= 23.43
        chap(i)     = 0.304
        ! LLi(i)      = 1.942079633622*parame(i,2)
        ! phi_criti(i)= 24.527323443155
        ! chap(i)     = 0.734064026277
        ! LLi(i)      = 1.515115760477*parame(i,2)
        ! phi_criti(i)= 17.682929717796
        ! chap(i)     = 0.335848717079
        IF (pol == 2) parame(i,11)= 8.2
        ! --- adjusted to Tc, Pc und omega ---
        ! mm(i)       = 58.1230000000
        ! parame(i,1) = 2.45352304112
        ! parame(i,2) = 3.74239117802
        ! parame(i,3) = 214.185157925
     CASE ('pentane')
        mm(i)       = 72.146
        parame(i,1) = mm(i)*  .03727896
        parame(i,2) = 3.77293174
        parame(i,3) = 231.197015
        IF (pol == 2) parame(i,11)= 9.99
     CASE ('hexane')
        mm(i)       = 86.177
        parame(i,1) = mm(i)*  .0354812325
        parame(i,2) = 3.79829291
        parame(i,3) = 236.769054
        lli(i)      = 2.24*parame(i,2)
        phi_criti(i)= 33.25
        chap(i)     = 0.205
        IF (pol == 2) parame(i,11)= 11.9
     CASE ('heptane')
        mm(i)       = 100.203
        parame(i,1) = mm(i)*  .034762384
        parame(i,2) = 3.80487025
        parame(i,3) = 238.400913
        lli(i)      = 2.35*parame(i,2)
        phi_criti(i)= 38.10
        chap(i)     = 0.173
        IF (pol == 2) parame(i,11)= 13.61
     CASE ('octane')
        mm(i)       = 114.231
        parame(i,1) = mm(i)*  .0334228038
        parame(i,2) = 3.83732677
        parame(i,3) = 242.775853
        ! LLi(i)      = 2.0*parame(i,2)
        ! phi_criti(i)= 18.75
        ! chap(i)     = 1.0
        lli(i)      = 2.63*parame(i,2)
        phi_criti(i)= 42.06
        chap(i)     = 0.155
        IF (pol == 2) parame(i,11)= 15.9
     CASE ('nonane')
        mm(i)       = 128.25
        parame(i,1) = mm(i)*  .0328062594
        parame(i,2) = 3.84483643
        parame(i,3) = 244.508457
     CASE ('decane')
        mm(i)       = 142.285
        parame(i,1) = mm(i)*  .03277373
        parame(i,2) = 3.8384498
        parame(i,3) = 243.866074
        lli(i)      = 1.845*parame(i,2)
        phi_criti(i)= 21.27
        chap(i)     = 1.0
        lli(i)      = 2.68*parame(i,2)
        phi_criti(i)= 45.0
        chap(i)     = 0.15
        IF (pol == 2) parame(i,11)= 19.1
        ! --- adjusted to Tc, Pc und omega ---
        !  parame(i,1) = 4.794137228322
        !  parame(i,2) = 4.030446690586
        !  parame(i,3) = 236.5884493386
     CASE ('dodecane')
        mm(i)       = 170.338
        parame(i,1) = mm(i)*  .0311484156
        parame(i,2) = 3.89589236
        parame(i,3) = 249.214532
     CASE ('hexadecane')
        mm(i)       = 226.446
        parame(i,1) = mm(i)*  .0293593045
        parame(i,2) = 3.95516743
        parame(i,3) = 254.700131
     CASE ('pentadecane')
        mm(i)       = 212.419
        parame(i,1) = 6.2855
        parame(i,2) = 3.9531
        parame(i,3) = 254.14
     CASE ('octadecane')
        mm(i)       = 254.5
        parame(i,1) = 7.3271
        parame(i,2) = 3.9668
        parame(i,3) = 256.20
        IF (pol == 2) parame(i,11)= 30.2
        ! --- adjusted to Tc, Pc und omega ---
        !  mm(i)       = 226.446000000000
        !  parame(i,1) = 6.66976520488694
        !  parame(i,2) = 4.25025597912511
        !  parame(i,3) = 249.582941976119
     CASE ('eicosane')
        mm(i)       = 282.553
        parame(i,1) = mm(i)*  .0282572812
        parame(i,2) = 3.98692612
        parame(i,3) = 257.747939
     CASE ('triacontane')
        ! mm(i)       = 422.822             ! polyethylene parameters
        ! parame(i,1) = mm(i)*2.622E-2
        ! parame(i,2) = 4.021767
        ! parame(i,3) = 252.0
        mm(i)       = 422.822             ! param. by extrapolation of n-alkanes
        parame(i,1) = mm(i)*  0.026922527
        parame(i,2) = 4.007608009
        parame(i,3) = 262.28622
     CASE ('octaeicosane')
        mm(i)       = 395.0             ! param. by extrapolation of n-alkanes (sloppy!!)
        parame(i,1) = mm(i)*  0.026922527
        parame(i,2) = 4.007608009
        parame(i,3) = 262.28622
     CASE ('tetracontane')
        ! mm(i)       = 563.1             ! polyethylene parameters
        ! parame(i,1) = mm(i)*2.622E-2
        ! parame(i,2) = 4.021767
        ! parame(i,3) = 252.0
        mm(i)       = 563.1             ! param. by extrapolation of n-alkanes
        parame(i,1) = mm(i)*0.026287593
        parame(i,2) = 4.023277
        parame(i,3) = 264.10466
     CASE ('isobutane')
        mm(i)       = 58.123
        parame(i,1) = mm(i)*  .0389105395
        parame(i,2) = 3.75735249
        parame(i,3) = 216.528584
     CASE ('methylbutane')      ! isopentane
        mm(i)       = 72.15
        parame(i,1) = 2.5620
        parame(i,2) = 3.8296
        parame(i,3) = 230.75
     CASE ('dimethylpropane')   ! neopentane
        mm(i)       = 72.15
        parame(i,1) = 2.3543
        parame(i,2) = 3.9550
        parame(i,3) = 225.69
     CASE ('2-methylpentane')
        mm(i)       = 86.177
        parame(i,1) = mm(i)*  .0340166994
        parame(i,2) = 3.85354665
        parame(i,3) = 235.5801
     CASE ('23-dimethylbutane')
        mm(i)       = 86.177
        parame(i,1) = mm(i)*  .0311599207
        parame(i,2) = 3.9544545
        parame(i,3) = 246.068188
     CASE ('ethylene')
        mm(i)       = 28.05
        parame(i,1) = mm(i)*  .0567939013
        parame(i,2) = 3.44499904
        parame(i,3) = 176.468725
        IF (pol >= 1) parame(i,1) = 1.5477      ! PCP-SAFT
        IF (pol >= 1) parame(i,2) = 3.4475
        IF (pol >= 1) parame(i,3) = 179.19
        IF (pol >= 1) parame(i,7) = 1.9155
        IF (pol == 2) parame(i,11)= 4.252
     CASE ('propylene')
        mm(i)       = 42.081
        parame(i,1) = mm(i)*  .0465710324
        parame(i,2) = 3.53559831
        parame(i,3) = 207.189309
        ! --- adjusted to Tc, Pc und omega ---
        ! mm(i)       = 42.081
        ! parame(i,1) = 2.086735327675
        ! parame(i,2) = 3.536779407969
        ! parame(i,3) = 198.3529810625
     CASE ('1-butene')
        mm(i)       = 56.107
        ! parame(i,1) = mm(i)*  0.40752372E-01
        parame(i,1) = 2.28649
        parame(i,2) = 3.64306
        parame(i,3) = 222.003
        IF (pol == 2) parame(i,11)= 7.97
     CASE ('1-pentene')
        mm(i)       = 70.134
        parame(i,1) = 2.6006
        parame(i,2) = 3.7399
        parame(i,3) = 231.99
     CASE ('1-hexene')
        mm(i)       = 84.616
        parame(i,1) = mm(i)*  .0352836857
        parame(i,2) = 3.77529612
        parame(i,3) = 236.810973
     CASE ('1-octene')
        mm(i)       = 112.215
        parame(i,1) = mm(i)*  .033345175
        parame(i,2) = 3.81329011
        parame(i,3) = 243.017587
     CASE ('2-butyne')
        mm(i)       = 54.090
        parame(i,1) = 2.59601041
        parame(i,2) = 3.34635221
        parame(i,3) = 237.77296813
        IF (pol >= 1) parame(i,1) = 2.877971
        IF (pol >= 1) parame(i,2) = 3.22362803
        IF (pol >= 1) parame(i,3) = 161.53695842
        IF (pol >= 1) parame(i,6) = 3.957228
     CASE ('cyclopentane')
        mm(i)       = 70.13
        parame(i,1) = mm(i)*  .0337262571
        parame(i,2) = 3.71139254
        parame(i,3) = 265.828755
     CASE ('cyclohexane')
        mm(i)       = 84.147
        parame(i,1) = mm(i)*  .0300695505
        parame(i,2) = 3.84990887
        parame(i,3) = 278.108786
        IF (pol == 2) parame(i,11)= 10.87
     CASE ('toluene')
        mm(i)       = 92.141
        parame(i,1) = mm(i)*  .0305499338
        parame(i,2) = 3.71689689
        parame(i,3) = 285.68996
        IF (pol == 2) parame(i,11)= 11.8
        ! --- adjusted to Tc, Pc und omega ---
        ! mm(i)       = 92.141
        ! parame(i,1) = 3.002119827762
        ! parame(i,2) = 3.803702734224
        ! parame(i,3) = 271.9428642880
     CASE ('m-xylene')
        mm(i)       = 106.167
        parame(i,1) = mm(i)*  .030011086
        parame(i,2) = 3.75625585
        parame(i,3) = 283.977525
     CASE ('o-xylene')
        mm(i)     =           106.167
        parame(i,1)   = mm(i)*  .0295409161
        parame(i,2) =           3.76000631
        parame(i,3) =           291.049123
     CASE ('thf')
        mm(i)       =   72.1057000000000
        ! parame(i,1) = mm(i)*  0.34311391E-01
        parame(i,1) =   2.47404685540709
        parame(i,2) =   3.51369375633677
        parame(i,3) =   274.181927093696
        parame(i,6) =   1.63100000000000
     CASE ('co2')
        mm(i)       = 44.01
        parame(i,1) = mm(i)*  .0470968503
        parame(i,2) = 2.7851954
        parame(i,3) = 169.207418
        IF (pol >= 1) parame(i,1) = mm(i)*  3.438191426159075E-002
        IF (pol >= 1) parame(i,2) = 3.18693935424469
        IF (pol >= 1) parame(i,3) = 163.333232725156
        IF (pol >= 1) parame(i,7) = 4.400000000000
        IF (pol >= 1) lli(i)      = 1.472215*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 17.706567
        IF (pol >= 1) chap(i)     = 0.5
        IF (pol == 2) parame(i,11)= 2.911
     CASE ('co')
        IF (pol /= 1) write (*,*) 'parameters for co missing'
        IF (pol /= 1) stop 5
        IF (pol >= 1) mm(i)       = 28.01
        IF (pol >= 1) parame(i,1) = mm(i)*  5.126059746332587E-002  ! 1.43580933494776
        IF (pol >= 1) parame(i,2) = 3.13556624711756
        IF (pol >= 1) parame(i,3) = 87.7191028693595
        IF (pol >= 1) parame(i,6) = 0.1098
     CASE ('n2')
        mm(i)       = 28.01
        parame(i,1) = mm(i)*  .0430301713
        parame(i,2) = 3.3129702
        parame(i,3) = 90.9606924
        IF (pol >= 1) parame(i,1) = mm(i)*  3.971157114787596E-002
        IF (pol >= 1) parame(i,2) = 3.42116853868336
        IF (pol >= 1) parame(i,3) = 92.3972606842862
        IF (pol >= 1) parame(i,7) = 1.52000000000000
        IF (pol >= 1) lli(i)      = 1.5188*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 19.9247
        IF (pol >= 1) chap(i)     = 0.375
        ! better RGT-results came later, with: 1.5822   21.201   0.3972
     CASE ('o2')
        mm(i)       = 32.05
        parame(i,1) = mm(i)*  .0353671563
        parame(i,2) = 3.19465166
        parame(i,3) = 114.430197
     CASE ('hydrogen')
        mm(i)       = 2.016
        parame(i,1) = mm(i)*  .258951975
        parame(i,2) = 4.43304935
        parame(i,3) = 29.6509579

        mm(i)       = 2.016
        parame(i,1) = 1.0
        parame(i,2) = 2.915
        parame(i,3) = 38.0

        ! mm(i)       = 2.016  !  Ghosh et al. 2003
        ! parame(i,1) = 1.0
        ! parame(i,2) = 2.986
        ! parame(i,3) = 19.2775
     CASE ('argon')
        ! mm(i)       = 39.948  ! adjusted m !!
        ! parame(i,1) = 0.9285
        ! parame(i,2) = 3.4784
        ! parame(i,3) = 122.23
        mm(i)       = 39.948 ! enforced m=1 !!
        parame(i,1) = 1.0
        parame(i,2) = 3.3658
        parame(i,3) = 118.34
        IF (pol == 2) parame(i,11)= 1.6411
     CASE ('xenon')
        mm(i)       =   131.29
        parame(i,1) =   1.0
        parame(i,2) =   3.93143
        parame(i,3) =   227.749
     CASE ('chlorine')  ! Cl2
        mm(i)       =   70.906
        parame(i,1) =   1.5514
        parame(i,2) =   3.3672
        parame(i,3) =   265.67
     CASE ('SF6')
        mm(i)       = 146.056  ! adjusted m !!
        parame(i,1) = 2.48191
        parame(i,2) = 3.32727
        parame(i,3) = 161.639
        ! mm(i)       = 146.056 ! enforced m=1 !!
        ! parame(i,1) = 1.0
        ! parame(i,2) = 4.55222
        ! parame(i,3) = 263.1356
     CASE ('benzene')
        mm(i)       = 78.114
        parame(i,1) = mm(i)*  .0315590546
        parame(i,2) = 3.64778975
        parame(i,3) = 287.354574
        IF (pol >= 1) mm(i)       = 78.114   ! PCP-SAFT with m=2 in QQ term
        IF (pol >= 1) parame(i,1) = mm(i)*  2.932783311E-2 !  = 2.29091435590515
        IF (pol >= 1) parame(i,2) =         3.7563854
        IF (pol >= 1) parame(i,3) =         294.06253
        IF (pol >= 1) parame(i,7) = 5.5907
     CASE ('ethylbenzene')
        mm(i)       = 106.167
        parame(i,1) = mm(i)*  .0290120497
        parame(i,2) = 3.79741116
        parame(i,3) = 287.348098
        IF (pol == 2) parame(i,11)= 13.3
     CASE ('propylbenzene')
        mm(i)     =           120.194
        parame(i,1)   = mm(i)*  .0278171627
        parame(i,2) =           3.8437772
        parame(i,3) =           288.128269
     CASE ('n-butylbenzene')
        mm(i)       = 134.221
        parame(i,1) = mm(i)*  .0280642225
        parame(i,2) = 3.87267961
        parame(i,3) = 283.072331
     CASE ('tetralin')
        mm(i)       = 132.205
        parame(i,1) = mm(i)*  .0250640795
        parame(i,2) = 3.87498866
        parame(i,3) = 325.065688
     CASE ('methylcyclohexane')
        mm(i)       = 98.182
        parame(i,1) = mm(i)*  .0271259953
        parame(i,2) = 3.99931892
        parame(i,3) = 282.334148
        IF (pol == 2) parame(i,11)= 13.1
     CASE ('methylcyclopentane')
        mm(i)       = 84.156
        parame(i,1) = mm(i)*  .0310459009
        parame(i,2) = 3.82534693
        parame(i,3) = 265.122799
     CASE ('acetone')
        mm(i)       = 58.0800000000000     ! PC-SAFT
        parame(i,1) = mm(i)*  4.870380408159182E-002  ! =2.82871694105885
        parame(i,2) = 3.24969003020675
        parame(i,3) = 250.262241927379
        lli(i)      = 2.0021*parame(i,2)
        phi_criti(i)= 21.336
        chap(i)     = 0.24931
        IF (pol >= 1) mm(i)       = 58.0800000000000     ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  4.725811736856114E-002  ! =2.74475145676603
        IF (pol >= 1) parame(i,2) = 3.27423145271184
        IF (pol >= 1) parame(i,3) = 232.990879135326
        IF (pol >= 1) parame(i,6) = 2.88000000000000
        IF (pol >= 1) lli(i)      = 2.0641*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 28.1783
        IF (pol >= 1) chap(i)     = 0.22695
        IF (pol >= 2) mm(i)       = 58.0800000000000     ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  4.902301475689938E-002  !  =2.84725669708072
        IF (pol >= 2) parame(i,2) = 3.23880349104868
        IF (pol >= 2) parame(i,3) = 220.884202656054
        IF (pol >= 2) parame(i,6) = 2.88000000000000
        IF (pol == 2) parame(i,11)= 6.40000000000000
        ! IF (pol >= 1) mm(i)       =   58.0800000000000     ! PCP-SAFT (with more emphasis on crit. pt.)
        ! IF (pol >= 1) parame(i,1) =   2.79925169706703
        ! IF (pol >= 1) parame(i,2) =   3.26496723061544
        ! IF (pol >= 1) parame(i,3) =   221.426206364920
        ! IF (pol >= 1) parame(i,6) =   3.22546611481946
     CASE ('butanone')
        mm(i)       = 72.1066               !  PC-SAFT
        parame(i,1) = mm(i)*  4.264192830122321E-002  !  =3.07476446724498
        parame(i,2) = 3.39324011060028
        parame(i,3) = 252.267273608975
        IF (pol >= 1) mm(i)       = 72.1066               !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  4.137668924230600E-002  !  =2.98353238051926
        IF (pol >= 1) parame(i,2) = 3.42393701353423
        IF (pol >= 1) parame(i,3) = 244.994381354681
        IF (pol >= 1) parame(i,6) = 2.78000000000000
        IF (pol >= 2) mm(i)       = 72.1066               !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  4.254697075199448E-002  !  =3.06791740122577
        IF (pol >= 2) parame(i,2) = 3.39138375903252
        IF (pol >= 2) parame(i,3) = 236.527763837528
        IF (pol >= 2) parame(i,6) = 2.78000000000000
        IF (pol == 2) parame(i,11)= 8.13000000000000
     CASE ('2-pentanone')
        ! mm(i)       = 86.134                !  PC-SAFT
        ! parame(i,1) = mm(i)*  3.982654501296355E-002  !  =3.43041962814660
        ! parame(i,2) = 3.46877976946838
        ! parame(i,3) = 249.834724442656
        ! mm(i)       = 86.134                !  PCP-SAFT
        ! parame(i,1) = mm(i)*  3.893594769994072E-002  !  =3.35370891918669
        ! parame(i,2) = 3.49417356096593
        ! parame(i,3) = 246.656329096835
        ! parame(i,6) = 2.70000000000000
        mm(i)       = 86.134                !  PCIP-SAFT
        parame(i,1) = mm(i)*  3.973160761515879E-002  !  =3.42224229032409
        parame(i,2) = 3.46827593107280
        parame(i,3) = 240.904278156822
        parame(i,6) = 2.70000000000000
        IF (pol == 2) parame(i,11)= 9.93000000000000
     CASE ('3-pentanone')
        mm(i)       = 86.134                !  PC-SAFT
        parame(i,1) = 3.36439508013322
        parame(i,2) = 3.48770251979329
        parame(i,3) = 252.695415552376
        IF (pol >= 1) mm(i)       = 86.134                !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 3.27863398611842
        IF (pol >= 1) parame(i,2) = 3.51592571835030
        IF (pol >= 1) parame(i,3) = 248.690775540981
        IF (pol >= 1) parame(i,6) = 2.82000000000000
        IF (pol == 2) mm(i)       = 86.134                !  PCIP-SAFT
        IF (pol == 2) parame(i,1) = 3.34821857026283
        IF (pol == 2) parame(i,2) = 3.48903345340516
        IF (pol == 2) parame(i,3) = 242.314578558329
        IF (pol == 2) parame(i,6) = 2.82000000000000
        IF (pol == 2) parame(i,11)= 9.93000000000000
     CASE ('2-octanone')
        mm(i)       =   128.212040000000
        parame(i,1) =   4.32004934535641
        parame(i,2) =   3.68827454027963
        parame(i,3) =   257.552579277391
     CASE ('cyclohexanone')          ! from DIPPR
        ! IF (pol.GE.1) mm(i)       = 98.1430                !  PCP-SAFT
        ! IF (pol.GE.1) parame(i,1) = 3.084202
        ! IF (pol.GE.1) parame(i,2) = 3.613681
        ! IF (pol.GE.1) parame(i,3) = 286.15865
        ! IF (pol.GE.1) parame(i,6) = 3.087862
        IF (pol >= 1) mm(i)       = 98.1500000000000
        IF (pol >= 1) parame(i,1) = 2.72291913132818
        IF (pol >= 1) parame(i,2) = 3.79018433908522
        IF (pol >= 1) parame(i,3) = 314.772193827344
        IF (pol >= 1) parame(i,6) = 3.24600000000000
        IF (pol /= 1) WRITE (*,*) 'no non-polar param. for cyclohexanone'
        IF (pol /= 1) STOP
     CASE ('propanal')
        mm(i)       = 58.08                             !  PC-SAFT
        parame(i,1) = 2.67564746980910
        parame(i,2) = 3.26295953984941
        parame(i,3) = 251.888982765626
        IF (pol >= 1) mm(i)       = 58.08               !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.60007872084995
        IF (pol >= 1) parame(i,2) = 3.28720732189761
        IF (pol >= 1) parame(i,3) = 235.205188090107
        IF (pol >= 1) parame(i,6) = 2.72000000000000
        IF (pol >= 2) mm(i)       = 58.08               !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.72471167411028
        IF (pol >= 2) parame(i,2) = 3.24781643022922
        IF (pol >= 2) parame(i,3) = 221.642071811094
        IF (pol >= 2) parame(i,6) = 2.72000000000000
        IF (pol >= 2) parame(i,11)= 6.50000000000000
     CASE ('butanal')
        mm(i)       = 72.1066000000000                  ! PC-SAFT
        parame(i,1) = 2.96824823599784
        parame(i,2) = 3.44068916025889
        parame(i,3) = 253.929404992884
        IF (pol >= 1) mm(i)       = 72.1066000000000    ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.86783706423953
        IF (pol >= 1) parame(i,2) = 3.47737904036296
        IF (pol >= 1) parame(i,3) = 247.543312127310
        IF (pol >= 1) parame(i,6) = 2.72000000000000
     CASE ('dmso')
        mm(i)       = 78.1300000000000                  ! PC-SAFT
        parame(i,1) = 2.92225114054231
        parame(i,2) = 3.27780791606297
        parame(i,3) = 355.688793038512
        IF (pol >= 1) mm(i)       = 78.1300000000000    ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = 3.02433694138348
        IF (pol >= 1) parame(i,2) = 3.24270742566613
        IF (pol >= 1) parame(i,3) = 309.357476696679
        IF (pol >= 1) parame(i,6) = 3.96000000000000
        IF (pol >= 2) mm(i)       = 78.1300000000000    ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 3.19078234633277
        IF (pol >= 2) parame(i,2) = 3.19778269816832
        IF (pol >= 2) parame(i,3) = 286.337981216861
        IF (pol >= 2) parame(i,6) = 3.96000000000000
        IF (pol >= 2) parame(i,11)= 7.97000000000000
     CASE ('acetone_JC')  ! Jog-Chapman
        ! mm(i)       = 58.0800000000000   ! Dominik et al.2005
        ! parame(i,1) = 2.221
        ! parame(i,2) = 3.607908
        ! parame(i,3) = 259.99
        ! parame(i,6) = 2.7
        ! parame(i,8) = 0.2258
        ! mm(i)       = 58.0800000000000
        ! parame(i,1) = mm(i)*  3.556617369195472E-002
        ! parame(i,2) = 3.58780367502515
        ! parame(i,3) = 273.025100470307
        ! parame(i,6) = 2.70000000000000
        ! parame(i,8) = 0.229800000000000

        mm(i)       = 58.08     !  Tumakaka Sadowski 2004
        parame(i,1) = mm(i)*  3.766E-2
        parame(i,2) = 3.6028
        parame(i,3) = 245.49
        parame(i,6) = 2.72
        parame(i,8) = 0.2969
        ! mm(i)       = 58.0800000000000     ! no adjust. DD-param.
        ! parame(i,1) = 1.87041620247774
        ! parame(i,2) = 3.79783535570774
        ! parame(i,3) = 208.885730881588
        ! parame(i,6) = 2.88000000000000
        ! parame(i,8) = 1.0/parame(i,1)
        WRITE (*,*) 'caution: parame(i,8) is now used for branching'
        STOP
        kij(1,2) = -0.005
        kij(2,1)=kij(1,2)
     CASE ('acetone_SF')  ! Saager-Fischer
        mm(i)       = 58.08
        parame(i,1) = mm(i)*  4.603296414764944E-002
        parame(i,2) = 3.29454924451643
        parame(i,3) = 221.052649057645
        parame(i,6) = 2.70000000000000
        parame(i,8) = 0.625410000000000
        mm(i)       = 58.08 ! form as expected from me - no DD-param adjusted.dat
        parame(i,1) = mm(i)*  4.364264724158790E-002  !  =2.53476495179143
        parame(i,2) = 3.37098670735567
        parame(i,3) = 254.366379701851
        parame(i,6) = 2.88000000000000
        ! mm(i)       = 58.08  ! form as expected but w/ sumseg/1.5 - no DD-param adjusted.dat
        ! parame(i,1) = mm(i)*  4.694644361257521E-002  !  =2.72664944501837
        ! parame(i,2) = 3.27842292595463
        ! parame(i,3) = 238.398883501772
        ! parame(i,6) = 2.88000000000000
        ! mm(i)       = 58.08  ! form as expected but w/ sumseg/1.5 and fdd*sumseg- no DD-param adjusted.dat
        ! parame(i,1) = mm(i)*  4.458214655521766E-002  !  =2.58933107192704
        ! parame(i,2) = 3.32050824493493
        ! parame(i,3) = 218.285994651271
        ! parame(i,6) = 2.88000000000000
        WRITE (*,*) 'caution: parame(i,8) is now used for branching'
        STOP
        kij(1,2) = 0.035
        kij(2,1)=kij(1,2)
     CASE ('ethylacetate_JC')  ! Jog-Chapman
        ! mm(i)       = 88.11
        ! parame(i,1) = 2.7481
        ! parame(i,2) = 3.6511
        ! parame(i,3) = 236.99
        ! parame(i,6) = 1.84
        ! parame(i,8) = 0.5458
        mm(i)       = 88.1060000000000
        parame(i,1) = mm(i)*  0.03117 ! 2.74626402
        parame(i,2) = 3.6493
        parame(i,3) = 236.75
        parame(i,6) = 1.8
        parame(i,8) = 0.5462
     CASE ('ethylacetate_SF')  ! Saager-Fischer
        mm(i)       = 88.106
        parame(i,1) = mm(i)*  3.564165384763394E-002
        parame(i,2) = 3.447379322
        parame(i,3) = 226.0930487
        parame(i,6) = 1.8
        parame(i,8) = 0.849967000000000
        WRITE (*,*) 'caution: parame(i,8) is now used for branching'
        STOP
     CASE ('12po_JC')  ! Jog-Chapman
        mm(i)       = 58.08
        parame(i,1) = 2.0105
        parame(i,2) = 3.6095
        parame(i,3) = 258.82
        parame(i,6) = 2.0
        parame(i,8) = 0.3979
        WRITE (*,*) 'caution: parame(i,8) is now used for branching'
        STOP
     CASE ('12po_SF')  ! Saager-Fischer
        mm(i)       = 58.08
        parame(i,1) = 2.1341
        parame(i,2) = 3.4739
        parame(i,3) = 252.95
        parame(i,6) = 2.0
        parame(i,8) = 0.916
        WRITE (*,*) 'caution: parame(i,8) is now used for branching'
        STOP
     CASE ('acrylonitrile')
        IF (pol >= 2) mm(i)       = 53.06     ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.168
        IF (pol >= 2) parame(i,2) = 3.575
        IF (pol >= 2) parame(i,3) = 214.83
        IF (pol >= 2) parame(i,6) = 3.91
        IF (pol == 2) parame(i,11)= 8.04
        IF (pol >= 2) mm(i)       = 53.0000000000000  ! second parameter set ??
        IF (pol >= 2) parame(i,1) = 2.45403467006041
        IF (pol >= 2) parame(i,2) = 3.41276825781723
        IF (pol >= 2) parame(i,3) = 195.194353082408
        IF (pol >= 2) parame(i,6) = 3.91000000000000
        IF (pol == 2) parame(i,11)= 8.04000000000000
     CASE ('butyronitrile')
        ! mm(i)       = 69.11
        ! parame(i,1) = 2.860
        ! parame(i,2) = 3.478
        ! parame(i,3) = 253.21
        ! parame(i,6) = 4.07
        mm(i)       = 69.11
        parame(i,1) = 2.989
        parame(i,2) = 3.441
        parame(i,3) = 234.04
        parame(i,6) = 4.07
        IF (pol == 2) parame(i,11)= 8.4
     CASE ('propionitrile')
        mm(i)       = 55.079                !  PC-SAFT
        parame(i,1) = 2.66211021227108
        parame(i,2) = 3.34032231132738
        parame(i,3) = 294.078737359580
        IF (pol >= 1) mm(i)       = 55.079                !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.50958981615666
        IF (pol >= 1) parame(i,2) = 3.39806320429568
        IF (pol >= 1) parame(i,3) = 239.152759066148
        IF (pol >= 1) parame(i,6) = 4.05000000000000
        IF (pol >= 2) mm(i)       = 55.079                !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.54684827683436
        IF (pol >= 2) parame(i,2) = 3.41240089912190
        IF (pol >= 2) parame(i,3) = 218.299491580335
        IF (pol >= 2) parame(i,6) = 4.05000000000000
        IF (pol == 2) parame(i,11)= 6.24000000000000
        ! IF (pol.GE.2) mm(i)       = 55.079                !  PCIP-SAFT my_DD adjusted
        ! IF (pol.GE.2) parame(i,1) = 2.61175151088002
        ! IF (pol.GE.2) parame(i,2) = 3.37194293181453
        ! IF (pol.GE.2) parame(i,3) = 233.346110749402
        ! IF (pol.GE.2) parame(i,6) = 3.74682245835235
        ! IF (pol.EQ.2) parame(i,11)= 6.24000000000000
     CASE ('nitromethane')
        mm(i)       = 61.04                 !  PC-SAFT
        parame(i,1) = 2.58429167547409
        parame(i,2) = 3.10839592337018
        parame(i,3) = 310.694151426943
        IF (pol >= 1) mm(i)       = 61.04                 !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.55847635262615
        IF (pol >= 1) parame(i,2) = 3.10129282495975
        IF (pol >= 1) parame(i,3) = 256.456941430554
        IF (pol >= 1) parame(i,6) = 3.46000000000000
        IF (pol >= 2) mm(i)       = 61.04                 !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.68229497771588
        IF (pol >= 2) parame(i,2) = 3.10654492320028
        IF (pol >= 2) parame(i,3) = 225.973607468282
        IF (pol >= 2) parame(i,6) = 3.46000000000000
        IF (pol >= 2) parame(i,11)= 7.37000000000000
     CASE ('nitroethane')
        mm(i)       = 75.067                !  PC-SAFT
        parame(i,1) = 3.01767629617259
        parame(i,2) = 3.21364231060938
        parame(i,3) = 286.571650044235
        IF (pol >= 1) mm(i)       = 75.067                !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.94901220582233
        IF (pol >= 1) parame(i,2) = 3.23117331990738
        IF (pol >= 1) parame(i,3) = 265.961000131109
        IF (pol >= 1) parame(i,6) = 3.23000000000000
        IF (pol >= 2) mm(i)       = 75.067                !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 3.09101689452968
        IF (pol >= 2) parame(i,2) = 3.19364569858756
        IF (pol >= 2) parame(i,3) = 246.676040248662
        IF (pol >= 2) parame(i,6) = 3.23000000000000
        IF (pol >= 2) parame(i,11)= 9.63000000000000
     CASE ('acetonitrile')
        mm(i)       = 41.052                !  PC-SAFT
        parame(i,1) = 2.32895689571957
        parame(i,2) = 3.18980108373791
        parame(i,3) = 311.307486044181
        IF (pol >= 1) mm(i)       = 41.052                !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.15721401484941
        IF (pol >= 1) parame(i,2) = 3.27301469369132
        IF (pol >= 1) parame(i,3) = 216.888948676921
        IF (pol >= 1) parame(i,6) = 3.92520000000000
        IF (pol >= 2) mm(i)       = 41.052                !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.10426253849664
        IF (pol >= 2) parame(i,2) = 3.39403305120647
        IF (pol >= 2) parame(i,3) = 199.070191065791
        IF (pol >= 2) parame(i,6) = 3.92520000000000
        IF (pol >= 2) parame(i,11)= 4.40000000000000
        ! IF (pol >= 2) mm(i)       = 41.052                !  PCIP-SAFT my_DD adjusted
        ! IF (pol >= 2) parame(i,1) = 2.36268539768398
        ! IF (pol >= 2) parame(i,2) = 3.18554306395900
        ! IF (pol >= 2) parame(i,3) = 225.143934506015
        ! IF (pol >= 2) parame(i,6) = 3.43151866932598
        ! IF (pol >= 2) parame(i,11)= 4.40000000000000
        ! mm(i)       = 41.053                   ! PCP-SAFT dipole and quadrupole
        ! parame(i,1) = 1.79993
        ! parame(i,2) = 3.47366
        ! parame(i,3) = 211.001
        ! parame(i,6) = 3.93800
        ! parame(i,7) = 2.44000
        ! parame(i,11)= 0.0
     CASE ('dmf')
        mm(i)       = 73.09              !  PC-SAFT
        parame(i,1) = 2.388
        parame(i,2) = 3.658
        parame(i,3) = 363.77
        IF (pol >= 1) mm(i)       = 73.09              !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.269
        IF (pol >= 1) parame(i,2) = 3.714
        IF (pol >= 1) parame(i,3) = 331.56
        IF (pol >= 1) parame(i,6) = 3.82
        IF (pol >= 2) mm(i)       = 73.09              !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.375
        IF (pol >= 2) parame(i,2) = 3.667
        IF (pol >= 2) parame(i,3) = 308.42
        IF (pol >= 2) parame(i,6) = 3.82
        IF (pol >= 2) parame(i,11)= 7.81
     CASE ('chloroform')
        mm(i)       = 119.377              !  PCIP-SAFT
        parame(i,1) = 2.5957
        parame(i,2) = 3.4299
        parame(i,3) = 264.664
        parame(i,6) = 1.04
        IF (pol == 2) parame(i,11)= 8.23
     CASE ('dimethyl-ether')
        mm(i)       = 46.069              ! PC-SAFT
        parame(i,1) = 2.26234331
        parame(i,2) = 3.27664053
        parame(i,3) = 212.934324
        IF (pol >= 1) mm(i)       = 46.06900000     ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.219164566
        IF (pol >= 1) parame(i,2) = 3.296939638
        IF (pol >= 1) parame(i,3) = 212.1048888
        IF (pol >= 1) parame(i,6) = 1.300000000
        IF (pol >= 2) mm(i)       = 46.06900000     ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.275432546
        IF (pol >= 2) parame(i,2) = 3.265847188
        IF (pol >= 2) parame(i,3) = 206.9045519
        IF (pol >= 2) parame(i,6) = 1.300000000
        IF (pol == 2) parame(i,11)= 5.290000000
     CASE ('methyl-ethyl-ether')
        mm(i)     =           60.096
        parame(i,1) = mm(i)*  .0442404671
        parame(i,2) =           3.37282595
        parame(i,3) =           216.010217
        IF (pol >= 1) mm(i)       = 60.096            ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.6425218483
        IF (pol >= 1) parame(i,2) = 3.3793846539
        IF (pol >= 1) parame(i,3) = 215.78717386
        IF (pol >= 1) parame(i,6) = 1.1700000000
        IF (pol >= 2) mm(i)       = 60.096            ! PICP-SAFT
        IF (pol >= 2) parame(i,1) = 2.6790914671
        IF (pol >= 2) parame(i,2) = 3.3610534228
        IF (pol >= 2) parame(i,3) = 212.87191199
        IF (pol >= 2) parame(i,6) = 1.1700000000
        IF (pol >= 2) parame(i,11)= 7.9300000000
     CASE ('diethyl-ether')
        mm(i)       = 74.123                          ! PC-SAFT
        parame(i,1) = 3.0368496
        parame(i,2) = 3.4856955
        parame(i,3) = 217.64113
        IF (pol >= 1) mm(i)       = 74.123            ! PCP-SAFT as published 2006
        IF (pol >= 1) parame(i,1) = 2.97256367
        IF (pol >= 1) parame(i,2) = 3.5126868769
        IF (pol >= 1) parame(i,3) = 219.52737657
        IF (pol >= 1) parame(i,6) = 1.1500000000
        IF (pol >= 1) mm(i)       = 74.123            ! PCP-SAFT after using data from DDB
        IF (pol >= 1) parame(i,1) = 2.9528473413
        IF (pol >= 1) parame(i,2) = 3.5205207103
        IF (pol >= 1) parame(i,3) = 220.29852329
        IF (pol >= 1) parame(i,6) = 1.1500000000
        IF (pol >= 2) mm(i)       = 74.123            ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.9956379
        IF (pol >= 2) parame(i,2) = 3.501724569
        IF (pol >= 2) parame(i,3) = 217.8941822
        IF (pol >= 2) parame(i,6) = 1.15
        IF (pol == 2) parame(i,11)= 8.73
     CASE ('dipropyl-ether')
        mm(i)       = 102.17            ! PC-SAFT
        parame(i,1) = 3.600
        parame(i,2) = 3.662
        parame(i,3) = 230.38
        IF (pol >= 1) mm(i)       = 102.17            ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = 3.597
        IF (pol >= 1) parame(i,2) = 3.663
        IF (pol >= 1) parame(i,3) = 230.32
        IF (pol >= 1) parame(i,6) = 1.21
        IF (pol >= 2) mm(i)       = 102.17            ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 3.613
        IF (pol >= 2) parame(i,2) = 3.657
        IF (pol >= 2) parame(i,3) = 229.32
        IF (pol >= 2) parame(i,6) = 1.21
        IF (pol == 2) parame(i,11)= 12.8
     CASE ('vinylacetate')
        mm(i)       = 86.092
        parame(i,1) = mm(i)*  .0374329292
        parame(i,2) = 3.35278602
        parame(i,3) = 240.492049
     CASE ('chloromethane')    ! R40
        mm(i)       = 50.488        ! PC-SAFT
        parame(i,1) = mm(i)*  0.039418879  !  1.9902
        parame(i,2) = 3.1974
        parame(i,3) = 237.27
        IF (pol >= 1) mm(i)       = 50.488        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  0.035790801  !  1.8070
        IF (pol >= 1) parame(i,2) = 3.3034
        IF (pol >= 1) parame(i,3) = 229.97
        IF (pol >= 1) parame(i,6) = 1.8963
        IF (pol >= 1) lli(i)      = 1.67703*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 20.75417
        IF (pol >= 1) chap(i)     = 0.5
        IF (pol >= 2) mm(i)       = 50.488        ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  3.68559992E-2  !  1.86078
        IF (pol >= 2) parame(i,2) = 3.275186
        IF (pol >= 2) parame(i,3) = 216.4621
        IF (pol >= 2) parame(i,6) = 1.8963
        IF (pol == 2) parame(i,11)= 4.72
     CASE ('fluoromethane')    ! R41
        IF (pol /= 1) write (*,*) 'non-polar parameters missing for fluoromethane'
        IF (pol /= 1) stop 5
        IF (pol >= 1) mm(i)       =   34.0329000000000
        IF (pol >= 1) parame(i,1) =   1.94494757526896
        IF (pol >= 1) parame(i,2) =   2.96858005012635
        IF (pol >= 1) parame(i,3) =   168.938697391009
        IF (pol >= 1) parame(i,6) =         1.57823038894029
     CASE ('dichloromethane')   ! R30
        mm(i)       = 84.932        ! PC-SAFT
        parame(i,1) = 2.3117
        parame(i,2) = 3.3161
        parame(i,3) = 270.98
        IF (pol >= 1) mm(i)       = 84.932        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = 2.2687
        IF (pol >= 1) parame(i,2) = 3.3373
        IF (pol >= 1) parame(i,3) = 269.08
        IF (pol >= 1) parame(i,6) = 1.6
        IF (pol >= 2) mm(i)       = 84.932        ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = 2.3435
        IF (pol >= 2) parame(i,2) = 3.2987
        IF (pol >= 2) parame(i,3) = 260.66
        IF (pol >= 2) parame(i,6) = 1.6
        IF (pol == 2) parame(i,11)= 6.48
     CASE ('difluoromethane')   ! R32
        IF (pol /= 1) write (*,*) 'non-polar parameters missing for difluoromethane'
        IF (pol /= 1) stop 5
        IF (pol >= 1) mm(i)       = 52.0236             ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  4.814700934384165E-002  !  2.50478075530028
        IF (pol >= 1) parame(i,2) = 2.79365980535456
        IF (pol >= 1) parame(i,3) = 160.893555378523
        IF (pol >= 1) parame(i,6) = 1.97850000000000
     CASE ('trifluoromethane')   ! R23
        IF (pol /= 1) write (*,*) 'non-polar parameters missing for trifluoromethane'
        IF (pol /= 1) stop 5
        IF (pol >= 1) mm(i)       =   70.0138000000000
        IF (pol >= 1) parame(i,1) =   2.66039274225485
        IF (pol >= 1) parame(i,2) =   2.82905884530501
        IF (pol >= 1) parame(i,3) =   149.527709542333
        IF (pol >= 1) parame(i,6) =   1.339963415253999E-002
     CASE ('tetrachloromethane')   ! R10
        mm(i)       = 153.822
        parame(i,1) = mm(i)*  .0150432213
        parame(i,2) = 3.81801454
        parame(i,3) = 292.838632
     CASE ('trichlorofluoromethane')   ! R11
        IF (pol /= 1) write (*,*) 'non-polar parameters missing for trichlorofluoromethane'
        IF (pol /= 1) stop 5
        IF (pol >= 1) mm(i)       =   137.368000000000
        IF (pol >= 1) parame(i,1) =   2.28793359008803
        IF (pol >= 1) parame(i,2) =   3.69013104930876
        IF (pol >= 1) parame(i,3) =   248.603173885090
        IF (pol >= 1) parame(i,6) =   0.23225538492979
     CASE ('chlorodifluoromethane')   ! R22   ( CHClF2 or CHF2Cl)
        IF (pol /= 1) write (*,*) 'non-polar parameters missing for chlorodifluoromethane'
        IF (pol /= 1) stop 5
        IF (pol >= 1) mm(i)       =   86.4684000000000
        IF (pol >= 1) parame(i,1) =   2.47218586047893
        IF (pol >= 1) parame(i,2) =   3.13845692489930
        IF (pol >= 1) parame(i,3) =   187.666355083434
        IF (pol >= 1) parame(i,6) =         1.04954264812860
     CASE ('chloroethane')
        mm(i)       = 64.514
        parame(i,1) = mm(i)*  .0350926868
        parame(i,2) = 3.41602397
        parame(i,3) = 245.42626
     CASE ('11difluoroethane')
        ! mm(i)       = 66.0500000000000         !  PC-SAFT
        ! parame(i,1) = mm(i)*  4.109944338817734E-002
        ! parame(i,2) = 3.10257444633546
        ! parame(i,3) = 192.177159144029
        ! mm(i)      = 66.05                             !  PC-SAFT assoc
        ! parame(i,1)= 2.984947188
        ! parame(i,2)= 2.978630027
        ! parame(i,3)= 137.8192282
        ! nhb_typ(i) = 2
        ! nhb_no(i,1)= 1
        ! nhb_no(i,2)= 1
        ! eps_hb(i,i,1,2)= 823.3478288
        ! eps_hb(i,i,2,1)= 823.3478288
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 0.96345994
        IF (pol >= 1) mm(i)       = 66.0500000000000        !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  3.949665745363346E-002  !  =2.60875422481249
        IF (pol >= 1) parame(i,2) = 3.13758353925036
        IF (pol >= 1) parame(i,3) = 179.517952627836
        IF (pol >= 1) parame(i,6) = 2.27000000000000
        IF (pol >= 1) lli(i)      = 2.03907*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 26.5
        IF (pol >= 1) chap(i)     = 0.4
        IF (pol >= 2) mm(i)       = 66.0500000000000        !  PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  4.093647666154238E-002  !  =2.70385428349487
        IF (pol >= 2) parame(i,2) = 3.10437129415885
        IF (pol >= 2) parame(i,3) = 170.464400902455
        IF (pol >= 2) parame(i,6) = 2.27000000000000
        IF (pol == 2) parame(i,11)= 5.01000000000000
     CASE ('1-chlorobutane')
        mm(i)       = 92.568
        parame(i,1) = mm(i)*  .0308793201
        parame(i,2) = 3.64240187
        parame(i,3) = 258.655298
     CASE ('chlorobenzene')
        ! mm(i)       = 112.558
        ! parame(i,1) = mm(i)*  .0235308686
        ! parame(i,2) = 3.75328494
        ! parame(i,3) = 315.039018
        mm(i)       = 112.558        ! PCIP-SAFT
        parame(i,1) = mm(i)*  0.023824167  !  =2.6816
        parame(i,2) = 3.7352
        parame(i,3) = 308.82
        parame(i,6) = 1.69
        IF (pol == 2) parame(i,11)= 14.1
     CASE ('styrene')
        mm(i)       = 104.150
        parame(i,1) = mm(i)*  2.9124104853E-2
        parame(i,2) = 3.760233548
        parame(i,3) = 298.51287564
     CASE ('methylmethanoate')
        mm(i)     = 60.053
        parame(i,1) = mm(i)*  .0446000264
        parame(i,2) = 3.08753499
        parame(i,3) = 242.626755
        IF (pol >= 1) mm(i)       = 60.053        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  4.366991153963102E-002  !  =2.62250919768946
        IF (pol >= 1) parame(i,2) = 3.10946396964
        IF (pol >= 1) parame(i,3) = 239.051951942
        IF (pol >= 1) parame(i,6) = 1.77
        IF (pol >= 2) mm(i)       = 60.053        ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  4.492572388931002E-2  !  2.69792449
        IF (pol >= 2) parame(i,2) = 3.078467837
        IF (pol >= 2) parame(i,3) = 232.1842551
        IF (pol >= 2) parame(i,6) = 1.77
        IF (pol == 2) parame(i,11)= 5.05
     CASE ('ethylmethanoate')
        mm(i)     =           74.079        ! PC-SAFT
        parame(i,1) = mm(i)*  .03898009
        parame(i,2) =           3.31087192
        parame(i,3) =           246.465646
        IF (pol >= 1) mm(i)       =         74.079        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  3.825407152074255E-002  !  =2.83382336418509
        IF (pol >= 1) parame(i,2) =         3.33160046679
        IF (pol >= 1) parame(i,3) =         244.495680932
        IF (pol >= 1) parame(i,6) =         1.93000000000
     CASE ('propylmethanoate')
        mm(i)     =           88.106
        parame(i,1) = mm(i)*  .0364206062
        parame(i,2) =           3.41679642
        parame(i,3) =           246.457732
        IF (pol >= 1) mm(i)       =         88.106
        IF (pol >= 1) parame(i,1) = mm(i)*  3.60050739149E-2  !  =3.17226304235232
        IF (pol >= 1) parame(i,2) =         3.42957609309
        IF (pol >= 1) parame(i,3) =         245.637644107
        IF (pol >= 1) parame(i,6) =         1.89
     CASE ('methylacetate')
        mm(i)       = 74.079        ! PC-SAFT
        parame(i,1) = mm(i)*  4.286817177E-2  !  =3.175631296
        parame(i,2) = 3.18722021277843
        parame(i,3) = 234.106931032456
        IF (pol >= 1) mm(i)       = 74.079        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  4.228922065E-2  !  =3.132743176
        IF (pol >= 1) parame(i,2) = 3.2011401688
        IF (pol >= 1) parame(i,3) = 233.17562886
        IF (pol >= 1) parame(i,6) = 1.72
        IF (pol >= 2) mm(i)       = 74.079        ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  4.298900538E-2  !  =3.18458252
        IF (pol >= 2) parame(i,2) = 3.180642322
        IF (pol >= 2) parame(i,3) = 229.3132680
        IF (pol >= 2) parame(i,6) = 1.72
        IF (pol == 2) parame(i,11)= 6.94
     CASE ('ethylacetate')
        mm(i)       = 88.106        ! PC-SAFT
        parame(i,1) = mm(i)*  .0401464427   !  =3.537142481
        parame(i,2) = 3.30789258
        parame(i,3) = 230.800689
        IF (pol >= 1) mm(i)       = 88.106        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  0.039792575   !  =3.505964572
        IF (pol >= 1) parame(i,2) = 3.317655188
        IF (pol >= 1) parame(i,3) = 230.2434769
        IF (pol >= 1) parame(i,6) = 1.78
        IF (pol >= 2) mm(i)       = 88.106        ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  0.040270267   !  =3.548052143
        IF (pol >= 2) parame(i,2) = 3.302097562
        IF (pol >= 2) parame(i,3) = 227.5026191
        IF (pol >= 2) parame(i,6) = 1.78
        IF (pol == 2) parame(i,11)= 8.62
     CASE ('ethyl-propanoate')
        mm(i)       = 102.133
        parame(i,1) = mm(i)*  .0375692464
        parame(i,2) = 3.40306071
        parame(i,3) = 232.778374
     CASE ('propyl-ethanoate')
        mm(i)       = 102.133
        parame(i,1) = mm(i)*  .0370721275
        parame(i,2) = 3.42272266
        parame(i,3) = 235.758378
        IF (pol >= 1) mm(i)       = 102.133        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  3.687149995200072E-2  !  =3.76579690459769
        IF (pol >= 1) parame(i,2) = 3.4289353421006
        IF (pol >= 1) parame(i,3) = 235.41679442910
        IF (pol >= 1) parame(i,6) = 1.78
        ! IF (pol.EQ.2) parame(i,11)= 10.41
     CASE ('nbutyl-ethanoate')
        mm(i)       = 116.16        ! PC-SAFT
        parame(i,1) = mm(i)*  .03427004
        parame(i,2) = 3.54269638
        parame(i,3) = 242.515768
        IF (pol >= 1) mm(i)       = 116.16        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  3.411585209773470E-002  !  =3.96289737967286
        IF (pol >= 1) parame(i,2) = 3.54821589228130
        IF (pol >= 1) parame(i,3) = 242.274388267447
        IF (pol >= 1) parame(i,6) = 1.87000000000000
        IF (pol >= 2) mm(i)       = 116.16        ! PCIP-SAFT
        IF (pol >= 2) parame(i,1) = mm(i)*  3.442139015733717E-002  !  =3.99838868067629
        IF (pol >= 2) parame(i,2) = 3.53576054452119
        IF (pol >= 2) parame(i,3) = 240.154409609249
        IF (pol >= 2) parame(i,6) = 1.87000000000000
        IF (pol == 2) parame(i,11)= 14.2000000000000
     CASE ('methyl-octanoate')
        mm(i)       = 158.24        ! PC-SAFT
        parame(i,1) = 5.2074
        parame(i,2) = 3.6069
        parame(i,3) = 244.12
     CASE ('methyl-decanoate')
        mm(i)       = 186.2912        ! PC-SAFT
        parame(i,1) = 5.8402
        parame(i,2) = 3.6871
        parame(i,3) = 248.27

        mm(i)       = 186.2912        ! PC-SAFT from GC-method Tim
        parame(i,1) = 7.716
        parame(i,2) = 3.337303029
        parame(i,3) = 204.250907

        mm(i)       = 186.2912        ! PC-SAFT from GC-method (tightly fit) Tim
        parame(i,1) = 7.728
        parame(i,2) = 3.334023322
        parame(i,3) = 206.9099379

        ! mm(i)       = 186.2912        ! PC-SAFT from fit to DIPPR
        ! parame(i,1) = 6.285005
        ! parame(i,2) = 3.594888
        ! parame(i,3) = 236.781461
        ! ! parame(i,6) = 2.08056

        ! mm(i)       =         186.291000000000
        ! parame(i,1) =    6.28500485898895
        ! parame(i,2) =         3.59488828061149
        ! parame(i,3) =         236.781461491921
        ! parame(i,6) =         2.08055996894836
        ! parame(i,8) =         1.00000000000000
        mm(i)       =         186.291000000000
        parame(i,1) =    6.14436331493372
        parame(i,2) =         3.61977264863944
        parame(i,3) =         242.071887817656

     CASE ('methyl-dodecanoate')
        mm(i)       = 214.344        ! PC-SAFT
        parame(i,1) = 6.5153
        parame(i,2) = 3.7406
        parame(i,3) = 250.7
     CASE ('methyl-tetradecanoate')
        mm(i)       = 242.398        ! PC-SAFT
        parame(i,1) = 7.1197
        parame(i,2) = 3.7968
        parame(i,3) = 253.77
     CASE ('methyl-hexadecanoate')
        mm(i)       = 270.451        ! PC-SAFT
        parame(i,1) = 7.891
        parame(i,2) = 3.814
        parame(i,3) = 253.71
     CASE ('methyl-octadecanoate')
        mm(i)       = 298.504        ! PC-SAFT
        parame(i,1) = 8.8759
        parame(i,2) = 3.7932
        parame(i,3) = 250.81
     CASE ('CH2F2')
        mm(i)       = 52.02
        parame(i,1) = 3.110084171
        parame(i,2) = 2.8145230485
        parame(i,3) = 158.98060151
     CASE ('naphthalene')
        ! mm(i)       = 128.174000000
        ! parame(i,1) = mm(i)*  2.4877834216412E-2
        ! parame(i,2) = 3.82355815011
        ! parame(i,3) = 341.560675334

        mm(i)       = 128.17400000000
        parame(i,1) = mm(i)*  2.6400924157729E-2
        parame(i,2) = 3.8102186020014
        parame(i,3) = 328.96792935903
     CASE ('h2s')
        mm(i)       =         34.0820000000000     ! PC-SAFT
        parame(i,1) = mm(i)*  4.838886696385162E-002   ! =    1.64918936386199
        parame(i,2) =         3.05478289838459
        parame(i,3) =         229.838873939562
        nhb_typ(i)     =              2
        nhb_no(i,1)    =              1
        nhb_no(i,2)    =              1
        eps_hb(i,i,1,2)=      536.634834731413
        eps_hb(i,i,2,1)=      536.634834731413
        eps_hb(i,i,1,1)=  0.0
        eps_hb(i,i,2,2)=  0.0
        kap_hb(i,i)    =     1.000000000000000E-003
        ! PC-SAFT from Xiaohua
        mm(i)       =         34.082     ! PC-SAFT
        parame(i,1) = 1.63677
        parame(i,2) =         3.06565
        parame(i,3) =         230.2121
        nhb_typ(i)     =              2
        nhb_no(i,1)    =              1
        nhb_no(i,2)    =              1
        eps_hb(i,i,1,2)=      275.1088
        eps_hb(i,i,2,1)=      275.1088
        eps_hb(i,i,1,1)=  0.0
        eps_hb(i,i,2,2)=  0.0
        kap_hb(i,i)    =     1.E-2
        ! IF (pol.GE.1) mm(i)       =         34.082       ! PCP-SAFT with quadrupole
        ! IF (pol.GE.1) parame(i,1) = mm(i)*  3.03171032558E-2  !  =1.03326751316478
        ! IF (pol.GE.1) parame(i,2) =         3.6868189914018
        ! IF (pol.GE.1) parame(i,3) =         246.862831266172
        ! IF (pol.GE.1) nhb_typ(i)     =              2
        ! IF (pol.GE.1) nhb_no(i,1)    =              1
        ! IF (pol.GE.1) nhb_no(i,2)    =              1
        ! IF (pol.GE.1) eps_hb(i,i,1,2)=      987.4927232
        ! IF (pol.GE.1) eps_hb(i,i,2,1)=      987.4927232
        ! IF (pol.GE.1) eps_hb(i,i,1,1)=  0.0
        ! IF (pol.GE.1) eps_hb(i,i,2,2)=  0.0
        ! IF (pol.GE.1) kap_hb(i,i)    =     5.5480659623d-4
        ! IF (pol.GE.1) parame(i,6) =        0.97833
        ! IF (pol.GE.1) parame(i,7) =         3.8623
        ! IF (pol.GE.1) LLi(i)      = 1.2737*parame(i,2)
        ! IF (pol.GE.1) phi_criti(i)= 14.316
        ! IF (pol.GE.1) chap(i)     = 0.4473
        IF (pol >= 1) mm(i)       =         34.0820000000000       ! PCP-SAFT no quadrupoLE
        IF (pol >= 1) parame(i,1) = mm(i)*  4.646468487062725E-002 ! 1.58360938976072
        IF (pol >= 1) parame(i,2) =         3.10111012646306
        IF (pol >= 1) parame(i,3) =         230.243457544889
        IF (pol >= 1) nhb_typ(i)     =              2
        IF (pol >= 1) nhb_no(i,1)    =              1
        IF (pol >= 1) nhb_no(i,2)    =              1
        IF (pol >= 1) eps_hb(i,i,1,2)=      584.367708701411
        IF (pol >= 1) eps_hb(i,i,2,1)=      584.367708701411
        IF (pol >= 1) eps_hb(i,i,1,1)=  0.0
        IF (pol >= 1) eps_hb(i,i,2,2)=  0.0
        IF (pol >= 1) kap_hb(i,i)    =     1.000000000000000E-003
        IF (pol >= 1) parame(i,6) =        0.978330000000000

        IF (pol >= 1) lli(i)      = 1.2737*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 14.316
        IF (pol >= 1) chap(i)     = 0.4473


        IF (pol == 2) parame(i,7) =        0.0
        IF (pol == 2) mm(i)       =         34.0820000000000     ! PCIP-SAFT
        IF (pol == 2) parame(i,1) = mm(i)*  4.806418212963168E-002 ! 1.63812345534211
        IF (pol == 2) parame(i,2) =         3.06556006883749
        IF (pol == 2) parame(i,3) =         221.746622243054
        IF (pol == 2) nhb_typ(i)     =              2
        IF (pol == 2) nhb_no(i,1)    =              1
        IF (pol == 2) nhb_no(i,2)    =              1
        IF (pol == 2) eps_hb(i,i,1,2)=      672.164783984789
        IF (pol == 2) eps_hb(i,i,2,1)=      672.164783984789
        IF (pol == 2) eps_hb(i,i,1,1)=  0.0
        IF (pol == 2) eps_hb(i,i,2,2)=  0.0
        IF (pol == 2) kap_hb(i,i)    =     1.000000000000000E-003
        IF (pol == 2) parame(i,6) =        0.978330000000000
        IF (pol == 2) parame(i,11) =         3.60200000000000
        IF (pol == 2) parame(i,7) =        0.0

        IF (pol >= 1)mm(i)       =    34.0820000000000     !PCP-SAFT D&Q
        IF (pol >= 1)parame(i,1) = mm(i)*  3.974667896078737E-002  !  = 1.35464631234155
        IF (pol >= 1)parame(i,2) =         3.30857082333438
        IF (pol >= 1)parame(i,3) =         234.248947273191
        IF (pol >= 1)nhb_typ(i)     =              2
        IF (pol >= 1)nhb_no(i,1)    =              1
        IF (pol >= 1)nhb_no(i,2)    =              1
        IF (pol >= 1)eps_hb(i,i,1,2)=      780.770936834770
        IF (pol >= 1)eps_hb(i,i,2,1)=      780.770936834770
        IF (pol >= 1)eps_hb(i,i,1,1)=  0.0
        IF (pol >= 1)eps_hb(i,i,2,2)=  0.0
        IF (pol >= 1)kap_hb(i,i)    =     1.000000000000000E-003
        IF (pol >= 1)parame(i,6) =        0.978330000000000
        IF (pol >= 1)parame(i,7) =         2.93750500000000

     CASE ('methanol')
        mm(i)       = 32.042        ! PC-SAFT
        parame(i,1) = mm(i)*  .0476100379
        parame(i,2) = 3.23000005
        parame(i,3) = 188.904644
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2899.49055
        eps_hb(i,i,2,1)= 2899.49055
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0351760892
        IF (pol >= 1) mm(i)       = 32.042       ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  7.213091821E-2  !  =2.31121888139672
        IF (pol >= 1) parame(i,2) = 2.8270129502
        IF (pol >= 1) parame(i,3) = 176.3760515
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2332.5845803
        IF (pol >= 1) eps_hb(i,i,2,1)= 2332.5845803
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 8.9248658086E-2
        IF (pol >= 1) parame(i,6) = 1.7
        IF (pol >= 1) lli(i)      = 1.75*parame(i,2)
        IF (pol >= 1) phi_criti(i)= 23.43
        IF (pol >= 1) chap(i)     = 0.304
        IF (pol == 2) mm(i)       = 32.042       ! PCIP-SAFT
        IF (pol == 2) parame(i,1) = 2.0693
        IF (pol == 2) parame(i,2) = 2.9547
        IF (pol == 2) parame(i,3) = 174.51
        IF (pol == 2) nhb_typ(i)     = 2
        IF (pol == 2) nhb_no(i,1)    = 1
        IF (pol == 2) nhb_no(i,2)    = 1
        IF (pol == 2) eps_hb(i,i,1,2)= 2418.5
        IF (pol == 2) eps_hb(i,i,2,1)= 2418.5
        IF (pol == 2) eps_hb(i,i,1,1)= 0.0
        IF (pol == 2) eps_hb(i,i,2,2)= 0.0
        IF (pol == 2) kap_hb(i,i)    = 0.06319
        IF (pol == 2) parame(i,6) = 1.7
        IF (pol == 2) parame(i,11)= 3.29
        ! mm(i)       =         32.0420000000000     ! PCP-SAFT with adjusted QQ
        ! parame(i,1) = mm(i)*  6.241807629559099E-002
        ! ! parame(i,1) =    2.00000000066333
        ! parame(i,2) =         2.97610169698593
        ! parame(i,3) =         163.268505098639
        ! nhb_typ(i)     =              2
        ! nhb_no(i,1)    =              1
        ! nhb_no(i,2)    =              1
        ! eps_hb(i,i,1,2)=      2449.55621933612
        ! eps_hb(i,i,2,1)=      2449.55621933612
        ! eps_hb(i,i,1,1)=  0.0
        ! eps_hb(i,i,2,2)=  0.0
        ! kap_hb(i,i)    =     8.431015160393653E-002
        ! parame(i,6) =         1.72000000000000
        ! parame(i,7) =         1.59810028824523

        mm(i)       = 32.042        ! PC-SAFT
        parame(i,1) = 1.507
        parame(i,2) = 3.325
        parame(i,3) = 211.6
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2519.7
        eps_hb(i,i,2,1)= 2519.7
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 0.03
        parame(i,6) = 1.7
     CASE ('ethanol')
        mm(i)       = 46.069
        parame(i,1) = mm(i)*  .0517195521
        parame(i,2) = 3.17705595
        parame(i,3) = 198.236542
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2653.38367
        eps_hb(i,i,2,1)= 2653.38367
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0323840159
        IF (pol >= 1) mm(i)       = 46.0690000000000
        IF (pol >= 1) parame(i,1) = mm(i)*  4.753626908781145E-002  !  =2.18994838060639
        IF (pol >= 1) parame(i,2) = 3.30120000000000
        IF (pol >= 1) parame(i,3) = 209.824555801706
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2584.53116785767
        IF (pol >= 1) eps_hb(i,i,2,1)= 2584.53116785767
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 2.349382956935725E-002
        IF (pol >= 1) parame(i,6) = 1.69000000000000
        ! mm(i)       = 46.0690000000000
        ! parame(i,1) = mm(i)*  5.117957752785066E-002   !  =2.357791957
        ! parame(i,2) = 3.24027031244304
        ! parame(i,3) = 175.657110615456
        ! nhb_typ(i)     = 2
        ! nhb_no(i,1)    = 1
        ! nhb_no(i,2)    = 1
        ! eps_hb(i,i,1,2)= 2273.62670516146
        ! eps_hb(i,i,2,1)= 2273.62670516146
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 7.030279197039477E-002
        ! parame(i,6) = 1.69000000000000
        ! parame(i,7) = 3.63701294195013

        ! IF (pol >= 1) mm(i)       = 46.068   !  von Elmar
        ! IF (pol >= 1) parame(i,1) = 2.372
        ! IF (pol >= 1) parame(i,2) = 3.196
        ! IF (pol >= 1) parame(i,3) = 203.8
        ! IF (pol >= 1) nhb_typ(i)     = 2
        ! IF (pol >= 1) nhb_no(i,1)    = 1
        ! IF (pol >= 1) nhb_no(i,2)    = 1
        ! IF (pol >= 1) eps_hb(i,i,1,2)= 2514.1
        ! IF (pol >= 1) eps_hb(i,i,2,1)= 2514.1
        ! IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        ! IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        ! IF (pol >= 1) kap_hb(i,i)    = 0.030
        ! IF (pol >= 1) parame(i,6) = 1.691

        IF (pol == 2) mm(i)       = 46.0690000000000
        IF (pol == 2) parame(i,1) = mm(i)*  4.733436280008321E-002  !  =2.18064676
        IF (pol == 2) parame(i,2) = 3.31260000000000
        IF (pol == 2) parame(i,3) = 207.594119926613
        IF (pol == 2) nhb_typ(i)  = 2
        IF (pol == 2) nhb_no(i,1) = 1
        IF (pol == 2) nhb_no(i,2) = 1
        IF (pol == 2) eps_hb(i,i,1,2)= 2589.68311382661
        IF (pol == 2) eps_hb(i,i,2,1)= 2589.68311382661
        IF (pol == 2) eps_hb(i,i,1,1)= 0.0
        IF (pol == 2) eps_hb(i,i,2,2)= 0.0
        IF (pol == 2) kap_hb(i,i) = 2.132561218631547E-002
        IF (pol == 2) parame(i,6) = 1.69000000000000
        IF (pol == 2) parame(i,7) = 0.0
        IF (pol == 2) parame(i,11)= 5.11000000000000
     CASE ('1-propanol')
        mm(i)       = 60.096
        parame(i,1) = mm(i)*  .0499154461
        parame(i,2) = 3.25221234
        parame(i,3) = 233.396705
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2276.77867
        eps_hb(i,i,2,1)= 2276.77867
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0152683094
     CASE ('1-butanol')
        mm(i)       = 74.123
        parame(i,1) = mm(i)*  .0341065046
        parame(i,2) = 3.72361538
        parame(i,3) = 269.798048
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2661.37119
        eps_hb(i,i,2,1)= 2661.37119
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .00489087833
        mm(i)       = 74.1230000000000
        parame(i,1) = mm(i)*  3.329202420547412E-002  !  =2.46770471018236
        parame(i,2) = 3.76179376417092
        parame(i,3) = 270.237284242002
        nhb_typ(i)     = 2
        nhb_no(i,1)    = 1
        nhb_no(i,2)    = 1
        eps_hb(i,i,1,2)= 2669.28754983370
        eps_hb(i,i,2,1)= 2669.28754983370
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 4.855584122733399E-003
        parame(i,6) = 1.66000000000000
     CASE ('1-pentanol')
        mm(i)       = 88.15               ! PC-SAFT
        parame(i,1) = mm(i)*  .041134139
        parame(i,2) = 3.45079143
        parame(i,3) = 247.278748
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2252.09237
        eps_hb(i,i,2,1)= 2252.09237
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0103189939
        IF (pol >= 1) mm(i)       = 88.1500000000000        ! PCP-SAFT
        IF (pol >= 1) parame(i,1) = mm(i)*  4.138903382168521E-002  !  =3.64844333138155
        IF (pol >= 1) parame(i,2) = 3.44250118689142
        IF (pol >= 1) parame(i,3) = 246.078034174947
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2236.72830142446
        IF (pol >= 1) eps_hb(i,i,2,1)= 2236.72830142446
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 1.040067895187016E-002
        IF (pol >= 1) parame(i,6) = 1.70000000000000
        IF (pol == 2) mm(i)       = 88.1500000000000        ! PCIP-SAFT
        IF (pol == 2) parame(i,1) = mm(i)*  4.161521814399406E-002  !  =3.66838147939308
        IF (pol == 2) parame(i,2) = 3.43496654431777
        IF (pol == 2) parame(i,3) = 244.177313808431
        IF (pol == 2) nhb_typ(i)     = 2
        IF (pol == 2) nhb_no(i,1)    = 1
        IF (pol == 2) nhb_no(i,2)    = 1
        IF (pol == 2) eps_hb(i,i,1,2)= 2241.27880639096
        IF (pol == 2) eps_hb(i,i,2,1)= 2241.27880639096
        IF (pol == 2) eps_hb(i,i,1,1)= 0.0
        IF (pol == 2) eps_hb(i,i,2,2)= 0.0
        IF (pol == 2) kap_hb(i,i)    = 1.049516309928397E-002
        IF (pol == 2) parame(i,6) = 1.70000000000000
        IF (pol == 2) parame(i,11)= 10.8000000000000
     CASE ('1-octanol')
        mm(i)       = 130.23
        parame(i,1) = mm(i)*  .0334446084
        parame(i,2) = 3.714535
        parame(i,3) = 262.740637
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2754.77272
        eps_hb(i,i,2,1)= 2754.77272
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .00219656803
     CASE ('1-nonanol')
        mm(i)       = 144.257
        parame(i,1) = mm(i)*  .0324692669
        parame(i,2) = 3.72924286
        parame(i,3) = 263.636673
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2941.9231
        eps_hb(i,i,2,1)= 2941.9231
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .00142696883
     CASE ('2-propanol')
        mm(i)       = 60.096
        parame(i,1) = mm(i)*  .0514663133
        parame(i,2) = 3.20845858
        parame(i,3) = 208.420809
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2253.91064
        eps_hb(i,i,2,1)= 2253.91064
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0246746934
        IF (pol >= 1) mm(i)       = 60.096           ! PCP-SAFT (Parameter von Elmar)
        IF (pol >= 1) parame(i,1) = 4.025
        IF (pol >= 1) parame(i,2) = 2.918
        IF (pol >= 1) parame(i,3) = 198.6
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 1871.9
        IF (pol >= 1) eps_hb(i,i,2,1)= 1871.9
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 0.030
        IF (pol >= 1) parame(i,6) = 1.661
     CASE ('2-methyl-2-butanol')
        mm(i)       = 88.15
        parame(i,1) = mm(i)*  .0289135026
        parame(i,2) = 3.90526707
        parame(i,3) = 266.011828
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2618.80124
        eps_hb(i,i,2,1)= 2618.80124
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .00186263367
     CASE ('ethanediol')
        mm(i)       = 62.07
        parame(i,1) = 3.64151628
        parame(i,2) = 2.78381944
        parame(i,3) = 262.53650804
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2476.87885
        eps_hb(i,i,2,1)= 2476.87885
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 0.1259948912
        IF (pol >= 1) parame(i,1) = 3.57454606 !  PCP-SAFT
        IF (pol >= 1) parame(i,2) = 2.80595405
        IF (pol >= 1) parame(i,3) = 255.64383587
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2439.633888
        IF (pol >= 1) eps_hb(i,i,2,1)= 2439.633888
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 0.1294968485
        IF (pol >= 1) parame(i,6) = 2.4103116    ! dipole moment
     CASE ('glycerol')
        mm(i)       = 92.09
        parame(i,1) = 2.09257199
        parame(i,2) = 3.79298757
        parame(i,3) = 492.74380031
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 3746.684679
        eps_hb(i,i,2,1)= 3746.684679
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 0.0014484419
        IF (pol >= 1) parame(i,1) = 1.95438379 !  PCP-SAFT
        IF (pol >= 1) parame(i,2) = 3.89445419
        IF (pol >= 1) parame(i,3) = 502.20994609
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 3852.103505
        IF (pol >= 1) eps_hb(i,i,2,1)= 3852.103505
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 0.0009570965
        IF (pol >= 1) parame(i,6) = 2.6801226
     CASE ('monoethanolamine')
        mm(i)       = 61.083
        IF (pol >= 1) parame(i,1) = 3.342
        IF (pol >= 1) parame(i,2) = 2.955
        IF (pol >= 1) parame(i,3) = 289.2
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2117.4
        IF (pol >= 1) eps_hb(i,i,2,1)= 2117.4
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 0.030
        IF (pol >= 1) parame(i,6) = 0.776
     CASE ('formamide')
        mm(i)       = 45.04
        parame(i,1) = 2.30026276
        parame(i,2) = 2.92582027
        parame(i,3) = 313.56210527
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2679.789451
        eps_hb(i,i,2,1)= 2679.789451
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 0.2373750869
        IF (pol >= 1) parame(i,1) = 1.73837559 !  PCP-SAFT
        IF (pol >= 1) parame(i,2) = 3.32497645
        IF (pol >= 1) parame(i,3) = 362.80434248
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2063.301444
        IF (pol >= 1) eps_hb(i,i,2,1)= 2063.301444
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 0.0192400111
        IF (pol >= 1) parame(i,6) = 3.717396
     CASE ('acetic-acid')
        mm(i)       = 60.053
        parame(i,1) = 1.36366520
        parame(i,2) = 3.79651163
        parame(i,3) = 199.225066
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 3092.40109
        eps_hb(i,i,2,1)= 3092.40109
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0870093874


        mm(i)       = 60.053
        parame(i,1) = 1.09174940
        parame(i,2) = 4.13711044
        parame(i,3) = 207.552969
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 3198.84362
        eps_hb(i,i,2,1)= 3198.84362
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0586552968

        ! mit gesetztem Dipol-Moment
!!$        mm(i)       = 60.053
!!$        parame(i,1) = 1.0427724
!!$        parame(i,2) = 4.2522071
!!$        parame(i,3) = 190.95725
!!$        nhb_typ(i)  = 2
!!$        nhb_no(i,1) = 1            ! no. of sites of type 1
!!$        nhb_no(i,2) = 1            ! no. of sites of type 2
!!$        eps_hb(i,i,1,2)= 3096.3619
!!$        eps_hb(i,i,2,1)= 3096.3619
!!$        eps_hb(i,i,1,1)= 0.0
!!$        eps_hb(i,i,2,2)= 0.0
!!$        kap_hb(i,i)    = 6.154307E-002
!!$        parame(i,6) = 3.50
     CASE ('propionic-acid')
        mm(i)       = 74.0800000000000
        parame(i,1) = mm(i)*  2.359519915877884E-002
        parame(i,2) = 3.99339224153844
        parame(i,3) = 295.947729838532
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2668.97826430874
        eps_hb(i,i,2,1)= 2668.97826430874
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 3.660242292423115E-002
     CASE ('acrylic-acid')
        mm(i)       = 72.0636
        parame(i,1) = mm(i)*  .0430585424
        parame(i,2) = 3.0545415
        parame(i,3) = 164.115604
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 3065.40667
        eps_hb(i,i,2,1)= 3065.40667
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .336261669
     CASE ('caproic-acid')
        mm(i)       =         116.16
        parame(i,1) = 5.87151
        parame(i,2) =         3.0694697
        parame(i,3) =         241.4569
        nhb_typ(i)     =      1
        eps_hb(i,i,1,1)=      2871.37
        kap_hb(i,i)    =      3.411613D-3
     CASE ('aniline')
        mm(i)       = 93.13
        parame(i,1) = mm(i)*  .0285695992
        parame(i,2) = 3.70214085
        parame(i,3) = 335.471062
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 1351.64234
        eps_hb(i,i,2,1)= 1351.64234
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0748830615

        mm(i)       = 93.1300000000000
        parame(i,1) = mm(i)*  2.834372610192228E-002  !  =2.63965121187202
        parame(i,2) = 3.71326867619433
        parame(i,3) = 332.253796842009
        nhb_typ(i)     = 2
        nhb_no(i,1)    = 1
        nhb_no(i,2)    = 1
        eps_hb(i,i,1,2)= 1392.14266886674
        eps_hb(i,i,2,1)= 1392.14266886674
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 7.424612087328866E-002
        parame(i,6) = 1.55000000000000
        IF (pol == 2) parame(i,11)= 12.1000000000000

     CASE ('methylisocyanate')
        mm(i)       =   57.0540
        parame(i,1) =   2.14783454354850
        parame(i,2) =   3.30276435689525
        parame(i,3) =   284.359877866415
        IF (pol >= 1) parame(i,1) =   2.54196677366949
        IF (pol >= 1) parame(i,2) =   3.11423929858242
        IF (pol >= 1) parame(i,3) =   215.712884707899
        IF (pol >= 1) parame(i,6) =   2.99864229366191
     CASE ('HF')
        ! mm(i)       = 20.006       !  PC-SAFT
        ! parame(i,1) = 0.87731
        ! parame(i,2) = 3.0006
        ! parame(i,3) = 112.24
        ! nhb_typ(i)     = 2
        ! nhb_no(i,1)    = 1
        ! nhb_no(i,2)    = 1
        ! eps_hb(i,i,1,2)= 2208.22
        ! eps_hb(i,i,2,1)= 2208.22
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 0.71265
        mm(i)       = 20.0060000000000    !  PCP-SAFT
        parame(i,1) = 1.00030000000000
        parame(i,2) = 3.17603622195029
        parame(i,3) = 331.133373208282
        nhb_typ(i)     = 2
        nhb_no(i,1)    = 1
        nhb_no(i,2)    = 1
        eps_hb(i,i,1,2)= 348.251433080979
        eps_hb(i,i,2,1)= 348.251433080979
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 2.868887975449893E-002
        parame(i,6) = 1.82600000000000
     CASE ('HCl')
        ! mm(i)       = 36.4610000000000
        ! parame(i,1) = mm(i)*  3.922046741026943E-002
        ! parame(i,2) = 3.08731180727493
        ! parame(i,3) = 203.898845304388
        ! nhb_typ(i)     = 2
        ! nhb_no(i,1)    = 1
        ! nhb_no(i,2)    = 1
        ! eps_hb(i,i,1,2)= 245.462773177367
        ! eps_hb(i,i,2,1)= 245.462773177367
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 0.256454187330899
        mm(i)       = 36.461        ! PCIP-SAFT
        parame(i,1) = 1.6335
        parame(i,2) = 2.9066
        parame(i,3) = 190.17
        parame(i,6) = 1.1086
        IF (pol == 2) parame(i,11)= 2.63
     CASE ('gen')
        mm(i)       = 302.0
        parame(i,1) = 8.7563
        parame(i,2) = 3.604243
        parame(i,3) = 255.6434
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 0.0
        eps_hb(i,i,2,1)= 0.0
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = 0.02
     CASE ('h2o')
        mm(i)       = 18.015
        parame(i,1) = mm(i)*  .05915
        parame(i,2) = 3.00068335
        parame(i,3) = 366.512135
        nhb_typ(i)  = 2
        nhb_no(i,1) = 1            ! no. of sites of type 1
        nhb_no(i,2) = 1            ! no. of sites of type 2
        eps_hb(i,i,1,2)= 2500.6706
        eps_hb(i,i,2,1)= 2500.6706
        eps_hb(i,i,1,1)= 0.0
        eps_hb(i,i,2,2)= 0.0
        kap_hb(i,i)    = .0348679836

        ! mm(i)       = 18.015
        ! parame(i,1) = 1.706
        ! parame(i,2) = 2.429
        ! parame(i,3) = 242.19
        ! nhb_typ(i)  = 2
        ! nhb_no(i,1) = 1            ! no. of sites of type 1
        ! nhb_no(i,2) = 1            ! no. of sites of type 2
        ! eps_hb(i,i,1,2)= 2644.2
        ! eps_hb(i,i,2,1)= 2644.2
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 0.153

        ! mm(i)       = 18.015
        ! parame(i,1) = mm(i)*  .0588185709
        ! parame(i,2) = 3.02483704
        ! parame(i,3) = 382.086672
        !c! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
        ! nhb_typ(i)  = 2
        ! nhb_no(i,1) = 1            ! no. of sites of type 1
        ! nhb_no(i,2) = 2            ! no. of sites of type 2
        !c! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
        ! eps_hb(i,i,1,2)= 2442.49782
        ! eps_hb(i,i,2,1)= 2442.49782
        ! eps_hb(i,i,1,1)=0.0
        ! eps_hb(i,i,2,2)=0.0
        ! kap_hb(i,i)    = .0303754635

        ! mit gefittetem Dipol-Moment - Haarlem-night
        ! mm(i)       = 18.015
        ! parame(i,1) = mm(i)*  7.0037160952278E-2
        ! parame(i,2) = 2.79276650240763
        ! parame(i,3) = 270.970053834558
        ! nhb_typ(i)  = 2
        ! nhb_no(i,1) = 1            ! no. of sites of type 1
        ! nhb_no(i,2) = 1            ! no. of sites of type 2
        ! eps_hb(i,i,1,2)= 1427.8287
        ! eps_hb(i,i,2,1)= 1427.8287
        ! eps_hb(i,i,1,1)=0.0
        ! eps_hb(i,i,2,2)=0.0
        ! kap_hb(i,i)    = 4.335167238E-2
        ! parame(i,6) = 3.968686856378

        IF (pol >= 1) mm(i)       = 18.015       !  PCP-SAFT
        IF (pol >= 1) parame(i,1) = 0.922688825223317
        IF (pol >= 1) parame(i,2) = 3.17562052023944
        IF (pol >= 1) parame(i,3) = 388.462197714696
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2000.67247409031
        IF (pol >= 1) eps_hb(i,i,2,1)= 2000.67247409031
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 2.040614952751225E-003
        IF (pol >= 1) parame(i,6) = 1.85500000000000
        IF (pol >= 1) parame(i,7) = 2.00000000000000

        IF (pol >= 1) mm(i) = 18.015           !PCP-SAFT with m=1 enforced and Q adjusted
        IF (pol >= 1) parame(i,1) = 1.0
        IF (pol >= 1) parame(i,2) = 3.10292190714748
        IF (pol >= 1) parame(i,3) = 308.033050211038
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 2166.02739727598
        IF (pol >= 1) eps_hb(i,i,2,1)= 2166.02739727598
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i) = 3.469980265130953E-004
        IF (pol >= 1) parame(i,6) = 1.85500000000000
        IF (pol >= 1) parame(i,7) = 3.50269655874463

        IF (pol >= 1) mm(i) = 18.015           !PCP-SAFT Q adjusted (excellent fit to exp. data)
        IF (pol >= 1) parame(i,1) = 1.19897029844512
        IF (pol >= 1) parame(i,2) = 2.87559712360227
        IF (pol >= 1) parame(i,3) = 283.758844011803
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 1
        IF (pol >= 1) nhb_no(i,2)    = 1
        IF (pol >= 1) eps_hb(i,i,1,2)= 1326.09825616021
        IF (pol >= 1) eps_hb(i,i,2,1)= 1326.09825616021
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i) = 5.680937058899043E-003
        IF (pol >= 1) parame(i,6) = 1.85500000000000
        IF (pol >= 1) parame(i,7) = 3.48575663070838

        IF (pol >= 1) mm(i) = 18.015           !PCP-SAFT Q adjusted, 4 sites (excellent fit to exp. data)
        IF (pol >= 1) parame(i,1) = 1.27084414211470
        IF (pol >= 1) parame(i,2) = 2.81979213575847
        IF (pol >= 1) parame(i,3) = 281.943549171383
        IF (pol >= 1) nhb_typ(i)     = 2
        IF (pol >= 1) nhb_no(i,1)    = 2
        IF (pol >= 1) nhb_no(i,2)    = 2
        IF (pol >= 1) eps_hb(i,i,1,2)= 952.657528272410
        IF (pol >= 1) eps_hb(i,i,2,1)= 952.657528272410
        IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
        IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
        IF (pol >= 1) kap_hb(i,i)    = 4.353859441452534E-003
        IF (pol >= 1) parame(i,6) = 1.85500000000000
        IF (pol >= 1) parame(i,7) = 3.43534786736533

        ! IF (pol.EQ.2) mm(i)       = 18.015       !  PCIP-SAFT - DQ with my=my_RPT
        ! IF (pol.EQ.2) parame(i,1) = 1.0
        ! IF (pol.EQ.2) parame(i,2) = 3.14540664928026
        ! IF (pol.EQ.2) parame(i,3) = 320.283823615925
        ! IF (pol.EQ.2) nhb_typ(i)     = 2
        ! IF (pol.EQ.2) nhb_no(i,1)    = 2
        ! IF (pol.EQ.2) nhb_no(i,2)    = 2
        ! IF (pol.EQ.2) eps_hb(i,i,1,2)= 1335.72887678032
        ! IF (pol.EQ.2) eps_hb(i,i,2,1)= 1335.72887678032
        ! IF (pol.EQ.2) eps_hb(i,i,1,1)= 0.0
        ! IF (pol.EQ.2) eps_hb(i,i,2,2)= 0.0
        ! IF (pol.EQ.2) kap_hb(i,i)    = 4.162619960844732E-003
        ! IF (pol.EQ.2) parame(i,6) = 1.85500000000000
        ! IF (pol.EQ.2) parame(i,7) = 2.00000000000000
        ! IF (pol.EQ.2) parame(i,11)= 1.45000000000000
        ! mm(i)       = 18.0150000000000
        ! parame(i,1) = 1.0
        ! parame(i,2) = 3.11505069470915
        ! parame(i,3) = 320.991387913502
        ! nhb_typ(i)     = 2
        ! nhb_no(i,1)    = 1
        ! nhb_no(i,2)    = 1
        ! eps_hb(i,i,1,2)= 2037.76329812542
        ! eps_hb(i,i,2,1)= 2037.76329812542
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 3.763148982832804E-003
        ! parame(i,6) = 1.85500000000000
        ! parame(i,7) = 2.00000000000000
        ! IF (pol.EQ.2) parame(i,11)= 1.45000000000000
        IF (pol == 2) mm(i)       = 18.015       !  PCIP-SAFT - DQ with my=my_0
        IF (pol == 2) parame(i,1) = 1.0
        IF (pol == 2) parame(i,2) = 3.11574491885322
        IF (pol == 2) parame(i,3) = 322.699984283499
        IF (pol == 2) nhb_typ(i)     = 2
        IF (pol == 2) nhb_no(i,1)    = 1
        IF (pol == 2) nhb_no(i,2)    = 1
        IF (pol == 2) eps_hb(i,i,1,2)= 2033.87777692450
        IF (pol == 2) eps_hb(i,i,2,1)= 2033.87777692450
        IF (pol == 2) eps_hb(i,i,1,1)= 0.0
        IF (pol == 2) eps_hb(i,i,2,2)= 0.0
        IF (pol == 2) kap_hb(i,i)    = 3.815764667176484E-003
        IF (pol == 2) parame(i,6) = 1.85500000000000
        IF (pol == 2) parame(i,7) = 2.00000000000000
        IF (pol == 2) parame(i,11)= 1.45000000000000
        ! mm(i)       = 18.015          ! Dortmund
        ! parame(i,1) = 0.11065254*mm(i)
        ! parame(i,2) = 2.36393615
        ! parame(i,3) = 300.288589
        ! nhb_typ(i)     = 2
        ! nhb_no(i,1)    = 1
        ! nhb_no(i,2)    = 1
        ! eps_hb(i,i,1,2)= 1193.45585
        ! eps_hb(i,i,2,1)= 1193.45585
        ! eps_hb(i,i,1,1)= 0.0
        ! eps_hb(i,i,2,2)= 0.0
        ! kap_hb(i,i)    = 0.091203519
        ! parame(i,6) = 1.8546
        ! parame(i,7) = 0.0
        ! parame(i,11)= 0.0
     CASE ('MBBA')
        mm(i)     =         267.37
        parame(i,1) =        12.194
        parame(i,2) =         3.064
        parame(i,3) =       270.7
        e_lc(i,i)   =        13.7      !Hino & Prausnitz
        s_lc(i,i)   =         0.176    !Hino & Prausnitz
     CASE ('PCH5')
        mm(i)     =         255.41
        parame(i,1) =        11.6
        parame(i,2) =         3.2
        parame(i,3) =       270.7
        ! E_LC(i,i)   =        16.7      !Hino & Prausnitz
        ! S_LC(i,i)   =         0.176    !Hino & Prausnitz
        e_lc(i,i)   =        8.95
        s_lc(i,i)   =         0.2

        ! mm(i)     =         255.41
        ! parame(i,1) =        11.6
        ! parame(i,2) =         3.2
        ! parame(i,3) =       290.7
        ! E_LC(i,i)   =        7.18
        ! S_LC(i,i)   =         0.2

     CASE ('Li')
        mm(i)       = 6.9
        parame(i,1) = 1.0
        parame(i,2) = 1.4
        parame(i,3) = 96.83
        parame(i,10)= 1.0
        !   the self-association is set to zero in routine F_EXPL for ions
        ! nhb_typ(i)  = 1
        ! nhb_no(i,1) = 3                      ! no. of sites of type 1
        ! eps_hb(i,i,1,1)= 2000.0
        ! kap_hb(i,i)    = 0.008
     CASE ('Na')
        mm(i)       = 23.0
        parame(i,1) = 1.0
        parame(i,2) = 1.9
        parame(i,3) = 147.38
        parame(i,10)= 1.0
        !   the self-association is set to zero in routine F_EXPL for ions
        nhb_typ(i)  = 1
        nhb_no(i,1) = 3                       ! no. of sites of type 1
        eps_hb(i,i,1,1)= 8946.28257     ! 25C, 3 sites, dG-ref-1
        kap_hb(i,i)    = 0.001648933
     CASE ('Ka')
        mm(i)       = 39.1
        parame(i,1) = 1.0
        parame(i,2) = 2.66
        parame(i,3) = 221.44
        parame(i,10)= 1.0
        ! the self-association is set to zero in routine F_EXPL for ions
        nhb_typ(i)  = 1
        nhb_no(i,1) = 3                      ! no. of sites of type 1
        eps_hb(i,i,1,1)= 3118.336216     ! 25C, 3 sites, dG-ref-1
        kap_hb(i,i)    = 0.00200559
     CASE ('Cs')
        mm(i)       = 132.9
        parame(i,1) = 1.0
        parame(i,2) = 3.38
        parame(i,3) = 523.28
        parame(i,10)= 1.0
        ! the self-association is set to zero in routine F_EXPL for ions
        ! nhb_typ(i)  = 1
        ! nhb_no(i,1) = 3                      ! no. of sites of type 1
        ! eps_hb(i,i,1,1)= 2000.0
        ! kap_hb(i,i)    = 0.00200559
     CASE ('Cl')
        mm(i)       = 35.5
        parame(i,1) = 1.0
        parame(i,2) = 3.62
        parame(i,3) = 225.44
        parame(i,10)= -1.0
        !   the self-association is set to zero in routine F_EXPL for ions
        nhb_typ(i)  = 1
        nhb_no(i,1) = 4                      ! no. of sites of type 1
        eps_hb(i,i,1,1)= 6744.12509     ! 25C, 3 sites, dG-ref-1
        kap_hb(i,i)    = 0.00155252
     CASE ('Br')
        mm(i)       = 79.9
        parame(i,1) = 1.0
        parame(i,2) = 3.9
        parame(i,3) = 330.82
        parame(i,10)= -1.0
        !   the self-association is set to zero in routine F_EXPL for ions
        nhb_typ(i)  = 1
        nhb_no(i,1) = 4                      ! no. of sites of type 1
        eps_hb(i,i,1,1)= 4516.033227     ! 25C, 3 sites, dG-ref-1
        kap_hb(i,i)    = 0.00200559
     CASE ('Io')
        mm(i)       = 126.9
        parame(i,1) = 1.0
        parame(i,2) = 4.4
        parame(i,3) = 380.60
        parame(i,10)= -1.0
        !   the self-association is set to zero in routine F_EXPL for ions
        nhb_typ(i)  = 1
        nhb_no(i,1) = 4                      ! no. of sites of type 1
        eps_hb(i,i,1,1)= 1631.203342     ! 25C, 3 sites, dG-ref-1
        kap_hb(i,i)    = 0.00200559
     CASE ('OH')
        mm(i)       = 17.0
        parame(i,1) = 1.0
        parame(i,2) = 3.06
        parame(i,3) = 217.26
        parame(i,10)= -1.0
        ! the self-association is set to zero in routine F_EXPL for ions
        nhb_typ(i)  = 1
        nhb_no(i,1) = 4                      ! no. of sites of type 1
        eps_hb(i,i,1,1)= 14118.68089    ! 25C, 3 sites, dG-ref-1
        kap_hb(i,i)    = 0.00200559
     CASE ('NO3')
        mm(i)       = 62.0
        parame(i,1) = 1.0
        parame(i,2) = 4.12
        parame(i,3) = 239.48
        parame(i,10)= -1.0
        ! the self-association is set to zero in routine F_EXPL for ions
        ! nhb_typ(i)  = 1
        ! nhb_no(i,1) = 4                      ! no. of sites of type 1
        ! eps_hb(i,i,1,1)= 2000.0
        ! kap_hb(i,i)    = 0.00200559
     CASE ('bf4')
        mm(i)       = 86.8
        parame(i,1) = 1.0
        parame(i,2) = 4.51 ! *0.85
        parame(i,3) = 164.7
        parame(i,10)= -1.0
     CASE ('pf6')
        mm(i)       = 145.0
        parame(i,1) = 1.0
        parame(i,2) = 5.06
        parame(i,3) = 224.9
        parame(i,10)= -1.0
     CASE ('emim')
        mm(i)       = 111.16
        parame(i,1) = 3.11
        parame(i,2) = 4.0
        parame(i,3) = 250.0
        parame(i,10)= 1.0
     CASE ('bmim')
        mm(i)       = 139.21
        ! parame(i,1) = 2.81
        ! parame(i,2) = 3.5
        parame(i,1) = 3.81
        parame(i,2) = 4.0
        parame(i,3) = 250.0
        parame(i,6) = 0.0
        parame(i,10)= 1.0
     CASE ('hmim')
        mm(i)       = 167.27
        parame(i,1) = 4.53
        parame(i,2) = 4.0
        parame(i,3) = 250.0
        parame(i,10)= 1.0
     CASE ('omim')
        mm(i)       = 195.32
        parame(i,1) = 5.30
        parame(i,2) = 4.0
        parame(i,3) = 250.0
        parame(i,10)= 1.0
     CASE ('sw')
        parame(i,1) = 1.0
        parame(i,2) = 1.0
        parame(i,3) = 100.0
        parame(i,4) = 0.0
        parame(i,5) = 0.0
        mm(i)       = 1.0
        parame(i,6) = 0.1175015839*2.0
        ! use Temp. in Kelvin in the input-file. For dimensionless quantities
        ! (P*=P*sig**3/epsilon, T*=T*kBol/epsilon, rho*=rho*sig**3) calculate
        ! P* = P *1E5 * (1.e-10)^3 / (100*8.31441/6.022045E+23)
        ! T* = (T+273.15)/100
        ! for rho* go to utilities.f (subroutine SI_DENS) and write
        ! density(ph) = dense(ph)*6.0/PI
     CASE ('c8-sim')
        mm(i)       =         114.231
        parame(i,1) = mm(i)*  4.095944E-2  !  =4.67883774717337
        parame(i,2) =         3.501769
        parame(i,3) =         163.8606
        ! mm(i)       =         114.231000000000
        ! parame(i,1) = mm(i)*  3.547001476437745E-002  !  =4.05177525654960
        ! parame(i,2) =         3.70988567055814
        ! parame(i,3) =         192.787548176343
     CASE ('argon_ge')
        mm(i)       = 39.948
        parame(i,1) = mm(i)*0.030327
        parame(i,2) = 3.149910
        parame(i,3) = 100.188975
     CASE ('argon_ge2')
        mm(i)       = 39.948
        parame(i,1) = mm(i)*0.030327
        parame(i,2) = 3.149910
        parame(i,3) = 0.8*100.188975
     CASE default
        WRITE (*,*) ' pure component parameters missing for ',compna(i)
        STOP
     END SELECT

     IF ( pol == 2 .AND. parame(i,11) == 0.0 ) THEN
        WRITE (*,*) ' polarizability missing for comp. ',compna(i)
        STOP
     END IF

     IF (nhb_typ(i) /= 0) THEN
        parame(i,12) = REAL(nhb_typ(i))
        parame(i,13) = kap_hb(i,i)
        no = 0
        DO j=1,nhb_typ(i)
           DO k=1,nhb_typ(i)
              parame(i,(14+no))=eps_hb(i,i,j,k)
              no=no+1
           END DO
        END DO
        DO j=1,nhb_typ(i)
           parame(i,(14+no))=REAL(nhb_no(i,j))
           no=no+1
        END DO
     ELSE
        DO k=12,25
           parame(i,k)=0.0
        END DO
     END IF

  END DO


  DO  i = 1,ncomp

     DO  j = 1,ncomp

        IF (compna(i) == 'ps'.AND.compna(j) == 'cyclohexane')THEN
           kij(i,j) = 0.0075
        ELSE IF(compna(i) == 'peva'.AND.compna(j) == 'ethylene')THEN
           ! --   0 Gew.% VA-------------
           ! kij(i,j) = 0.039
           ! -- 7.5 Gew.% VA-------------
           ! kij(i,j) = 0.0377325
           ! kij(i,j) = 0.0374867
           ! ---12.7 Gew.% VA------------
           ! kij(i,j) = 0.036854
           ! kij(i,j) = 0.0366508
           ! ---27.3 Gew.% VA------------
           ! kij(i,j) = 0.034386
           ! kij(i,j) = 0.0352375
           ! ---31.8 Gew.% VA------------
           kij(i,j) = 0.033626
           ! kij(i,j) = 0.0350795
           ! ---42.7 Gew.% VA------------
           ! kij(i,j) = 0.031784
           ! kij(i,j) = 0.035239
        ELSE IF(compna(i) == 'gen'.AND.compna(j) == 'h2o')THEN
           kij(i,j) = - 0.2
        ELSE IF(compna(i) == 'peva'.AND.compna(j) == 'vinylacetate')THEN
           kij(i,j) = 0.019
        ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'co2') THEN
           IF ( pol == 0 ) kij(i,j) = 0.195
           IF ( pol == 1 ) kij(i,j) = 0.06
        ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'acetone') THEN
           kij(i,j) = 0.021
        ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'hexane') THEN
           kij(i,j) = 0.012
        ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'pentane')THEN
           kij(i,j) = 0.005
        ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'methylcyclohexane') THEN
           kij(i,j) = 0.0073
        ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'ethylbenzene')THEN
           kij(i,j) = 0.008
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = 0.181
           kij(i,j) = 0.088
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'propane') THEN
           kij(i,j) = 0.0206
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'butane') THEN
           kij(i,j) = 0.01
        ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'argon') THEN
           kij(i,j) = 0.01
        ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'butane') THEN
           kij(i,j) = 0.026
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'pentane') THEN
           ! kij(i,j) = -0.0195
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'hexane') THEN
           ! kij(i,j) = 0.008
           kij(i,j) = 0.004
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'ethylene') THEN
           kij(i,j) = 0.0404
           ! kij(i,j) = 0.0423
           ! kij(i,j) = 0.044
        ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'cyclohexane') THEN
           kij(i,j) = -0.1
        ELSE IF(compna(i) == 'ldpe'.AND.compna(j) == 'cyclopentane')THEN
           kij(i,j) = -0.016
        ELSE IF(compna(i) == 'pp'.AND.compna(j) == 'propane') THEN
           kij(i,j) = 0.0242
        ELSE IF(compna(i) == 'pp'.AND.compna(j) == 'pentane') THEN
           kij(i,j) = 0.0137583176
        ELSE IF(compna(i) == 'pp'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.1767         ! without quadrupol-term
           kij(i,j) = 0.063         ! with quadrupol-term
        ELSE IF(compna(i) == 'pba'.AND.compna(j) == 'ethylene') THEN
           kij(i,j) = 0.026
        ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'ethane') THEN
           kij(i,j) = 0.017
        ELSE IF(compna(i) == 'n2'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = -0.04
        ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.051875  ! PC-SAFT
           IF (pol == 1) kij(i,j) = -0.0353125  ! PCP-SAFT
        ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'co') THEN
           ! IF (pol == 1) kij(i,j) = -0.003  ! PCP-SAFT
           IF (pol == 1) kij(i,j) = 0.018  ! PCP-SAFT
        ELSE IF(compna(i) == 'ethane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.095
           kij(i,j) = 0.021
           ! kij(i,j) = 0.024
        ELSE IF(compna(i) == 'propane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.042
        ELSE IF(compna(i) == 'argon_ge'.AND.compna(j) == 'argon_ge2') THEN
           read (*,*) kij(i,j)
        ELSE IF(compna(i) == 'butane'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = 0.115
           ! kij(i,j) = 0.048
           kij(i,j) = 0.036
        ELSE IF(compna(i) == 'pentane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.143         ! without quadrupol-term
           kij(i,j) = 0.0            ! with quadrupol-term
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = 0.125         ! without quadrupol-term
           kij(i,j) = 0.0495
        ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.11         ! without quadrupol-term
           ! kij(i,j) = 0.05
           ! kij(i,j) = 0.039         ! with quadrupol-term
        ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = 0.128         ! without quadrupol-term
           kij(i,j) = 0.053         ! with quadrupol-term
        ELSE IF(compna(i) == 'dodecane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.12         ! without quadrupol-term
           kij(i,j) = 0.0508         ! with quadrupol-term
        ELSE IF(compna(i) == 'benzene'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = 0.087968750000    ! without quadrupol-term
           ! kij(i,j) =  0.008203125   ! only co2 quadrupol
           kij(i,j) = 0.042     ! both quadrupol
           ! kij(i,j) = 0.003     ! both quadrupol
        ELSE IF(compna(i) == 'toluene'.AND.compna(j) == 'co2') THEN
           ! kij(i,j) = 0.110784912  ! without quadrupol-term
           kij(i,j) = 0.0305        ! with quadrupol-term
        ELSE IF(compna(i) == 'cyclohexane'.AND.compna(j) == 'co2') THEN
           kij(i,j) = 0.13
           lij(i,j) = - 0.00
           ! kij(i,j) = 0.045
        ELSE IF(compna(i) == 'chloromethane'.AND.compna(j) == 'co2')THEN
           ! kij(i,j) = 0.04  ! PC-SAFT
           kij(i,j) = 0.025  ! PCP-SAFT
           ! kij(i,j) = 0.083  ! PCIP-SAFT
        ELSE IF(compna(i) == 'acetone'.AND.compna(j) == 'n2')THEN
           kij(i,j) = 0.035211  ! PCP-SAFT
           lij(i,j) = + 0.013225  ! PCP-SAFT
           kij(i,j) = 0.023  ! PCP-SAFT
           lij(i,j) = + 0.013225  ! PCP-SAFT
           kij(i,j) = 0.037676  ! PCP-SAFT (with more emphasis on crit. pt.)
           lij(i,j) = 0.0  ! PCP-SAFT (with more emphasis on crit. pt.)
           !kij(i,j) = 1.722238535635467E-002  ! PCP-SAFT
           !lij(i,j) = 2.815974678394451E-003  ! PCP-SAFT
           !kij(i,j) = 1.931522058164026E-002  ! PCP-SAFT
           !lij(i,j) = 0.0  ! PCP-SAFT
           !kij(i,j) =  1.641053794134795E-002  ! PCP-SAFT
           !lij(i,j) = -5.850421759950764E-003  ! PCP-SAFT
           if ( num == 0 ) write (*,*) 'calculation with lij only possible with num=1'
           if ( num == 0 ) stop 5
        ELSE IF(compna(i) == 'acetone'.AND.compna(j) == 'co2')THEN
           kij(i,j) = 0.015  ! PC-SAFT
           IF (pol == 1) kij(i,j) = -0.02  ! PCP-SAFT
           IF (pol == 2) kij(i,j) = -0.005  ! PCIP-SAFT  where DQ with my=my_vacuum
           ! IF (pol.EQ.2) kij(i,j) = 0.0  ! PCIP-SAFT  where DQ with my=my_RPT
        ELSE IF(compna(i) == 'methanol'.AND.compna(j) == 'co2')THEN
           ! kij(i,j) =  0.0288  ! PC-SAFT
           ! kij(i,j) = - 0.035  ! PCP-SAFT for co2 and PC-SAFT methanol
           ! kij(i,j) = - 0.035  ! PCP-SAFT
           ! lij(i,j) =  0.3  ! PCP-SAFT
        ELSE IF(compna(i) == 'dimethyl-ether'.AND.compna(j) == 'co2')THEN
           kij(i,j) = 0.00896894  ! PC-SAFT
           ! kij(i,j) = - 0.015  ! PCP-SAFT
        ELSE IF(compna(i) == 'dimethyl-ether'.AND.compna(j) == 'h2o')THEN
           ! kij(i,j) = -0.134  ! PC-SAFT
        ELSE IF(compna(i) == 'dichloromethane'.AND.compna(j) == 'co2')THEN
           ! kij(i,j) = 0.06881725  ! PC-SAFT
           ! kij(i,j) = 0.05839145  ! PCP-SAFT
           kij(i,j) = -0.00944346  ! PCP-SAFT co2 dichloromethane PC-SAFT
           ! kij(i,j) = 0.06  ! PCIP-SAFT
        ELSE IF(compna(i) == 'h2s'.AND.compna(j) == 'methane')THEN
           ! kij(i,j) =  0.0414 ! PC-SAFT
           kij(i,j) =  0.0152 ! PCP-SAFT Dipole momnet (d with Q)
        ELSE IF(compna(i) == 'butane'.AND.compna(j) == 'h2s')THEN
           kij(i,j) = -0.002  ! PCP-SAFT
        ELSE IF(compna(i) == 'methanol'.AND.compna(j) == 'h2s')THEN
           kij(i,j) = 0.0  ! PCP-SAFT
           lij(i,j) = 0.0  ! PCP-SAFT
        ELSE IF(compna(i) == 'co2'.AND.compna(j) == 'hydrogen') THEN
           ! lij(i,j) = -0.08 !!!!! Lij not kij !!!!
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'n2') THEN
           kij(i,j) = 0.11
        ELSE IF(compna(i) == 'propane'.AND.compna(j) == 'n2') THEN
           kij(i,j) = 0.0251171875
        ELSE IF(compna(i) == 'co2'.AND.compna(j) == 'hexadecane') THEN
           ! kij(i,j) = 0.1194   ! PC-SAFT ohne QQ
           kij(i,j) = 0.0588
        ELSE IF(compna(i) == 'ethane'.AND.compna(j) == 'acetone') THEN
           ! kij(i,j) = 0.065  ! no DD
           kij(i,j) = 0.038   ! DD non-polarizable
           ! kij(i,j) = 0.025   ! DD polarizable
        ELSE IF(compna(i) == 'butane'.AND.compna(j) == 'acetone') THEN
           ! kij(i,j) = 0.065  ! no DD
           kij(i,j) = 0.037  ! DD non-polarizable
           ! kij(i,j) = 0.025  ! DD polarizable
        ELSE IF(compna(i) == 'pentane'.AND.compna(j) == 'acetone') THEN
           ! kij(i,j) = 0.072  ! no DD
           ! kij(i,j) = 0.041  ! DD non-polarizable
           kij(i,j) = 0.039  ! DD polarizable
           ! kij(i,j) = 0.035  ! DD polarizable
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'acetone') THEN
           ! kij(i,j) = 0.063
           kij(i,j) = 0.038  ! PCP-SAFT
           ! kij(i,j) = 0.036  ! PCIP-SAFT
        ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'acetone') THEN
           kij(i,j) = 0.035  ! PCP-SAFT
        ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'acetone') THEN
           ! kij(i,j) = 0.059  ! no DD
           ! kij(i,j) =  0.03281250  ! DD non-polarizable
           kij(i,j) = 0.028  ! DD polarizable
        ELSE IF(compna(i) == 'hexadecane'.AND.compna(j) == 'acetone') THEN
           ! kij(i,j) = 0.07  ! no DD
           ! kij(i,j) =  0.0415  ! DD non-polarizable
           kij(i,j) = 0.035  ! DD polarizable
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'butanone')THEN
           kij(i,j) = 0.027   ! PCP-SAFT
           ! kij(i,j) = 0.033   ! PCP-SAFT with lij
           ! lij(i,j) = 0.13    ! PCP-SAFT
           ! kij(i,j) = 0.042    ! PC-SAFT
        ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'butanone')THEN
           kij(i,j) = 0.042  ! no DD
           ! kij(i,j) = 0.027  ! DD non-polarizable
        ELSE IF(compna(i) == 'heptane'.AND.compna(j) == '2-pentanone')THEN
           kij(i,j) = 0.041  ! no DD
           ! kij(i,j) = 0.031  ! DD non-polarizable
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == '3-pentanone')THEN
           kij(i,j) = 0.0
        ELSE IF(compna(i) == 'pentane'.AND.compna(j) == 'propanal')THEN
           ! kij(i,j) = 0.055  ! no DD
           kij(i,j) = 0.027  ! DD non-polarizable
           ! kij(i,j) = 0.026  ! DD polarizable 22
        ELSE IF(compna(i) == 'cyclohexane'.AND.compna(j) == 'propanal')THEN
           ! kij(i,j) = 0.06  ! no DD
           kij(i,j) = 0.036  ! DD non-polarizable
           kij(i,j) = 0.035  ! DD polarizable 22
        ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'butanal')THEN
           kij(i,j) = 0.041  ! no DD
           ! kij(i,j) = 0.025  ! DD non-polarizable
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'thf')THEN
           kij(i,j) = 0.012  ! PCP-SAFT
        ELSE IF(compna(i) == 'octane'.AND.compna(j) == 'thf')THEN
           kij(i,j) = 0.012  ! PCP-SAFT
        ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'thf')THEN
           kij(i,j) = 0.012  ! PCP-SAFT
        ELSE IF(compna(i) == 'toluene'.AND.compna(j) == 'dmso')THEN
           ! kij(i,j) = 0.025  ! no DD
           kij(i,j) = - 0.0105  ! DD non-polarizable
           ! kij(i,j) = - 0.019  ! DD polarizable
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'acrylonitrile')THEN
           kij(i,j) = - 0.05  ! DD polarizable
        ELSE IF(compna(i) == 'heptane' .AND. compna(j) == 'butyronitrile')THEN
           kij(i,j) = - 0.002  ! DD polarizable 11
           kij(i,j) = 0.002  ! DD polarizable 22
        ELSE IF(compna(i) == '1-butene'.AND.compna(j) == 'dmf')THEN
           ! kij(i,j) = 0.04  ! no DD
           ! kij(i,j) = 0.004  ! DD non-polarizable
           kij(i,j) = 0.005  ! DD polarizable 22
        ELSE IF(compna(i) == 'cyclohexane'.AND.compna(j) == 'dmf')THEN
           kij(i,j) = 0.0135  ! DD polarizable 11
           kij(i,j) = 0.022  ! DD polarizable 22
        ELSE IF(compna(i) == 'ethylene'.AND.compna(j) == 'dmf')THEN
           ! kij(i,j) = - 0.0215  ! DD polarizable 11
           kij(i,j) = - 0.01  ! DD polarizable 22
        ELSE IF(compna(i) == 'nbutyl-ethanoate'.AND.compna(j) == 'dmf')THEN
           ! kij(i,j) = 0.016  ! no DD
           ! kij(i,j) = -0.01  ! DD non-polarizable
           kij(i,j) = - 0.015  ! DD polarizable 22
        ELSE IF(compna(i) == 'methylacetate' .AND. compna(j) == 'cyclohexane')THEN
           kij(i,j) = 0.066  ! PC-SAFT
           ! kij(i,j) = 0.061  ! PCP-SAFT
           ! kij(i,j) = 0.0625  ! PCIP-SAFT
        ELSE IF(compna(i) == 'methylacetate'.AND.compna(j) == 'decane')THEN
           kij(i,j) = 0.0625  ! PCIP-SAFT
        ELSE IF(compna(i) == 'methylacetate' .AND. compna(j) == 'methanol')THEN
           ! kij(i,j) = -0.07  ! PCIP-SAFT
        ELSE IF(compna(i) == 'pentane' .AND. compna(j) == 'propionitrile')THEN
           kij(i,j) =  0.0498
           IF (pol >= 1) kij(i,j) = -0.01
           IF (pol >= 2) kij(i,j) = -0.027
        ELSE IF(compna(i) == 'hexane' .AND. compna(j) == 'propionitrile')THEN
           kij(i,j) = 0.05
           IF (pol >= 1) kij(i,j) = 0.0
           IF (pol >= 2) kij(i,j) = -0.03
        ELSE IF(compna(i) == 'octane' .AND. compna(j) == 'propionitrile')THEN
           kij(i,j) = 0.0  ! DD polarizable 22
        ELSE IF(compna(i) == 'cyclohexane' .AND. compna(j) == 'nitromethane')THEN
           kij(i,j) =  0.14  ! no DD
           ! kij(i,j) =   0.07  ! DD non-polarizable
           ! kij(i,j) =   0.055  ! DD polarizable 22
        ELSE IF(compna(i) == 'cyclohexane' .AND. compna(j) == 'nitroethane')THEN
           ! kij(i,j) =  0.06  ! no DD
           kij(i,j) =   0.03  ! DD polarizable 22
        ELSE IF(compna(i) == 'acetone' .AND. compna(j) == 'nitromethane')THEN
           ! kij(i,j) = - 0.017  ! no DD
           kij(i,j) = - 0.021  ! DD non-polarizable
           ! kij(i,j) = - 0.023  ! DD polarizable 22
        ELSE IF(compna(i) == 'acetone' .AND. compna(j) == 'h2o')THEN
           kij(i,j) = - 0.2    ! PCP-SAFT (no cross-association)
        ELSE IF(compna(i) == 'methylcyclohexane' .AND. compna(j) == 'acetonitrile')THEN
           ! kij(i,j) =   0.09  ! no DD
           ! kij(i,j) =   0.033  ! DD non-polarizable
           ! kij(i,j) =   0.025  ! DD polarizable 22
           kij(i,j) =   0.04  ! DD polarizable 22 und my angepasst
        ELSE IF(compna(i) == 'ethylacetate' .AND. compna(j) == 'h2o')THEN
           if ( pol == 1) kij(i,j) =  - 0.08
        ELSE IF(compna(i) == 'ethylacetate' .AND. compna(j) == 'acetonitrile')THEN
           kij(i,j) =   0.007  ! no DD
           ! kij(i,j) =   -0.045  ! DD polarizable 22
        ELSE IF(compna(i) == 'dimethyl-ether' .AND. compna(j) == 'propane')THEN
           ! kij(i,j) = 0.03  ! no DD
           kij(i,j) =  0.0225  ! DD non-polarizable
        ELSE IF(compna(i) == 'dimethyl-ether' .AND. compna(j) == 'butane')THEN
           kij(i,j) =  0.0295  ! DD non-polarizable
        ELSE IF(compna(i) == 'benzene' .AND. compna(j) == 'pentane')THEN
           ! kij(i,j) =  0.012968750   ! ohne QQ
           ! kij(i,j) =  0.004921875   ! mit QQ=5.6D (gefittet)
           ! kij(i,j) = -0.006406250   ! mit QQ=7.45D (Literatur)
        ELSE IF(compna(i) == 'benzene' .AND. compna(j) == 'heptane')THEN
           kij(i,j) = 0.01328125
           ! kij(i,j) = 0.0038
        ELSE IF(compna(i) == 'benzene' .AND. compna(j) == '1-hexene')THEN
           kij(i,j) = 0.0067
        ELSE IF(compna(i) == 'ethylene' .AND. compna(j) == 'co2') THEN
           kij(i,j) = 0.04
           kij(i,j) = -0.029
        ELSE IF(compna(i) == 'ethylene' .AND. compna(j) == 'vinylacetate') THEN
           kij(i,j) = - 0.013847
        ELSE IF(compna(i) == 'triacontane' .AND. compna(j) == 'ethylene') THEN
           kij(i,j) =  0.028
        ELSE IF(compna(i) == 'triacontane' .AND. compna(j) == 'methane')THEN
           ! kij(i,j) =  0.061953125   ! polyethylene parameters
           kij(i,j) =  0.039609375   ! param. by extrapolation of n-alkanes
        ELSE IF(compna(i) == 'tetracontane' .AND. compna(j) == 'methane')THEN
           ! kij(i,j) =  0.06515625    ! polyethylene parameters
           kij(i,j) =  0.04453125    ! param. by extrapolation of n-alkanes
           ! --- kij and lij adjusted -------
           ! kij(i,j) =  0.045786119   ! param. by extrapolation of n-alkanes
           ! lij(i,j) =  +8.53466437d-4   ! param. by extrapolation of n-alkanes
        ELSE IF(compna(i) == 'eicosane' .AND. compna(j) == 'methane')THEN
           kij(i,j) =  0.0360134457445
        ELSE IF(compna(i) == 'tetracontane' .AND. compna(j) == 'methane')THEN
           kij(i,j) =  0.0360134457445  ! assumed equal to eicosane-C1
        ELSE IF(compna(i) == 'chlorobenzene' .AND. compna(j) == 'cyclohexane') THEN
           kij(i,j) = 0.013
        ELSE IF(compna(i) == 'chloroethane' .AND. compna(j) == 'butane')THEN
           kij(i,j) = 0.025
        ELSE IF(compna(i) == 'tetrachloromethane' .AND. compna(j) == '2-methylpentane') THEN
           kij(i,j) = 0.0070105
        ELSE IF(compna(i) == 'tetrachloromethane' .AND. compna(j) == 'hexane') THEN
           kij(i,j) = 0.004
        ELSE IF(compna(i) == 'hydrogen' .AND. compna(j) == 'hexane') THEN
           kij(i,j) = 0.1501
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'co2') THEN
           ! kij(i,j) = -0.08
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'propane') THEN
           kij(i,j) = - 0.07
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'ethane') THEN
           kij(i,j) = 0.0
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'butane') THEN
           kij(i,j) = 0.028
           kij(i,j) = 0.016
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'cyclohexane')THEN
           kij(i,j) = 0.037
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == '2-methylpentane') THEN
           kij(i,j) = 0.028
        ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == '1-octanol')THEN
           kij(i,j) = 0.06
        ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'cyclohexane') THEN
           kij(i,j) = 0.0508
           ! kij(i,j) = 0.03
        ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'heptane')THEN
           kij(i,j) = 0.034
        ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'octane')THEN
           kij(i,j) =  1.79428322E-002
           lij(i,j)= 4.022964E-005
           if ( i > j ) lij(i,j) = - lij(i,j)
           if ( i > j ) write (*,*) 'sign of lij entry depends on index. Changed sign !'
           if ( i > j ) read (*,*)
        ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'decane')THEN
           ! kij(i,j) = 0.042   !  PC-SAFT
           kij(i,j) = 0.0126865   !  PCP-SAFT
           ! kij(i,j) = 0.000    !  PCIP-SAFT
        ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'isobutane') THEN
           kij(i,j) = 0.051
        ELSE IF(compna(i) == 'methanol' .AND. compna(j) == '1-octanol') THEN
           kij(i,j) = 0.02
        ELSE IF(compna(i) == '1-butanol' .AND. compna(j) == 'butane') THEN
           kij(i,j) = 0.015
        ELSE IF(compna(i) == '1-nonanol' .AND. compna(j) == 'co2') THEN
           kij(i,j) = 0.076
           kij(i,j) = 0.01443
        ELSE IF(compna(i) == '1-propanol' .AND. compna(j) == 'ethylbenzene') THEN
           kij(i,j) = 0.0251757813
        ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'ethanol') THEN
           kij(i,j) = 0.085
        ELSE IF(compna(i) == 'hexane'.AND.compna(j) == '1-chlorobutane') THEN
           kij(i,j) = 0.017
        ELSE IF(compna(i) == 'aniline'.AND.compna(j) == 'methylcyclopentane') THEN
           kij(i,j) = 0.0153
        ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'nbutyl-ethanoate')THEN
           kij(i,j) = 0.027
        ELSE IF(compna(i) == '1-hexene'.AND.compna(j) == 'ethyl-ethanoate')THEN
           kij(i,j) = 0.026
        ELSE IF(compna(i) == 'co2'.AND.compna(j) == '1-butanol')THEN
           ! kij(i,j) = 0.075
           kij(i,j) = 0.0
        ELSE IF(compna(i) == 'acrylic-acid'.AND.compna(j) == 'co2')THEN
           kij(i,j) = -0.1
        ELSE IF(compna(i) == 'bmim'.AND.compna(j) == 'hydrogen')THEN
           lij(i,j) =  0.55   !!!!! Lij not kij !!!!
        ELSE IF(compna(i) == 'bf4'.AND.compna(j) == 'hydrogen')THEN
           lij(i,j) =  0.55   !!!!! Lij not kij !!!!
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'butane')THEN !K.R.Hall FPE 2007 254 112-125 kij=0.1850
           kij(i,j) = -0.07
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == '1-butanol')THEN
           kij(i,j) = -0.015
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'aniline')THEN
           kij(i,j) = 0.0
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'methanol') THEN
           kij(i,j) = -0.02
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'ethanol') THEN
           kij(i,j) = -0.027
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'styrene') THEN
           kij(i,j) = 0.1
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'propyl-ethanoate') THEN
           kij(i,j) = -0.205
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'ethyl-propanoate') THEN
           kij(i,j) = 0.0
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == '1-pentanol') THEN
           kij(i,j) = 0.0165
           kij(i,j) = -0.015
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == '1-nonanol') THEN
           kij(i,j) = 0.0
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'methane') THEN
           ! kij(i,j) = +0.06
           kij(i,j) = -0.08
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'propane') THEN
           kij(i,j) = 0.02
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'hexane') THEN
           kij(i,j) = 0.2
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'acetic-acid') THEN
           kij(i,j) = -0.0
        ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'co2') THEN
           if (pol == 0) kij(i,j) = 0.0030625  ! for T=50C, according to X.Tang
           ! stop 5                                ! very T-dependent
        ELSE IF(compna(i) == 'toluene'.AND.compna(j) == 'acetic-acid') THEN
           kij(i,j) = -0.1
        ELSE IF(compna(i) == 'caproic-acid'.AND.compna(j) == 'cyclohexane') THEN
           kij(i,j) = 0.041531
        ELSE IF(compna(i) == '1-octanol'.AND.compna(j) == 'h2o')THEN
           kij(i,j) = 0.07
        ELSE IF(compna(i) == 'acetone'.AND.compna(j) == 'benzene')THEN
           kij(i,j) = 0.02132466    ! PC-SAFT
           ! kij(i,j) = 0.01495148    ! PCP-SAFT
           ! kij(i,j) = -0.00735575    ! PCP-SAFT but non-polar benzene
        ELSE IF(compna(i) == '1-propanol'.AND.compna(j) == 'benzene')THEN
           kij(i,j) = 0.02017
        ELSE IF(compna(i) == '2-propanol'.AND.compna(j) == 'benzene')THEN
           kij(i,j) = 0.021386
        ELSE IF(compna(i) == '1-pentanol'.AND.compna(j) == 'benzene')THEN
           kij(i,j) = 0.0129638671875
        ELSE IF(compna(i) == 'ethanediol'.AND.compna(j) == 'dodecane')THEN
           kij(i,j) = 0.08
        ELSE IF(compna(i) == 'CH2F2' .AND. compna(j) == 'co2') THEN
           kij(i,j) =  2.2548828125E-2
        ELSE IF(compna(i) == 'dmso' .AND. compna(j) == 'co2') THEN
           kij(i,j) =  -0.00
        ELSE IF(compna(i) == 'dmf'.AND.compna(j) == 'h2o')THEN
           kij(i,j) =  -0.0
        ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'h2o')THEN
           kij(i,j) =  0.11
        ELSE IF(compna(i) == '11difluoroethane'.AND.compna(j) == 'HF')THEN
           kij(i,j) =  0.032 !  association:   eps_kij = 0.16
        ELSE IF(compna(i) == '11difluoroethane'.AND.compna(j) == 'co2')THEN
           kij(i,j) =  -0.004 !  PCP-SAFT (taken from simulation)
        ELSE IF(compna(i) == 'difluoromethane'.AND.compna(j) == 'HF')THEN
           kij(i,j) = -0.02
        ELSE IF(compna(i) == 'naphthalene'.AND.compna(j) == 'co2')THEN
           kij(i,j) =  0.137  ! angepasst an SVLE-Linie (tradition.CO2-Parameter)
           kij(i,j) =  0.09  ! angepasst an SVLE-Linie (tradition.CO2-Parameter)
        ELSE IF(compna(i) == 'pg2'.AND.compna(j) == 'methanol')THEN
           kij(i,j) =  0.03
        ELSE IF(compna(i) == 'pg2'.AND.compna(j) == 'co2')THEN
           ! kij(i,j) =  0.05
        ELSE IF(compna(i) == 'PCH5'.AND.compna(j) == 'co2')THEN
           ! kij(i,j) =  -0.047
           kij(i,j) =  +0.055
           ! kij(i,j) =  +0.036
       ELSE IF (compna(i) == 'pu' .AND. compna(j) == 'cyclopentane') THEN
           kij(i,j) = 8.7268282271172337E-002 !adjusted to solubility
        ELSE
        END IF

        kij(j,i) = kij(i,j)
        lij(j,i) = -lij(i,j)

     END DO

  END DO

END SUBROUTINE pcsaft_par
