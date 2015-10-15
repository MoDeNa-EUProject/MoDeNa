!> @file
!! Main algorithm of for diffusion of gases in the foam.
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging

!*******************************************************************************************
!0 pointer allocated,
!> model supplied to the integrator
subroutine modelPU(neq, time, ystate, yprime)	! ODEPACK call
! 12/07/23 - change of model to describe diffusion on H2 and N2 throug PS, wall
!	is discretized into finite elements, pressure in cell is in equilibrium to the
!	pressure on the wall, solubility given by Henry's law cpol = H*cgas
!	- not ordinary definition of H but the definition of H fits
!	!!!this model is for filling the foam
! 2015/03/12 change in model with respect to ODEPACK calling
! michal.vonka@seznam.cz
    implicit none
    integer neq							! = 1 (X-direction)

    double precision :: time
    double precision :: ystate(*)		!
    double precision :: yprime(*)		!

    integer i, j
    integer job
    integer nFV, onecell, ncell

    double precision Dair, DCO2, Dpent     ! difusivities in polymer
    double precision Sair, SCO2, Spent     ! solubilities
    double precision HHair, HHCO2, HHpent   ! henry law cPOL = S*p*cGAS = H*cGAS
    double precision D              ! diffusivity in gas
    double precision pBCair, pBCCO2, pBCpent 	! boundary conds

    double precision dcell, dwall, hwall
    double precision cwg
    double precision RT, pressure

    double precision, allocatable :: cAIR(:) , cCO2(:), cPENT(:) !concentrations
    double precision, allocatable :: dcAIR(:), dcCO2(:), dcPENT(:) ! changes of partial pressures
    double precision, allocatable :: jAIR(:) , jCO2(:),jPENT(:)  ! fluxes in x, dir

    !c
    dcell   = ystate(nEQ + 1 )
    dwall   = ystate(nEQ + 2 )
    !tend    = ystate(nEQ + 3 )
    Dair    = ystate(nEQ + 4 ) != rpar(4) != Dair average air
    DCO2    = ystate(nEQ + 5 ) != rpar(5) != DCO2
    Sair    = ystate(nEQ + 6 ) != rpar(6) !  Sair
    SCO2    = ystate(nEQ + 7 ) != rpar(7) != SCO2
    pressure= ystate(nEQ + 8 ) != rpar(8) != pressure
    !ystate(nEQ + 9 ) = rpar(9) != initpressure
    RT      = ystate(nEQ + 10) != rpar(10)!= R*T

    ncell = dint(ystate(nEQ + 11))
    onecell = dint(ystate(nEQ + 12)) != dble(ipar(4))

    Dpent = ystate(nEQ + 13)
    Spent = ystate(nEQ + 14)

    D = ystate(nEQ + 15)    ! = rpar(13) != Dgas	  = 2.0d-5      !
    nFV  = onecell*ncell

    ! boundary conditions

    pBCair = ystate(nEQ + 16) != rpar(14) != pBCair
    pBCCO2 = ystate(nEQ + 17) != rpar(15) != pBCCO2
    pBCpent= ystate(nEQ + 18) != rpar(16) != pBCpent

    continue
    !c
    !c ------------------------------------------------------------
    !c
    if (.NOT. allocated(cAIR)) then
        allocate ( cAIR(0:nFV))
        allocate (dcAIR(1:nFV))
        allocate (jAIR (0:nFV))
        allocate ( cCO2(0:nFV))
        allocate (dcCO2(1:nFV))
        allocate (jCO2 (0:nFV))	! fluxes
        allocate ( cPENT(0:nFV))
        allocate (dcPENT(1:nFV))
        allocate (jPENT (0:nFV))
    endif

    !c        zero fluxes
    jAIR(0:nFV) = 0.0d0
    jCO2(0:nFV) = 0.0d0
    jPENT(0:nFV) = 0.0d0

    dcAIR(1:nFV) = 0.0d0
    dcCO2(1:nFV) = 0.0d0
    dcPENT(1:nFV) = 0.0d0
    !c	unpack state vars
    job = 1

    cAIR(1:nFV) = ystate(1:nFV)
    cCO2(1:nFV) = ystate(nFV+1:2*nFV)
    cPENT(1:nFV) = ystate(2*nFV+1:3*nFV)
    !c ************************************************************************
    !c

    hwall = dwall/dfloat(onecell-1)	! size of FV in wall, ocell-1 elements in wall

    !c
    !c ************************************************************************
    !c	Diffusion of AIR
    !c
    cAIR(0) = pBCair/RT !pressure/RT		! boudary condition
    do j = 0, nFV-1
        if (mod(j,onecell).eq.0) then   ! from gas to polymer
            HHair = cAIR(j)*RT*Sair     ! H = p * S
            cwg = (D*hwall*cAIR(j) + Dair*dcell*cAIR(j+1))/(D*hwall+Dair*dcell*HHair)
            jAIR(j) = D*(cAIR(j)-cwg)/dcell*2.0d0
            continue
        endif
        if ((mod(j,onecell).ne.0).and.(mod(j,onecell).ne.(onecell-1))) then ! through polymer
            jAIR(j) = Dair*(cAIR(j)-cAIR(j+1))/hwall ! j = D*(Left-Right)/h = - D*(Right-Left)/h
            continue
        endif
        if (mod(j,onecell).eq.(onecell-1)) then ! to gas in cell
            HHair = cAIR(j+1)*RT*Sair   ! H = p * S, p = c(in next cell)* RT
            cwg= (Dair*dcell*cAIR(j)+D*hwall*cAIR(j+1))/(Dair*HHair*dcell+D*hwall)
            jAIR(j) = D*(cwg-cAIR(j+1))/dcell*2.0d0
            continue
        endif
    enddo
    jAIR(nFV) = 0.0d0		! last cell, zero flux
    continue

    !
    !c BALANCES AIR
    do j = 1, nFV
        if (mod(j,onecell).eq.0) then   ! balance in cells
            dcAIR(j) = (jAIR(j-1) - jAIR(j))/dcell
        else                            ! balance in walls
            dcAIR(j) = (jAIR(j-1) - jAIR(j))/hwall
        endif
    enddo
    ! save into the yprime
    yprime(1:nFV) = dcAIR(1:nFV)

    !c ************************************************************************
    !!c	Diffusion of CO2
    !!
    cCO2(0) = pBCCO2/RT  ! boudary condition
    do j = 0, nFV-1
        if (mod(j,onecell).eq.0) then   ! from gas to polymer
            HHCO2 = cCO2(j)*RT*SCO2     ! H = p * S
            cwg = (D*hwall*cCO2(j) + DCO2*dcell*cCO2(j+1))/(D*hwall+DCO2*dcell*HHCO2)
            jCO2(j) = D*(cCO2(j)-cwg)/dcell*2.0d0
            continue
        endif
        if ((mod(j,onecell).ne.0).and.(mod(j,onecell).ne.(onecell-1))) then ! through polymer
            jCO2(j) = DCO2*(cCO2(j)-cCO2(j+1))/hwall ! j = D*(Left-Right)/h = - D*(Right-Left)/h
            continue
        endif
        if (mod(j,onecell).eq.(onecell-1)) then ! to gas in cell
            HHCO2 = cCO2(j+1)*RT*SCO2     ! H = p * S, p = c(in next cell)* RT
            cwg= (DCO2*dcell*cCO2(j)+D*hwall*cCO2(j+1))/(DCO2*HHCO2*dcell+D*hwall)
            jCO2(j) = D*(cwg-cCO2(j+1))/dcell*2.0d0
            continue
        endif
    enddo
    jCO2(nFV) = 0.0d0		! last cell, zero flux
    continue

    !
    !c BALANCES AIR
    do j = 1, nFV
        if (mod(j,onecell).eq.0) then   ! balance in cells
            dcCO2(j) = (jCO2(j-1) - jCO2(j))/dcell
        else                            ! balance in walls
            dcCO2(j) = (jCO2(j-1) - jCO2(j))/hwall
        endif
    enddo
    ! save into the yprime
    yprime(nFV+1:2*nFV) = dcCO2(1:nFV)


    !c
    !c ************************************************************************
    !c	Diffusion of PENTANE
    !c
    cPENT(0) = pBCpent/RT	! boudary condition
    do j = 0, nFV-1
        if (mod(j,onecell).eq.0) then   ! from gas to polymer
            HHpent = cPENT(j)*RT*Spent ! H = p * S
            cwg = (D*hwall*cPENT(j) + Dpent*dcell*cPENT(j+1))/(D*hwall+Dpent*dcell*HHpent)
            jPENT(j) = D*(cPENT(j)-cwg)/dcell*2.0d0
            continue
        endif
        if ((mod(j,onecell).ne.0).and.(mod(j,onecell).ne.(onecell-1))) then ! through polymer
            jPENT(j) = Dpent*(cPENT(j)-cPENT(j+1))/hwall ! j = D*(Left-Right)/h = - D*(Right-Left)/h
            continue
        endif
        if (mod(j,onecell).eq.(onecell-1)) then ! to gas in cell
            HHpent = cPENT(j+1)*RT*Spent   ! H = p * S, p = c(in next cell)* RT
            cwg= (Dpent*dcell*cPENT(j)+D*hwall*cPENT(j+1))/(Dpent*HHpent*dcell+D*hwall)
            jPENT(j) = D*(cwg-cPENT(j+1))/dcell*2.0d0
            continue
        endif
    enddo
    jPENT(nFV) = 0.0d0		! last cell, zero flux
    continue

    !
    !c BALANCES PENTANE
    do j = 1, nFV
        if (mod(j,onecell).eq.0) then   ! balance in cells
            dcPENT(j) = (jPENT(j-1) - jPENT(j))/dcell
        else                            ! balance in walls
            dcPENT(j) = (jPENT(j-1) - jPENT(j))/hwall
        endif
    enddo
    ! save into the yprime
    yprime(2*nFV+1:3*nFV) = dcPENT(1:nFV)

    continue
return
end


subroutine jdem (neq, t, y, j, ia, ja, pdj)
!-----------------------------------------------------------------------
! Should introduce own sparcity model, must fill
! be careful about mf for ODE pack
!-----------------------------------------------------------------------
!
integer neq, j, ia, ja
double precision t, y, pdj
dimension y(neq), ia(*), ja(*), pdj(neq)

return
end

!do j = 1, ncell			! more complex structure
!    !c	! CONTINUITY OF FLUXES on wall gas/polymer - expressed by henry s law on side of gas
!    cwg = (Dair*dcell*cAIR((j-1)*(onecell)+1) + D*hwall*cAIR((j-1)*(onecell)  ))/(D*hwall+Dair*dcell*HHair)
!    !c
!    !c	! flux through the gas and wall - constant flux
!    !c	! flux_intensity = D * (LEFT-RIGHT)
!    jAIR((j-1)*(onecell)) = D*(cAIR((j-1)*(onecell)) - cwg)/(dcell/2.0d0)
!
!    !c	diffusion through the wall
!    do i = 1, onecell-2	! fluxes through the wall
!        jAIR((j-1)*(onecell)+i) =                                 &
!                    Dair/hwall * (cAIR((j-1)*(onecell) + i)      &
!                                - cAIR((j-1)*(onecell) + i+1))
!!        write(*,*) (j-1)*(onecell)
!!        continue
!    enddo
!    !c	! CONTINUITY OF FLUXES on wall polymer/gas, GAS SIDE
!    cwg = (Dair *dcell*cAIR((j-1)*onecell+(onecell-1))        &
!        +  D*    hwall*cAIR((j-1)*onecell+(onecell  )))      &
!        /(D*hwall+Dair*dcell*HHair)
!    !c
!    jAIR((j-1)*onecell+(onecell-1))=Dair*(cAIR((j-1)*onecell+(onecell-1))  &
!                                 - cwg*HHair) /(hwall / 2.0d0)
!enddo

!jAIR(nFV) = 0.0d0		! last cell, zero flux
!continue

!c	balances in walls
!do j = 1, ncell
!    do i = 1, onecell-1
!        dcAIR((j-1)*onecell+i) = (jAIR((j-1)*onecell + i-1)         &
!                               -  jAIR((j-1)*onecell + i  ))/hwall
!    enddo
!enddo
!continue
!!c
!!c	ballance of the cells
!do j = 1, ncell
!    dcAIR((j-1)*onecell + onecell) = (jAIR((j-1)*onecell + onecell-1)      &
!                                 - jAIR((j-1)*onecell+onecell))/dcell
!enddo
