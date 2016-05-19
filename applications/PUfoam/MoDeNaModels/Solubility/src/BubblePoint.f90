
!> \file BubblePoint.f90
!! \brief This Subroutine performs a bubble point calculation.



Module BubblePoint

 Implicit None

 Private
 Public :: BubblePointCalculation

 Contains


 !> This subroutine finds the equilibrium pressure and vapor composition
 !> for a given liquid composition and fixed temperature.
 !> By setting the variable iterate_t = 1, pressure can be fixed and
 !> temperature will be an output.

 Subroutine BubblePointCalculation

  Use BASIC_VARIABLES
  Use STARTING_VALUES, Only: bubble_point_rachford_rice
  Use utilities, Only: file_open

  !local variables
  Integer :: converg,iterate_t
  Real    :: rhoi1(nc),rhoi2(nc)
  Real    :: H
  CHARACTER (LEN=50)                     :: filename


  Call READ_INPUT
  nphas = 2
  iterate_t = 0 !iterate pressure to reach the given liquid composition
  !initial values of liaquid and gas packing fraction
  dense(1) = 0.4
  dense(2) = 1.E-05
  !set initial mole fractions for gas and liquid phase
  xi(1,1:ncomp) = xiF(1:ncomp)  !this is the fixed liquid composition
  xi(2,1)       = 0.999
  xi(2,2:ncomp) = (1. - xi(2,1)) / REAL(ncomp-1)

  !rhoi1 and rhoi2 are output values
  call bubble_point_rachford_rice(iterate_t,converg,rhoi1,rhoi2)

!   write(*,*)'converg p', converg, p
!   write(*,*)'xl',xi(1,1:ncomp)
!   write(*,*)'xv',xi(2,1:ncomp)

  !Calculate solubility coefficient H of component 1
  !according to x1*H = y1*p
  H = xi(2,1) * p / xi(1,1)
  filename = './out.txt'
  CALL file_open(filename,78)
  ! write(78,*) H / 100000.
  ! write(*,*)'Solubility Coefficient /bar', H/100000.0
  if (ncomp == 2) then
      H=xi(1,1)*mm(1)/(xi(1,2)*mm(2))/p*1e5
      write(*,*)'Solubility Coefficient g/g/bar', H
  elseif (ncomp == 3) then
      H=xi(1,1)*mm(1)/(xi(1,2)*mm(2)+xi(1,3)*mm(3))/p
      H=H*1100/mm(1)/1e-3
      write(*,*)'Solubility Coefficient mol/m3/Pa', H,mm(1)
  endif
  write(78,*) H


 End Subroutine BubblePointCalculation


End Module BubblePoint
