!> @file
!! subroutines for special mathematical functions
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module specfun
    use constants, only: dp,pi,iu,eulergamma
    implicit none
    private
    public factorial,digamma,bessel,hankel
contains
!********************************BEGINNING*************************************
!> calculates Bessel functions of the first and second kind and
!! their derivatives for complex argument
subroutine bessel(n,x,tol,j,dj,y,dy)
    use besselj, only: ZBESJ
    use bessely, only: ZBESY
    integer, intent(in) :: n  !maximum order calculated
    complex(dp), intent(in) :: x  !argument
    real(dp), intent(in) :: tol  !relative tolerance, currently not used
    complex(dp), dimension(0:n), intent(out) :: &
        j,&    !bessel function of the first kind
        dj,&    !derivative of bessel function of the first kind
        y,&    !bessel function of the second kind
        dy    !derivative of bessel function of the second kind
    real(dp), dimension(n+1) :: cyr,cyi,cwr,cwi
    integer :: i,m,nz,ierr
    call ZBESJ(realpart(x),imagpart(x),0.e0_dp,1,n+1,cyr,cyi,nz,ierr)
    if (nz>0) then
        write(*,*) 'underflow encountered'
        write(*,*) 'bessel function J not calculated',n,x
        stop
    endif
    if (ierr /= 0) then
        write(*,*) 'did not converge'
        write(*,*) 'bessel function J not calculated',n,x
        stop
    endif
    do i=0,n
        j(i)=cmplx(cyr(i+1),cyi(i+1),kind=dp)
    enddo
!    write(*,*) x
!    write(*,*) j
    call ZBESY(realpart(x),imagpart(x),0.e0_dp,1,n+1,cyr,cyi,nz,cwr,cwi,ierr)
    if (nz>0) then
        write(*,*) 'underflow encountered'
        write(*,*) 'bessel function Y not calculated',n,x
        stop
    endif
    if (ierr /= 0) then
        write(*,*) 'did not converge'
        write(*,*) 'bessel function Y not calculated',n,x
        stop
    endif
    do i=0,n
        y(i)=cmplx(cyr(i+1),cyi(i+1),kind=dp)
    enddo
!    write(*,*) x
!    write(*,*) y
!    stop
    !http://keisan.casio.com/exec/system/1180573474
    dj(0)=-j(1)
    dy(0)=-y(1)
    do m=1,n
        dj(m)=j(m-1)-m/x*j(m)
        dy(m)=y(m-1)-m/x*y(m)
    enddo
end subroutine bessel
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates Hankel functions of the first and second kind and
!! their derivatives for complex argument
subroutine hankel(n,x,tol,h1,dh1,h2,dh2)
    integer, intent(in) :: n  !maximum order calculated
    complex(dp), intent(in) :: x  !argument
    real(dp), intent(in) :: tol  !relative tolerance
    complex(dp), dimension(0:n), intent(out) :: &
        h1,&    !hankel function of the first kind
        dh1,&    !derivative of hankel function of the first kind
        h2,&    !hankel function of the second kind
        dh2    !derivative of hankel function of the second kind
    integer :: i
    complex(dp), dimension(0:n) :: j,dj,y,dy
    !http://keisan.casio.com/has10/SpecExec.cgi?id=system/2006/1222514923
    call bessel(n,x,tol,j,dj,y,dy)
    h1=j+iu*y
    h2=j-iu*y
    dh1(0)=-h1(1)
    dh2(0)=-h2(1)
    do i=1,n
        dh1(i)=h1(i-1)-i/x*h1(i)
        dh2(i)=h2(i-1)-i/x*h2(i)
    enddo
end subroutine hankel
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates factorial
recursive function factorial(n) result(fact)
    integer, intent(in) :: n  !argument
    real(dp) :: fact
    if (n<0) then
        stop 'factorial defined only for nonnegative integers'
    elseif (n==0) then
        fact=1
    else
        fact=n*factorial(n-1)
    endif
end function factorial
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculates digamma function
!! taken from http://mathworld.wolfram.com/DigammaFunction.html
real(dp) function digamma(n)
    integer, intent(in) :: n  !argument
    integer :: i
    digamma=-eulergamma
    do i=1,n-1
        digamma=digamma+1.0_dp/i
    enddo
end function digamma
!***********************************END****************************************


!SUBROUTINE JNCOMP(Z,NS,JN)
!    !taken from Daniel Mackowski code for scattering of radiation by long cylinders
!    !it might be more effective, I did not test it thoroughly
!    IMPLICIT none
!    COMPLEX*16 JN(0:*),A,Z
!    integer :: n,nd,ns
!
!!  CALCULATES THE INTEGER-ORDER BESSEL FUNCTION JN(Z)
!!  N=0,1,...NS, FOR COMPLEX ARGUMENT Z.  REFER TO BOHREN & HUFFMAN
!!  FOR DETAILS.  CODED 9/20/90
!
!    ND=NINT((101+CDABS(Z))**.499+NS)
!    ND=2*NINT(ND/2.)
!    JN(ND)=0.
!    JN(ND-1)=1.D-32
!    A=0.
!    DO N=ND-1,3,-2
!        JN(N-1)=2.*N*JN(N)/Z-JN(N+1)
!        JN(N-2)=2.*(N-1)*JN(N-1)/Z-JN(N)
!        A=A+JN(N-1)
!    enddo
!    JN(0)=2.*JN(1)/Z-JN(2)
!    A=JN(0)+2.*A
!    DO N=0,NS
!        JN(N)=JN(N)/A
!    enddo
!RETURN
!END subroutine


!!********************************BEGINNING*************************************
!!calculates Bessel functions of the first and second kind and their derivatives for complex argument
!!doesn't work properly for high orders
!subroutine bessel(n,x,tol,j,dj,y,dy)
!    integer, intent(in) :: n  !maximum order calculated
!    complex(dp), intent(in) :: x  !argument
!    real(dp), intent(in) :: tol  !relative tolerance
!    complex(dp), dimension(0:n), intent(out) :: j    !bessel function of the first kind
!    complex(dp), dimension(0:n), intent(out) :: dj    !derivative of bessel function of the first kind
!    complex(dp), dimension(0:n), intent(out) :: y    !bessel function of the second kind
!    complex(dp), dimension(0:n), intent(out) :: dy    !derivative of bessel function of the second kind
!    integer :: l,m,lmax=200
!    complex(dp) :: jold,yold
!    jold=0;yold=0
!    j=0;dj=0;y=0;dy=0
!!    write(*,*) x
!    do m=0,n
!        do l=0,lmax
!!            write(*,*) l,m,j(m)
!!            j(m)=j(m)+(-1)**l/(2**(2*l+m)*factorial(l)*factorial(m+l)+0.0_dp)*x**(2*l+m) !http://mathworld.wolfram.com/BesselFunctionoftheFirstKind.html (it works just as well)
!            j(m)=j(m)+(-1)**l/(factorial(l)*factorial(m+l)+0.0_dp)*(x/2)**(2*l+m) !http://keisan.casio.com/exec/system/1180573474
!            if (abs(j(m)-jold)/abs(j(m))<tol) then
!!                write(*,*) m,l
!                exit
!            elseif (l==lmax) then
!                write(*,*) 'bessel function J not calculated',n,x
!                stop
!            endif
!            jold=j(m)
!        enddo
!        write(*,*) j(m)
!        !http://mathworld.wolfram.com/BesselFunctionoftheSecondKind.html
!        y(m)=2/pi*log(x/2)*j(m) !second term
!        do l=0,m-1
!            y(m)=y(m)-(x/2)**(-m)/pi*factorial(m-l-1)/factorial(l)*(x**2/4)**l !first term
!        enddo
!        do l=0,lmax
!!            write(*,*) l,m,y(m)
!            y(m)=y(m)-(x/2)**m/pi*(digamma(l+1)+digamma(m+l+1))*(-x**2/4)**l/(factorial(l)*factorial(m+l)) !third term
!!            write(*,*) l
!            if (abs(y(m)-yold)/abs(y(m))<tol) then
!!                write(*,*) l
!                exit
!            elseif (l==lmax) then
!                write(*,*) 'bessel function Y not calculated',n,x
!                stop
!            endif
!            yold=y(m)
!        enddo
!    enddo
!    !http://keisan.casio.com/exec/system/1180573474
!    dj(0)=-j(1)
!    dy(0)=-y(1)
!    do m=1,n
!        dj(m)=j(m-1)-m/x*j(m)
!        dy(m)=y(m-1)-m/x*y(m)
!    enddo
!end subroutine bessel
!!***********************************END****************************************
end module specfun
