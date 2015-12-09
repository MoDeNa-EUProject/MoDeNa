      subroutine vandermonde(a,mc,n1,cc,ipiv,info)
      
      integer n1,cc,info,ipiv(n1)
      double precision a(n1,n1), mc(cc,n1)
      double precision mmc(n1,cc)
      
      a=transpose(a)
      mmc=transpose(mc)
      

      call dgesv(n1,cc,a,n1,ipiv,mmc,n1,info)
      
      mc = transpose(mmc)
            
      end subroutine
