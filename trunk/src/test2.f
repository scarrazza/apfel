************************************************************************
*
*      test2.f:
*
*      Program to test the stucture of the matrices of integrals needed
*      to construct the evolution matrices.
*
************************************************************************
      program test2
*
      implicit none
*
      include "../commons/scales.h"
      include "../commons/grid.h"
      include "../commons/Nf_FF.h"
      include "../commons/wrap.h"
*
      integer alpha,beta,kk
      double precision mu2
      double precision integralsQCD
      double precision integralsMatching
      double precision t1,t2
      double precision coup,a_QCD
      real*4 output(0:nint_max,0:nint_max)
*
      call SetNumberOfGrids(1)
      call SetGridParameters(1,20,3,1d-5)
*
*     Initialize
*
      call InitializeAPFEL
*
      mu2 = 100d0
      coup = a_QCD(mu2)
*
      call cpu_time(t1)
      igrid = 1
c      do kk=1,7
      do kk=7,7
         do alpha=0,nin(1)
            do beta=0,nin(1)
               output(alpha,beta) = integralsQCD(alpha,beta,coup,kk)
            enddo
            write(6,*) (output(alpha,beta),beta=0,nin(1))
         enddo
         write(6,*) "  "
      enddo
      call cpu_time(t2)
*
      write(6,*) "Time elapsed =",t2-t1," s"
      write(6,*) "  "
      stop
*
      call cpu_time(t1)
      do kk=1,5
         do alpha=0,nin(1)
            do beta=0,nin(1)
               call RSLintegralsMatching(alpha,beta)
               output(alpha,beta)=integralsMatching(alpha,beta,coup,kk)
            enddo
            write(6,*) (output(alpha,beta),beta=0,nin(1))
         enddo
         write(6,*) "  "
      enddo
      call cpu_time(t2)
*
      write(6,*) "Time elapsed =",t2-t1," s"
      write(6,*) "  "
*
      end
