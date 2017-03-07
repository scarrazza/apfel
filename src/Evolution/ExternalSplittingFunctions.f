************************************************************************
*
*     ExternalSplittingFunctions.f:
*
*     This function returns the "external" splitting function matrix on the
*     interpolation grid. This is possible only if the internal grids are 
*     used or if one single external grid is provided.
*
************************************************************************
      subroutine ComputeExternalSplittingFunctions(Bs2Bs,pt,nf,x,beta)
*
      implicit none
*
      include "../commons/EvolOp.h"
      include "../commons/grid.h"
      include "../commons/ipt.h"
      include "../commons/integrals.h"
      include "../commons/transQCD.h"
**
*     Input Variables
*
      integer pt,nf
      integer beta
      double precision x
      character*5 Bs2Bs
**
*     Internal Variables
*
      integer k,l,g,h
      integer n
      integer alpha
      integer lx,ux
      double precision w_int_gen,wg
      double precision tol
      parameter(tol=1d-10)

      double precision sfEv2Ev(0:13,0:13)
      double precision sfEv2Ph(-7:6,0:13)
      double precision sfPh2Ph(-7:6,-7:6)
      common / ExtSplittingFuncsAPFEL / sfEv2Ev,sfEv2Ph,sfPh2Ph

      character*5 Bs2BsC
      common / ExtSplittingFuncsBasisAPFEL / Bs2BsC
*
*     Check whether the evolution operator has been actually computed
*
      if(.not.EvolOp)then
         write(6,*) "The evolution operator computation is disabled."
         write(6,*) "The 'ExternalSplittingFunctions' function cannot",
     1              " be used."
         write(6,*) "   "
         call exit(-10)
      endif
*
*     Check consistency of the input variables
*
      if(Bs2Bs.ne."Ev2Ev".and.
     1   Bs2Bs.ne."Ev2Ph".and.
     2   Bs2Bs.ne."Ph2Ph")then
         write(6,*) "In ExternalSplittingFunctions.f:"
         write(6,*) "Invalid Basis flag, Bs2Bs = ",Bs2Bs
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'Ev2Ev'"
         write(6,*) "- 'Ev2Ph'"
         write(6,*) "- 'Ph2Ph'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(pt.lt.0.or.pt.gt.ipt)then
         write(6,*) "In ExternalSplittingFunctions.f:"
         write(6,*) "Perturbative order out of range, pt = ",pt
         write(6,*) "  "
         write(6,'(a,i2,a)') " pt must be in the range [0,",ipt,"]"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
         write(6,*) "In ExternalSplittingFunctions.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
      if (x.lt.xmin(1)) x = xmin(1)
      if (x.gt.xmax) x = 1d0
      if(beta.lt.0.or.beta.gt.nin(0))then
         write(6,*) "In ExternalSplittingFunctions.f:"
         write(6,*) "Invalid index, beta =",beta
         call exit(-10)
      endif
*
*     Construct splitting function matrix in the Ev2Ev bases
*
      Bs2BsC = Bs2Bs
      do k=0,13
         do l=0,13
            sfEv2Ev(k,l) = 0d0
         enddo
      enddo
*
      n = inter_degree(0)
      do alpha=0,beta
         if(xg(0,alpha).gt.x) goto 11
      enddo
 11   lx = alpha - 1
      ux = lx + n + 1
      do alpha=lx,ux
         wg = w_int_gen(n,alpha,x)
         if(wg.eq.0d0) cycle
*
         sfEv2Ev(1,1)   = sfEv2Ev(1,1)   + SP(0,nf,4,pt,alpha,beta) * wg
         sfEv2Ev(1,2)   = sfEv2Ev(1,2)   + SP(0,nf,5,pt,alpha,beta) * wg
         sfEv2Ev(2,1)   = sfEv2Ev(2,1)   + SP(0,nf,6,pt,alpha,beta) * wg
         sfEv2Ev(2,2)   = sfEv2Ev(2,2)   + SP(0,nf,7,pt,alpha,beta) * wg
         sfEv2Ev(3,3)   = sfEv2Ev(3,3)   + SP(0,nf,3,pt,alpha,beta) * wg
         sfEv2Ev(4,4)   = sfEv2Ev(4,4)   + SP(0,nf,2,pt,alpha,beta) * wg
         sfEv2Ev(5,5)   = sfEv2Ev(5,5)   + SP(0,nf,2,pt,alpha,beta) * wg
         sfEv2Ev(9,9)   = sfEv2Ev(9,9)   + SP(0,nf,1,pt,alpha,beta) * wg
         sfEv2Ev(10,10) = sfEv2Ev(10,10) + SP(0,nf,1,pt,alpha,beta) * wg
         if(nf.eq.3)then
            sfEv2Ev(6,3)   = sfEv2Ev(6,3)
     1                     + SP(0,nf,3,pt,alpha,beta) * wg
            sfEv2Ev(7,3)   = sfEv2Ev(7,3)
     1                     + SP(0,nf,3,pt,alpha,beta) * wg
            sfEv2Ev(8,3)   = sfEv2Ev(8,3)
     1                     + SP(0,nf,3,pt,alpha,beta) * wg
            sfEv2Ev(11,1)  = sfEv2Ev(11,1)
     1                     + SP(0,nf,4,pt,alpha,beta) * wg
            sfEv2Ev(11,2)  = sfEv2Ev(11,2)
     1                     + SP(0,nf,5,pt,alpha,beta) * wg
            sfEv2Ev(12,1)  = sfEv2Ev(12,1)
     1                     + SP(0,nf,4,pt,alpha,beta) * wg
            sfEv2Ev(12,2)  = sfEv2Ev(12,2)
     1                     + SP(0,nf,5,pt,alpha,beta) * wg
            sfEv2Ev(13,1)  = sfEv2Ev(13,1)
     1                     + SP(0,nf,4,pt,alpha,beta) * wg
            sfEv2Ev(13,2)  = sfEv2Ev(13,2)
     1                     + SP(0,nf,5,pt,alpha,beta) * wg
         elseif(nf.eq.4)then
            sfEv2Ev(6,6)   = sfEv2Ev(6,6)
     1                     + SP(0,nf,2,pt,alpha,beta) * wg
            sfEv2Ev(7,3)   = sfEv2Ev(7,3)
     1                     + SP(0,nf,3,pt,alpha,beta) * wg
            sfEv2Ev(8,3)   = sfEv2Ev(8,3)
     1                     + SP(0,nf,3,pt,alpha,beta) * wg
            sfEv2Ev(11,11) = sfEv2Ev(11,11)
     1                     + SP(0,nf,1,pt,alpha,beta) * wg
            sfEv2Ev(12,1)  = sfEv2Ev(12,1)
     1                     + SP(0,nf,4,pt,alpha,beta) * wg
            sfEv2Ev(12,2)  = sfEv2Ev(12,2)
     1                     + SP(0,nf,5,pt,alpha,beta) * wg
            sfEv2Ev(13,1)  = sfEv2Ev(13,1)
     1                     + SP(0,nf,4,pt,alpha,beta) * wg
            sfEv2Ev(13,2)  = sfEv2Ev(13,2)
     1                     + SP(0,nf,5,pt,alpha,beta) * wg
         elseif(nf.eq.5)then
            sfEv2Ev(6,6)   = sfEv2Ev(6,6)
     1                     + SP(0,nf,2,pt,alpha,beta) * wg
            sfEv2Ev(7,7)   = sfEv2Ev(7,7)
     1                     + SP(0,nf,2,pt,alpha,beta) * wg
            sfEv2Ev(8,3)   = sfEv2Ev(8,3)
     1                     + SP(0,nf,3,pt,alpha,beta) * wg
            sfEv2Ev(11,11) = sfEv2Ev(11,11)
     1                     + SP(0,nf,1,pt,alpha,beta) * wg
            sfEv2Ev(12,12) = sfEv2Ev(12,12)
     1                     + SP(0,nf,1,pt,alpha,beta) * wg
            sfEv2Ev(13,1)  = sfEv2Ev(13,1)
     1                     + SP(0,nf,4,pt,alpha,beta) * wg
            sfEv2Ev(13,2)  = sfEv2Ev(13,2)
     1                     + SP(0,nf,5,pt,alpha,beta) * wg
         elseif(nf.eq.6)then
            sfEv2Ev(6,6)   = sfEv2Ev(6,6)
     1                     + SP(0,nf,2,pt,alpha,beta) * wg
            sfEv2Ev(7,7)   = sfEv2Ev(7,7)
     1                     + SP(0,nf,2,pt,alpha,beta) * wg
            sfEv2Ev(8,8)   = sfEv2Ev(8,8)
     1                     + SP(0,nf,2,pt,alpha,beta) * wg
            sfEv2Ev(11,11) = sfEv2Ev(11,11)
     1                     + SP(0,nf,1,pt,alpha,beta) * wg
            sfEv2Ev(12,12) = sfEv2Ev(12,12)
     1                     + SP(0,nf,1,pt,alpha,beta) * wg
            sfEv2Ev(13,13) = sfEv2Ev(13,13)
     1                     + SP(0,nf,1,pt,alpha,beta) * wg
         endif
      enddo
*
*     Rotate if needed
*
      if(Bs2Bs.eq."Ev2Ph")then
         do k=0,13
            do l=0,13
               sfEv2Ph(k-7,l) = 0d0
               do h=1,13
                  sfEv2Ph(k-7,l) = sfEv2Ph(k-7,l)
     1                           + Tev2phQCD(nf,k,h) * sfEv2Ev(h,l)
               enddo
            enddo
         enddo
      elseif(Bs2Bs.eq."Ph2Ph")then
         do k=0,13
            do l=0,13
               sfPh2Ph(k-7,l-7) = 0d0
               do h=1,13
                  do g=1,13
                     sfPh2Ph(k-7,l-7) = sfPh2Ph(k-7,l-7)
     1                                + Tev2phQCD(nf,k,h)
     2                                * sfEv2Ev(h,g)
     3                                * Tph2evQCD(nf,g,l)
                  enddo
               enddo
            enddo
         enddo
      endif
*
      return
      end
*
************************************************************************
      function ExternalSplittingFunctions(i,j)
*
      implicit none
**
*     Input Variables
*
      integer i,j
**
*     Internal Variables
*
      double precision sfEv2Ev(0:13,0:13)
      double precision sfEv2Ph(-7:6,0:13)
      double precision sfPh2Ph(-7:6,-7:6)
      common / ExtSplittingFuncsAPFEL / sfEv2Ev,sfEv2Ph,sfPh2Ph
*
      character*5 Bs2BsC
      common / ExtSplittingFuncsBasisAPFEL / Bs2BsC
**
*     Output Variables
*
      double precision ExternalSplittingFunctions
*
      if(Bs2BsC.eq."Ev2Ev")then
         ExternalSplittingFunctions = sfEv2Ev(i,j)
      elseif(Bs2BsC.eq."Ev2Ph")then
         ExternalSplittingFunctions = sfEv2Ph(i,j)
      elseif(Bs2BsC.eq."Ph2Ph")then
         ExternalSplittingFunctions = sfPh2Ph(i,j)
      endif
*
      return
      end
