************************************************************************
*
*     xPDFxQ.f:
*
*     This function returns the value of the i-th PDF in the physical
*     basis at the scale Q in GeV and for the bjorken variable x using
*     the interpolation of the cached PDFs.
*
************************************************************************
      function xPDFxQ(id,x,Q)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/CacheParams.h"
      include "../commons/scales.h"
      include "../commons/fph.h"
**
*     Input Variables
*
      integer id
      double precision x,Q
**
*     Internal Variables
*
      integer i
      integer idx,idQ
      integer alpha,tau
      integer isg,iQtrans
      double precision Q2
      double precision w_int_xQ
      double precision tol
      parameter(tol=1d-10)
**
*     Output Variables
*
      double precision xPDFxQ
*
*     Check that the APFEL evolution has been initialized
*
      if(InCachePDFs.ne."done")then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Impossible to use this function because PDFs"
         write(6,*) "have not been cached."
         write(6,*) "Call 'CachePDFsAPFEL(Q0)' before calling xPDFxQ."
         write(6,*) "   "
         call exit(-10)
      endif
*
      Q2 = Q * Q
*
*     Check consistency of the input variables
*
      if(abs(id).le.6)then
         i = id
      elseif(id.eq.21)then
         i = 0
      elseif(abs(id).ne.11.or.abs(id).ne.13.or.abs(id).ne.15)then
         i = ( id - sign(1,id) * 9 ) / 2
      elseif(id.eq.22)then
         i = id
      else
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Invalid PDF index, id =",id
         call exit(-10)
      endif
*
*     Check that PDFs have been called in the allowed region
*
      if(x.lt.xmin(1)*(1d0-tol).or.x.gt.xmax*(1d0+tol))then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Value of x out of range, x =",x
         call exit(-10)
      endif
*
      if(Q2.lt.Q2min*(1d0-tol).or.Q2.gt.Q2max*(1d0+tol))then
         write(6,*) "In xPDFxQ.f:"
         write(6,*) "Value of Q out of range, Q =",Q
         call exit(-10)
      endif
*
*     trim away numerical fluctiations
*
      if(x.lt.xmin(1)) x  = xmin(1)
      if(x.gt.xmax)    x  = 1d0
*
      if(Q2.lt.Q2min)  Q2 = Q2min
      if(Q2.gt.Q2max)  Q2 = Q2max
*
*     If the value of Q2 happens to fall in the tiny gap
*     between two subsequent subgrids, push it on the first
*     node of the second grid.
*
      iQtrans = nQ(nfin)
      do isg=nfin,nffi-1
         if(Q2.gt.Q2g(iQtrans).and.Q2.lt.Q2g(iQtrans+1))then
            Q2 = Q2g(iQtrans+1)
            goto 101
         endif
         iQtrans = iQtrans + nQ(isg+1)
      enddo
*
*     Interpolation
*
 101  xPDFxQ = 0d0
      idx = inter_degree(0)
      idQ = inter_degreeQ
      if(abs(id).le.6)then
         do tau=1,nQ2g
            do alpha=0,nin(0)
               xPDFxQ = xPDFxQ + w_int_xQ(idx,idQ,alpha,tau,x,Q2)
     1                * fphxQ(i,alpha,tau)
            enddo
         enddo
      elseif(abs(id).eq.11.or.abs(id).eq.13.or.abs(id).eq.15)then
         do tau=1,nQ2g
            do alpha=0,nin(0)
               xPDFxQ = xPDFxQ + w_int_xQ(idx,idQ,alpha,tau,x,Q2)
     1                * fleptonxQ(i,alpha,tau)
            enddo
         enddo
      elseif(id.eq.22)then
         do tau=1,nQ2g
            do alpha=0,nin(0)
               xPDFxQ = xPDFxQ + w_int_xQ(idx,idQ,alpha,tau,x,Q2)
     1                * fgammaxQ(alpha,tau)
            enddo
         enddo
      endif
      if(dabs(xPDFxQ).le.1d-12) xPDFxQ = 0d0
*
      return
      end
*
************************************************************************
      function w_int_xQ(kx,kQ,beta,tau,x,Q2)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/CacheParams.h"
**
*     Input Variables
*
      integer kx,kQ,beta,tau
      double precision x,Q2
**
*     Internal Variables
*
      integer j
      integer delta,xbound,Qbound
**
*     Output Variables
*
      double precision w_int_xQ
*
      w_int_xQ = 0d0
      xbound = beta - kx
      Qbound = tau  - kQ
*
      if(kx.gt.beta) xbound = 0
      if(kQ.gt.tau)  Qbound = 0
*
      if(x.lt.xg(0,xbound).or.x.ge.xg(0,beta+1)) return
      if(Q2.lt.Q2g(Qbound).or.Q2.ge.Q2g(tau+1)) return
*
*     x-grid
*
      do j=0,beta-xbound
         if(x.ge.xg(0,beta-j).and.x.lt.xg(0,beta-j+1))then
            w_int_xQ = 1d0
            do delta=0,kx
               if(delta.ne.j) w_int_xQ = w_int_xQ
     1              * dlog(x/xg(0,beta-j+delta)) 
     2              / dlog(xg(0,beta)/xg(0,beta-j+delta))
            enddo
         endif
      enddo
*
*     Q2 grid
*
      do j=0,tau-Qbound
         if(Q2.ge.Q2g(tau-j).and.Q2.lt.Q2g(tau-j+1))then
            do delta=0,kQ
               if(delta.ne.j) w_int_xQ = w_int_xQ
     1              * dlog(Q2/Q2g(tau-j+delta)) 
     2              / dlog(Q2g(tau)/Q2g(tau-j+delta))
            enddo
         endif
      enddo
*
      return
      end
